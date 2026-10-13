#!/usr/bin/env Rscript
# ============================================================================
# CMIP6 kriging downscaling - multi-covariate version
# ============================================================================
# Built for multiple variable-specific covariates.
# Example: relative humidity using DEM + temperature + vapor pressure.
# Uses month-specific covariate files and a custom kriging equation.
# Two-step downscaling: 28km -> 9km -> 1km.
#
# Engine: automap::autoKrige (universal kriging), the same engine as the
# single-covariate templates. This template does not use KrigR, so no KrigR
# source files (Kriging.R, misc.R, and so on) are needed.
# ============================================================================

# Suppress warnings
options(warn = -1)

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

args <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=", "", args[grep("--file=", args)]))
if (length(script_path) > 0) {
    setwd(script_path)
}

# ============================================================================
# CONFIGURATION (Replaced by Python script)
# ============================================================================

input_file <- "{{INPUT_FILE}}"
output_file <- "{{OUTPUT_FILE}}"
se_file <- "{{SE_FILE}}"
variable_name <- "{{VARIABLE}}"
target_date <- "{{DATE}}"
target_res <- as.numeric("{{TARGET_RESOLUTION}}")
timestep_idx <- "{{TIMESTEP_INDEX}}"
nmax <- as.numeric("{{NMAX}}")
cores <- as.numeric("{{CORES}}")
iteration <- as.numeric("{{ITERATION}}")  # 1 or 2
month_number <- as.numeric("{{MONTH_NUMBER}}")  # 1-12

# ----------------------------------------------------------------------------
# Coverage gate parameters (enforced by the Python driver; reported here)
# ----------------------------------------------------------------------------
# min_coverage: minimum fraction of the area that must carry data for the
#   downscale to pass. 0.9988 area covered == 0.0012 (0.12%) NaN allowed.
# coverage_denominator: what "the area" means.
#   "domain" = cells expected to have data (all covariates finite) -> measures
#              kriging holes/failures; ocean/out-of-domain cells are excluded.
#   "grid"   = every cell in the target grid (whole bounding box).
# BOTH are computed and printed. The driver decides pass/fail on the selected
# one. Fallbacks apply if the placeholders are left unsubstituted.
min_coverage <- suppressWarnings(as.numeric("{{MIN_COVERAGE}}"))
if (is.na(min_coverage)) min_coverage <- 0.9988
coverage_denominator <- "{{COVERAGE_DENOMINATOR}}"
if (!(coverage_denominator %in% c("domain", "grid"))) coverage_denominator <- "domain"

# ----------------------------------------------------------------------------
# Result-scan bounds (physical sanity check on each kriging attempt).
# A result is "meaningful" only if almost all predicted values fall within the
# physical range for the variable. For hurs (relative humidity) that is 0-100 %.
# valid_tol_frac widens the band to allow small kriging overshoot before the
# result is rejected; an accepted field is finally clamped to [valid_min, valid_max].
valid_min <- suppressWarnings(as.numeric("{{VALID_MIN}}"))
valid_max <- suppressWarnings(as.numeric("{{VALID_MAX}}"))
if (is.na(valid_min)) valid_min <- 0
if (is.na(valid_max)) valid_max <- 100
valid_tol_frac <- 0.10                     # 10% of the range as overshoot tolerance

# GMTED and covariates
gmted_file <- "{{GMTED_FILE}}"  # GMTED covariate file for this iteration
kriging_equation <- "{{KRIGING_EQUATION}}"  # e.g., "Data ~ DEM + layer.1 + layer.2"

# Multiple covariate paths (JSON array as string)
covariate_configs <- "{{COVARIATE_CONFIGS}}"  # Will be parsed

# Spatial extent
extent_lonmin <- as.numeric("{{EXTENT_LONMIN}}")
extent_lonmax <- as.numeric("{{EXTENT_LONMAX}}")
extent_latmin <- as.numeric("{{EXTENT_LATMIN}}")
extent_latmax <- as.numeric("{{EXTENT_LATMAX}}")

# Shapefile (optional)
shapefile_path <- "{{SHAPEFILE_PATH}}"

# Working directory
working_dir <- "{{WORKING_DIR}}"
setwd(working_dir)

# ============================================================================
# LOAD LIBRARIES
# ============================================================================

cat("======================================================================\n")
cat("CMIP6 KRIGING DOWNSCALING - MULTI-COVARIATE VERSION\n")
cat("======================================================================\n")
cat("Variable:", variable_name, "\n")
cat("Date:", target_date, "\n")
cat("Iteration:", iteration, "\n")
cat("Month:", month_number, "\n")
cat("Target Resolution:", target_res, "degrees\n")
cat("Input:", input_file, "\n")
cat("Output:", output_file, "\n")
cat("======================================================================\n\n")

suppressPackageStartupMessages({
    library(raster)
    library(ncdf4)
    library(sp)
    library(rgdal)
    library(automap)
    library(gstat)
    library(foreach)
    library(doSNOW)
    library(stringr)
    library(jsonlite)
})

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

log_info <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

get_utm_zone <- function(lon) {
    zone <- floor((lon + 180) / 6) + 1
    list(zone  = zone,
         proj4 = sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", zone))
}

# Fatal error: print to BOTH stdout (visible on interactive runs) and stderr
# (captured by the Python driver, which logs result.stderr), then exit non-zero.
fail <- function(...) {
    msg <- paste0(...)
    cat("ERROR:", msg, "\n")
    message("ERROR: ", msg)
    quit(status = 1)
}

# ============================================================================
# STEP 0: Kriging engine
# ============================================================================
# No KrigR sourcing. Universal kriging is done with automap::autoKrige below.
cat("Step 0: Kriging engine = automap::autoKrige (no KrigR sourcing required)\n")

# ============================================================================
# STEP 1: Setup spatial extent
# ============================================================================

cat("\nStep 1: Setting up spatial extent...\n")

# Create extent object
cat("  Using Texas extent:\n")
cat("    Longitude:", extent_lonmin, "to", extent_lonmax, "\n")
cat("    Latitude:", extent_latmin, "to", extent_latmax, "\n")

e <- extent(c(extent_lonmin, extent_lonmax, extent_latmin, extent_latmax))
Extent <- as(e, 'SpatialPolygons')
crs(Extent) <- "+proj=longlat +datum=WGS84 +no_defs"

# Load and crop shapefile if provided
Extent_crop <- NULL
if (shapefile_path != "NULL" && file.exists(shapefile_path)) {
    cat("  Loading shapefile:", shapefile_path, "\n")
    Shape_shp <- readOGR(shapefile_path, verbose = FALSE)
    Extent_crop <- crop(x = Shape_shp, y = Extent)
    cat("  Shapefile loaded and cropped successfully\n")
} else {
    Extent_crop <- Extent
}

# ============================================================================
# STEP 2: Load GMTED covariate base
# ============================================================================

cat("\nStep 2: Loading GMTED covariate file...\n")
cat("  Loading:", gmted_file, "\n")

GMTED_DE <- readRDS(gmted_file)
cat("  Loaded GMTED with", length(GMTED_DE), "resolution levels\n")

# ============================================================================
# STEP 3: Load and prepare multiple month-specific covariates
# ============================================================================

cat("\nStep 3: Loading month-specific covariates...\n")
cat("  Month:", month_number, "\n")

# Parse covariate configurations from JSON (keep a list of named lists)
covariate_list <- fromJSON(covariate_configs, simplifyVector = FALSE)
cat("  Number of covariates to load:", length(covariate_list), "\n")

# Store covariate layers
covariate_layers_coarse <- list()
covariate_layers_fine <- list()

for (i in seq_along(covariate_list)) {
    cov_config <- covariate_list[[i]]
    cov_type <- cov_config$type
    cov_name <- cov_config$name

    cat("\n  Covariate", i, ":", cov_name, "(", cov_type, ")\n")

    if (cov_type == "monthly_file") {
        # Month-specific single-field file (e.g., vapor pressure climatology)
        file_pattern <- cov_config$path
        cov_file <- sprintf(file_pattern, month_number)

        # Resolve wildcards in the path (e.g., r* realization member)
        if (grepl("\\*", cov_file)) {
            matched_files <- Sys.glob(cov_file)
            if (length(matched_files) == 0) {
                fail("No covariate file matches pattern: ", cov_file)
            }
            if (length(matched_files) > 1) {
                cat("    WARNING: Multiple files match pattern; using first:\n")
                for (mf in matched_files) cat("      -", mf, "\n")
            }
            cov_file <- matched_files[1]
        }

        cat("    Loading:", cov_file, "\n")

        if (!file.exists(cov_file)) {
            fail("Covariate file not found: ", cov_file)
        }

        # Load and process
        vp <- mean(stack(cov_file))
        vp_cropped <- mean(crop(vp, Extent_crop))
        vp_masked <- mask(vp_cropped, Extent_crop)

        # Resample to GMTED resolutions
        cov_coarse <- resample(vp_masked, GMTED_DE[[1]], method = 'bilinear')
        cov_fine <- resample(vp_masked, GMTED_DE[[2]], method = 'bilinear')

        cat("    Processed and resampled\n")

    } else if (cov_type == "monthly_nc_stack") {
        # Pre-downscaled monthly temperature. The path now points at a
        # PRECOMPUTED monthly mean (cdo timmean), so this is a single field.
        # The mean() below is a no-op on a 1-layer file and stays correct if a
        # multi-layer file is ever supplied.
        file_pattern <- cov_config$path
        cov_file <- sprintf(file_pattern, month_number)

        # Resolve wildcards in the path (e.g., r* realization member)
        if (grepl("\\*", cov_file)) {
            matched_files <- Sys.glob(cov_file)
            if (length(matched_files) == 0) {
                fail("No covariate file matches pattern: ", cov_file)
            }
            if (length(matched_files) > 1) {
                cat("    WARNING: Multiple files match pattern; using first:\n")
                for (mf in matched_files) cat("      -", mf, "\n")
            }
            cov_file <- matched_files[1]
        }

        cat("    Loading:", cov_file, "\n")

        if (!file.exists(cov_file)) {
            fail("Covariate file not found: ", cov_file)
        }

        # Load and average (single-field mean file -> instant)
        tas_stack <- stack(cov_file)
        tas_mean <- mean(tas_stack, na.rm = TRUE)
        tas_masked <- mask(tas_mean, Extent_crop)

        # Resample to GMTED resolutions
        cov_coarse <- resample(tas_masked, GMTED_DE[[1]], method = 'bilinear')
        cov_fine <- resample(tas_masked, GMTED_DE[[2]], method = 'bilinear')

        cat("    Processed and resampled\n")

    } else {
        fail("Unknown covariate type: ", cov_type)
    }

    # Store in lists
    covariate_layers_coarse[[i]] <- cov_coarse
    covariate_layers_fine[[i]] <- cov_fine

    cat("    Resolution (coarse):", res(cov_coarse), "\n")
    cat("    Resolution (fine):", res(cov_fine), "\n")
}

# ============================================================================
# STEP 4: Stack covariates with GMTED
# ============================================================================

cat("\nStep 4: Stacking covariates...\n")

# Start with GMTED
cat("  Adding GMTED DEM...\n")
covs_1_list <- list(GMTED_DE[[1]])
covs_2_list <- list(GMTED_DE[[2]])

# Add all other covariates
for (i in seq_along(covariate_layers_coarse)) {
    cat("  Adding covariate", i, "...\n")
    covs_1_list[[i + 1]] <- covariate_layers_coarse[[i]]
    covs_2_list[[i + 1]] <- covariate_layers_fine[[i]]
}

# Create stacks
covs_1 <- stack(covs_1_list)
covs_2 <- stack(covs_2_list)

cat("  Coarse stack layers:", nlayers(covs_1), "\n")
cat("  Fine stack layers:", nlayers(covs_2), "\n")

# ============================================================================
# STEP 5: Load input climate data
# ============================================================================

cat("\nStep 5: Loading input data...\n")
cat("  Input file:", input_file, "\n")

# Load data using stack
data_start <- stack(input_file)
cat("  Loaded", nlayers(data_start), "layers\n")

# Crop and mask to extent
data_start <- crop(data_start, Extent_crop)
data_start <- mask(data_start, Extent_crop)
cat("  Cropped and masked to extent\n")
cat("  Final layers:", nlayers(data_start), "\n")

# ============================================================================
# STEP 6: Setup kriging equation
# ============================================================================

cat("\nStep 6: Setting up kriging equation...\n")
cat("  Equation:", kriging_equation, "\n")

# Parse kriging equation
KrigingEquation <- as.formula(kriging_equation)
cat("  Kriging equation ready\n")

# Verify equation matches number of covariates
cat("  Number of covariate layers:", nlayers(covs_1), "\n")
cat("  Expected: DEM (GMTED) + ", length(covariate_list), " additional covariates\n")

# ============================================================================
# STEP 7: Perform kriging (automap::autoKrige, universal kriging)
# ============================================================================

cat("\nStep 7: Performing kriging downscaling (adaptive cascade)...\n")
log_info("  Target resolution:", target_res, "degrees")
log_info("  NMAX:", nmax)

n_cov <- length(covariate_list)                 # extra covariates (tas, vp, ...)
predictor_names <- c("DEM", paste0("layer.", seq_len(n_cov)))
log_info("  Predictors:", paste(predictor_names, collapse = ", "))

# --- Training frame: source/coarse grid = GMTED_DE[[1]] ---
log_info("  Building training data at source resolution...")
data_coarse <- resample(data_start[[1]], GMTED_DE[[1]], method = "bilinear")

train_data <- as.data.frame(data_coarse, xy = TRUE)
names(train_data) <- c("x", "y", "Data")

covs_1_df <- as.data.frame(covs_1, xy = TRUE)
names(covs_1_df) <- c("x", "y", predictor_names)

training_df <- merge(train_data, covs_1_df, by = c("x", "y"))
training_df <- training_df[complete.cases(training_df), ]
log_info("  Training points:", nrow(training_df))
if (nrow(training_df) < 10) {
    fail("Too few training points after NA removal: ", nrow(training_df))
}

# --- Prediction grid: target/fine grid = GMTED_DE[[2]] ---
log_info("  Building prediction grid at target resolution...")
pred_df_full <- as.data.frame(covs_2, xy = TRUE)
names(pred_df_full) <- c("x", "y", predictor_names)
pred_df <- pred_df_full[complete.cases(pred_df_full[, predictor_names, drop = FALSE]), ]
log_info("  Prediction points:", nrow(pred_df))
if (nrow(pred_df) == 0) {
    fail("No prediction cells have complete covariates.")
}

# --- Covariate diagnostics (verbose, like the working templates) ---
for (p in predictor_names) {
    r <- suppressWarnings(cor(training_df$Data, training_df[[p]], use = "complete.obs"))
    log_info("  cor(Data,", p, "):", ifelse(is.na(r), "NA", round(r, 4)))
}

# --- Standardize response and predictors (Z-scores from training stats) ---
log_info("  Standardizing response and predictors...")
data_mean <- mean(training_df$Data, na.rm = TRUE)
data_sd   <- sd(training_df$Data, na.rm = TRUE)
if (is.na(data_sd) || data_sd == 0) fail("Response has zero variance; cannot krige.")
log_info("  Data mean:", round(data_mean, 4), "sd:", round(data_sd, 4))

pred_mean <- sapply(training_df[predictor_names], mean, na.rm = TRUE)
pred_sd   <- sapply(training_df[predictor_names], sd,   na.rm = TRUE)
pred_sd[is.na(pred_sd) | pred_sd == 0] <- 1     # guard a constant covariate

training_df$Data <- (training_df$Data - data_mean) / data_sd
for (p in predictor_names) {
    training_df[[p]] <- (training_df[[p]] - pred_mean[[p]]) / pred_sd[[p]]
    pred_df[[p]]     <- (pred_df[[p]]     - pred_mean[[p]]) / pred_sd[[p]]
}

# --- Reproject to UTM (metres) for a valid distance-based variogram ---
lon_center <- mean(c(extent_lonmin, extent_lonmax))
utm <- get_utm_zone(lon_center)
log_info("  UTM zone:", utm$zone, "(center lon", round(lon_center, 2), ")")

train_sp <- training_df
coordinates(train_sp) <- ~x + y
proj4string(train_sp) <- CRS("+proj=longlat +datum=WGS84")
train_sp <- spTransform(train_sp, CRS(utm$proj4))

pred_sp <- pred_df
coordinates(pred_sp) <- ~x + y
proj4string(pred_sp) <- CRS("+proj=longlat +datum=WGS84")
pred_sp <- spTransform(pred_sp, CRS(utm$proj4))

# --- Adaptive multi-covariate kriging cascade (always keeps the covariates) --
# autofitVariogram can fit a Gaussian model with a zero nugget, which makes the
# kriging system singular and yields all-NA predictions. The Gaussian model is
# kept; a small enforced nugget floor keeps the same model solvable.
#
# Order of attempts (each uses the covariates: Data ~ DEM + layer.1 + ...):
#   1. Full autofit (Sph/Exp/Gau/Mat/Ste), no floor, adaptive nmax.
#   2. If that fails or the scan rejects it: the same model set with an enforced
#      nugget floor (2% then 5% of sill), stepping nmax 75 -> 60 -> 40.
#   3. Otherwise, a covariate regression trend, lm(Data ~ DEM + layer.1 +
#      layer.2), predicted on the target grid.
# After each attempt the result is scanned for meaning (finite + physically in
# range + not degenerate); a rejected attempt is discarded and the next tried.
# All attempts are local (nmax); global kriging is not used.

# --- correlation-based adaptive parameters ---
cors <- sapply(predictor_names, function(p)
    suppressWarnings(cor(training_df$Data, training_df[[p]], use = "complete.obs")))
max_cor <- max(abs(cors), na.rm = TRUE)
if (!is.finite(max_cor)) max_cor <- 0
if (max_cor >= 0.50) {
    adaptive_nugget <- 0.01; nmax_ad <- nmax;                 cor_cat <- "Strong"
} else if (max_cor >= 0.30) {
    adaptive_nugget <- 0.02; nmax_ad <- max(20, round(nmax * 0.75)); cor_cat <- "Moderate"
} else {
    adaptive_nugget <- 0.05; nmax_ad <- max(15, round(nmax * 0.50)); cor_cat <- "Weak"
}
log_info(sprintf("  Adaptive: max|cor|=%.3f (%s) -> nugget=%.3f, adaptive nmax=%d",
                 max_cor, cor_cat, adaptive_nugget, nmax_ad))

nmax_steps <- unique(c(nmax, 60, 40))
nmax_steps <- nmax_steps[nmax_steps > 0 & nmax_steps <= nmax]   # e.g. 75,60,40

f_uk <- KrigingEquation

# Fit a variogram over the full model set. floor_frac > 0 raises the nugget
# (psill[1]) to at least floor_frac * sill, so the kriging system cannot be
# singular; floor_frac = 0 uses the variogram as autofit returned it.
fit_model <- function(formula, floor_frac) {
    models <- c("Sph", "Exp", "Gau", "Mat", "Ste")   # Gaussian always available
    vm <- tryCatch(autofitVariogram(formula, train_sp, model = models)$var_model,
                   error = function(e) NULL)
    if (!is.null(vm) && floor_frac > 0) {
        total_sill <- sum(vm$psill, na.rm = TRUE)
        floor_val  <- floor_frac * total_sill
        if (length(vm$psill) >= 1 && is.finite(vm$psill[1]) && vm$psill[1] < floor_val)
            vm$psill[1] <- floor_val
    }
    vm
}

# One kriging attempt -> standardized pred+var (+ the model used), or NULL if it
# errors / all-NA.
try_krige <- function(formula, nmax_val, floor_frac) {
    vm <- fit_model(formula, floor_frac)
    if (is.null(vm)) return(NULL)
    out <- tryCatch(krige(formula, train_sp, pred_sp, model = vm,
                          nmax = nmax_val, debug.level = 0),
                    error = function(e) { message("  krige error: ", conditionMessage(e)); NULL })
    if (is.null(out)) return(NULL)
    if (sum(is.finite(out$var1.pred)) == 0) return(NULL)
    list(pred = out$var1.pred, var = out$var1.var,
         model = as.character(vm$model[nrow(vm)]), nugget = vm$psill[1])
}

# Scan a candidate field (hurs units) for meaning: enough finite cells, almost
# all within the physical band, and not degenerate (non-zero spread).
scan_meaningful <- function(pred_hurs, tag) {
    n_all <- length(pred_hurs)                      # = prediction-domain cells
    fin   <- is.finite(pred_hurs)
    nfin  <- sum(fin)
    cov_frac <- if (n_all > 0) nfin / n_all else 0  # domain coverage of THIS result
    if (nfin == 0) { log_info(paste0("  scan[", tag, "]: REJECT (no finite values)")); return(FALSE) }
    vals <- pred_hurs[fin]
    tol  <- valid_tol_frac * (valid_max - valid_min)
    out_frac <- mean(vals < (valid_min - tol) | vals > (valid_max + tol))
    spread   <- sd(vals)
    ok_cov    <- cov_frac >= min_coverage           # enough of the domain filled
    ok_bounds <- out_frac <= 0.01                   # values physically sensible
    ok_spread <- is.finite(spread) && spread > 1e-6 # not degenerate
    ok <- ok_cov && ok_bounds && ok_spread
    log_info(sprintf("  scan[%s]: coverage %.2f%%, range %.2f..%.2f, out-of-band %.2f%%, sd %.3f -> %s",
                     tag, 100 * cov_frac, min(vals), max(vals), 100 * out_frac, spread,
                     ifelse(ok, "OK", "REJECT")))
    ok
}

krig_start <- Sys.time()
predictions <- NULL; se_values <- NULL; method_used <- NA_character_

# Ordered attempts. Each: label, function returning list(pred_hurs, se_hurs).
attempt_result <- function(res) {
    if (is.null(res)) return(NULL)
    list(pred = res$pred * data_sd + data_mean,
         se   = sqrt(pmax(res$var, 0)) * data_sd)
}

# 1. Full autofit (incl. Gaussian), adaptive nmax, no nugget floor.
log_info(sprintf("  Attempt 1: adaptive UK (full models, no floor, nmax=%d)", nmax_ad))
res1 <- try_krige(f_uk, nmax_ad, floor_frac = 0)
cand <- attempt_result(res1)
if (!is.null(cand) && scan_meaningful(cand$pred, "adaptive")) {
    predictions <- cand$pred; se_values <- cand$se
    method_used <- sprintf("adaptive UK full-models nmax=%d (model %s)", nmax_ad, res1$model)
}

# 2. Full model set with an enforced nugget floor; escalate the floor and step
#    nmax down. The floor guarantees a non-zero nugget so the system is solvable.
if (is.null(predictions)) {
    for (ff in c(0.02, 0.05)) {
        for (nm in nmax_steps) {
            log_info(sprintf("  Attempt: UK full-models + nugget floor %.0f%%, nmax=%d",
                             100 * ff, nm))
            res <- try_krige(f_uk, nm, floor_frac = ff)
            cand <- attempt_result(res)
            if (!is.null(cand) && scan_meaningful(cand$pred, sprintf("floor%.0f_nmax%d", 100*ff, nm))) {
                predictions <- cand$pred; se_values <- cand$se
                method_used <- sprintf("UK full-models + floor %.0f%%, nmax=%d (model %s)",
                                       100 * ff, nm, res$model)
                break
            }
        }
        if (!is.null(predictions)) break
    }
}

# 3. Covariate regression trend (lm).
if (is.null(predictions)) {
    log_info("  Attempt: covariate regression trend (lm ~ covariates)")
    lm_fit <- tryCatch(lm(f_uk, data = as.data.frame(training_df)),
                       error = function(e) { message("  lm error: ", conditionMessage(e)); NULL })
    if (!is.null(lm_fit)) {
        pr_std <- tryCatch(predict(lm_fit, newdata = as.data.frame(pred_df)),
                           error = function(e) NULL)
        if (!is.null(pr_std)) {
            cand <- list(pred = pr_std * data_sd + data_mean,
                         se   = rep(NA_real_, length(pr_std)))
            # last covariate-based resort: accept if finite (clamp handles range)
            if (sum(is.finite(cand$pred)) > 0) {
                scan_meaningful(cand$pred, "lm_trend")   # logged, not gating here
                predictions <- cand$pred; se_values <- cand$se
                method_used <- "covariate regression trend (lm)"
            }
        }
    }
}

krig_end <- Sys.time()
log_info("  Cascade done in", round(difftime(krig_end, krig_start, units = "mins"), 2), "min")

if (is.null(predictions)) fail("All covariate methods failed to produce a meaningful field.")

# Clamp accepted field to the physical range (removes small kriging overshoot).
predictions[is.finite(predictions) & predictions < valid_min] <- valid_min
predictions[is.finite(predictions) & predictions > valid_max] <- valid_max

log_info("  Method used:", method_used)
log_info("  Prediction range:",
         round(min(predictions, na.rm = TRUE), 3), "to",
         round(max(predictions, na.rm = TRUE), 3))

# --- Place predictions back on the exact target grid (GMTED_DE[[2]]) ---
pred_raster <- GMTED_DE[[2]]
values(pred_raster) <- NA_real_
se_raster <- pred_raster

cells <- cellFromXY(pred_raster, cbind(pred_df$x, pred_df$y))
pred_raster[cells] <- predictions
se_raster[cells]   <- se_values

n_finite <- sum(is.finite(values(pred_raster)))
log_info("  Valid pixels:", n_finite, "/", ncell(pred_raster))
if (n_finite == 0) fail("No finite predictions produced.")

# ----------------------------------------------------------------------------
# COVERAGE GATE (report both denominators; driver enforces on the selected one)
# ----------------------------------------------------------------------------
#   domain = cells where all covariates are finite (data is expected here)
#   grid   = every cell in the target bounding box
# coverage = finite predictions / denominator ; NaN fraction = 1 - coverage.
n_domain <- nrow(pred_df)            # expected-data cells (complete covariates)
n_grid   <- ncell(pred_raster)       # all cells in the target grid
cov_domain <- if (n_domain > 0) n_finite / n_domain else 0
cov_grid   <- if (n_grid   > 0) n_finite / n_grid   else 0
cov_domain <- min(cov_domain, 1)     # guard tiny rounding above 1
sel_cov <- if (coverage_denominator == "grid") cov_grid else cov_domain
cov_pass <- (sel_cov + 1e-9) >= min_coverage

# Machine-readable tokens for the driver (do not change these key names).
cat(sprintf("COVERAGE_DOMAIN=%.6f\n", cov_domain))
cat(sprintf("COVERAGE_GRID=%.6f\n",   cov_grid))
cat(sprintf("COVERAGE_SELECTED=%.6f (denominator=%s)\n", sel_cov, coverage_denominator))
cat(sprintf("COVERAGE_MIN=%.6f\n", min_coverage))
cat(sprintf("COVERAGE_PASS=%s\n", ifelse(cov_pass, "TRUE", "FALSE")))
log_info(sprintf("  Coverage [%s]: %.2f%% (NaN %.2f%%); threshold %.2f%% -> %s",
                 coverage_denominator, 100 * sel_cov, 100 * (1 - sel_cov),
                 100 * min_coverage, ifelse(cov_pass, "PASS", "FAIL")))
log_info(sprintf("  (domain coverage %.2f%%, grid coverage %.2f%%)",
                 100 * cov_domain, 100 * cov_grid))

crs(pred_raster) <- CRS("+proj=longlat +datum=WGS84")
crs(se_raster)   <- CRS("+proj=longlat +datum=WGS84")

# --- Write outputs (NetCDF, variable named for the target) ---
log_info("  Writing output:", basename(output_file))
names(pred_raster) <- variable_name
writeRaster(pred_raster, output_file, overwrite = TRUE,
            format = "CDF", varname = variable_name)

if (se_file != "NULL" && nchar(se_file) > 0) {
    log_info("  Writing SE:", basename(se_file))
    names(se_raster) <- paste0("se_", variable_name)
    writeRaster(se_raster, se_file, overwrite = TRUE,
                format = "CDF", varname = paste0("se_", variable_name))
}

# ============================================================================
# STEP 8: Verify outputs
# ============================================================================

cat("\nStep 8: Verifying outputs...\n")

if (file.exists(output_file)) {
    file_size <- file.info(output_file)$size / 1024^2
    cat("  Downscaled data:", basename(output_file), "(", round(file_size, 2), "MB )\n")
} else {
    fail("Output file not created: ", output_file)
}

if (se_file != "NULL" && file.exists(se_file)) {
    file_size <- file.info(se_file)$size / 1024^2
    cat("  Standard error:", basename(se_file), "(", round(file_size, 2), "MB )\n")
}

cat("\n======================================================================\n")
cat("DOWNSCALING COMPLETE\n")
cat("======================================================================\n")
cat("Variable:", variable_name, "\n")
cat("Month:", month_number, "\n")
cat("Iteration:", iteration, "\n")
cat("Covariates used:", length(covariate_list) + 1, "(DEM +", length(covariate_list), ")\n")
cat("Output files:\n")
cat("  -", output_file, "\n")
if (se_file != "NULL") {
    cat("  -", se_file, "\n")
}
cat("======================================================================\n")
