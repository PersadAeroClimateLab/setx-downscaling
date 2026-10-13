#!/usr/bin/env Rscript
#==============================================================================
# ADAPTIVE KRIGING TEMPLATE - ELEVATION FAMILY
#==============================================================================
# Single-covariate adaptive kriging for the elevation family:
#   tas, tasmin, tasmax, sfcWind   (covariate = GMTED elevation .RData)
# The covariate type is auto-detected from covariate_base_path (.RData ->
# elevation, .nc -> monthly NetCDF covariate).
#
# Two-step downscaling: 28km -> 9km -> 1km. Engine: automap::autoKrige + gstat.
#
# Reports coverage on two denominators as machine-readable tokens
# (COVERAGE_DOMAIN=, COVERAGE_GRID=, COVERAGE_PASS=) that the driver parses and
# enforces; the driver owns deletion on a gate failure. Threshold and
# denominator are config-driven ({{MIN_COVERAGE}}, {{COVERAGE_DENOMINATOR}};
# default 0.9988 / domain). The cascade keeps the best field it produces and
# reports it even below threshold.
#
# Cascading strategies (each tried until one meets the coverage gate):
#   Strategy 1: autoKrige adaptive     Strategy 2: explicit variogram models
#   Strategy 3: autoKrige boosted nmax Strategy 4: autoKrige reduced nmax + high nugget
#==============================================================================

suppressMessages({
  library(sp)
  library(raster)
  library(automap)
  library(gstat)
})

#==============================================================================
# CONFIGURATION
#==============================================================================

input_file <- "{{INPUT_FILE}}"
output_file <- "{{OUTPUT_FILE}}"
se_file <- "{{SE_FILE}}"
variable_name <- "{{VARIABLE}}"
target_resolution <- {{TARGET_RESOLUTION}}
nmax_base <- {{NMAX}}
cores <- {{CORES}}
iteration <- {{ITERATION}}
month_number <- {{MONTH_NUMBER}}
covariate_base_path <- "{{COVARIATE_BASE_PATH}}"
gmted_file <- "{{GMTED_FILE}}"
kriging_equation <- "{{KRIGING_EQUATION}}"

extent_lonmin <- {{EXTENT_LONMIN}}
extent_lonmax <- {{EXTENT_LONMAX}}
extent_latmin <- {{EXTENT_LATMIN}}
extent_latmax <- {{EXTENT_LATMAX}}

# Production coverage gate (config-driven; driver enforces on the tokens below).
min_coverage <- suppressWarnings(as.numeric("{{MIN_COVERAGE}}"))
if (is.na(min_coverage)) min_coverage <- 0.9988
coverage_denominator <- "{{COVERAGE_DENOMINATOR}}"
if (!(coverage_denominator %in% c("domain", "grid"))) coverage_denominator <- "domain"

#==============================================================================
# HELPER FUNCTIONS
#==============================================================================

get_utm_zone <- function(lon) {
  zone <- floor((lon + 180) / 6) + 1
  proj4 <- sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", zone)
  return(list(zone = zone, proj4 = proj4))
}

# Build rasters from predictions and evaluate coverage on BOTH denominators.
build_and_check <- function(predictions, variance, pred_coords, data_sd, data_mean,
                            label, n_domain_cells) {
  pred_raster <- rasterFromXYZ(data.frame(x = pred_coords$x, y = pred_coords$y, z = predictions))
  se_raster <- rasterFromXYZ(data.frame(x = pred_coords$x, y = pred_coords$y, z = sqrt(variance) * data_sd))
  crs(pred_raster) <- CRS("+proj=longlat +datum=WGS84")
  crs(se_raster) <- CRS("+proj=longlat +datum=WGS84")

  n_finite <- sum(is.finite(values(pred_raster)))
  n_grid <- ncell(pred_raster)
  cov_grid <- if (n_grid > 0) n_finite / n_grid else 0
  cov_domain <- if (n_domain_cells > 0) n_finite / n_domain_cells else 0
  cov_domain <- min(cov_domain, 1)
  sel_cov <- if (coverage_denominator == "grid") cov_grid else cov_domain
  passed <- (sel_cov + 1e-9) >= min_coverage

  cat(sprintf("  [%s] %d/%d valid (domain %.2f%%, grid %.2f%%) -> selected[%s] %.2f%% : %s\n",
              label, n_finite, n_grid, 100 * cov_domain, 100 * cov_grid,
              coverage_denominator, 100 * sel_cov, ifelse(passed, "PASSED", "below-threshold")))
  cat(sprintf("  [%s] Prediction range: %.3f to %.3f\n",
              label, min(predictions, na.rm = TRUE), max(predictions, na.rm = TRUE)))

  return(list(pred = pred_raster, se = se_raster, n_finite = n_finite, n_grid = n_grid,
              n_domain = n_domain_cells, cov_domain = cov_domain, cov_grid = cov_grid,
              sel_cov = sel_cov, passed = passed))
}

emit_coverage_tokens <- function(res) {
  cat(sprintf("COVERAGE_DOMAIN=%.6f\n", res$cov_domain))
  cat(sprintf("COVERAGE_GRID=%.6f\n", res$cov_grid))
  cat(sprintf("COVERAGE_SELECTED=%.6f (denominator=%s)\n", res$sel_cov, coverage_denominator))
  cat(sprintf("COVERAGE_MIN=%.6f\n", min_coverage))
  cat(sprintf("COVERAGE_PASS=%s\n", ifelse(isTRUE(res$passed), "TRUE", "FALSE")))
}

emit_coverage_fail_none <- function() {
  cat("COVERAGE_DOMAIN=0.000000\n")
  cat("COVERAGE_GRID=0.000000\n")
  cat(sprintf("COVERAGE_SELECTED=0.000000 (denominator=%s)\n", coverage_denominator))
  cat(sprintf("COVERAGE_MIN=%.6f\n", min_coverage))
  cat("COVERAGE_PASS=FALSE\n")
}

#==============================================================================
# MAIN PROCESSING
#==============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
if (iteration == 3) {
  cat("ADAPTIVE KRIGING (GATED, ELEVATION) - DIRECT FALLBACK (source -> 1km)\n")
} else {
  cat("ADAPTIVE KRIGING (GATED, ELEVATION) - ITERATION", iteration, "\n")
}
cat(strrep("=", 70), "\n")
cat("Variable:", variable_name, "\n")
cat("Month:", month_number, "\n")
cat("Iteration:", iteration, "\n")
cat("Target resolution:", target_resolution, "\n")
cat("nmax (base):", nmax_base, "\n")
cat(sprintf("Coverage gate: %.2f%% on the %s denominator\n", 100 * min_coverage, coverage_denominator))
cat("Input:", basename(input_file), "\n")
cat("Output:", basename(output_file), "\n")
cat("\n")

#==============================================================================
# LOAD INPUT DATA
#==============================================================================

cat("Loading input data...\n")
input_raster <- raster(input_file)
input_vals <- values(input_raster)
n_input_valid <- sum(is.finite(input_vals))
n_input_total <- length(input_vals)
cat(sprintf("Input valid cells: %d / %d (%.1f%%)\n", n_input_valid, n_input_total,
            n_input_valid / n_input_total * 100))

if (n_input_valid == 0) {
  cat("ERROR: Input raster has no valid values!\n")
  emit_coverage_fail_none()
  quit(status = 1)
}

lon_center <- (extent(input_raster)@xmin + extent(input_raster)@xmax) / 2
utm_info <- get_utm_zone(lon_center)
cat("UTM zone:", utm_info$zone, "(center lon:", round(lon_center, 2), ")\n\n")

#==============================================================================
# LOAD COVARIATE (GMTED ELEVATION, or a monthly NetCDF covariate)
#==============================================================================

is_elevation_covariate <- grepl("\\.RData$", covariate_base_path, ignore.case = TRUE)

cat("Loading GMTED:", basename(gmted_file), "\n")
gmted_list <- readRDS(gmted_file)
source_gmted <- gmted_list[[1]]
target_gmted <- gmted_list[[2]]
cat("Target GMTED:", nrow(target_gmted), "x", ncol(target_gmted), "=", ncell(target_gmted), "cells\n")

if (is_elevation_covariate) {
  cat("Covariate type: ELEVATION (GMTED)\n")
  covariate_raster <- source_gmted
} else {
  cat("Covariate type: MONTHLY NC COVARIATE\n")
  covariate_file <- sprintf(covariate_base_path, month_number)
  cat("Loading covariate:", basename(covariate_file), "\n")
  if (!file.exists(covariate_file)) {
    cat("ERROR: Covariate file not found:", covariate_file, "\n")
    emit_coverage_fail_none()
    quit(status = 1)
  }
  covariate_raster <- raster(covariate_file)
}
cat("\n")

#==============================================================================
# PREPARE TRAINING DATA
#==============================================================================

cat("Preparing training data...\n")
if (iteration == 1 || iteration == 3) {
  input_resampled <- resample(input_raster, source_gmted, method = "bilinear")
  covariate_resampled <- resample(covariate_raster, source_gmted, method = "bilinear")
} else {
  input_resampled <- input_raster
  covariate_resampled <- resample(covariate_raster, input_raster, method = "bilinear")
}

training_df <- as.data.frame(input_resampled, xy = TRUE, na.rm = TRUE)
names(training_df) <- c("x", "y", "Data")
covariate_df <- as.data.frame(covariate_resampled, xy = TRUE, na.rm = TRUE)
names(covariate_df) <- c("x", "y", "layer")
training_df <- merge(training_df, covariate_df, by = c("x", "y"))
training_df <- training_df[complete.cases(training_df), ]
cat("Training points:", nrow(training_df), "\n")

if (nrow(training_df) < 10) {
  cat("ERROR: Too few training points:", nrow(training_df), "\n")
  emit_coverage_fail_none()
  quit(status = 1)
}
cat("\n")

#==============================================================================
# PREPARE PREDICTION GRID
#==============================================================================

cat("Preparing prediction grid...\n")
prediction_df <- as.data.frame(target_gmted, xy = TRUE, na.rm = FALSE)
names(prediction_df) <- c("x", "y", "layer_template")
if (is_elevation_covariate) {
  covariate_pred <- as.data.frame(target_gmted, xy = TRUE, na.rm = FALSE)
  names(covariate_pred) <- c("x", "y", "layer")
} else {
  covariate_target <- resample(covariate_raster, target_gmted, method = "bilinear")
  covariate_pred <- as.data.frame(covariate_target, xy = TRUE, na.rm = FALSE)
  names(covariate_pred) <- c("x", "y", "layer")
}
prediction_df <- merge(prediction_df, covariate_pred, by = c("x", "y"))
prediction_df$layer_template <- NULL

# DOMAIN denominator: prediction cells with a finite covariate.
n_domain_cells <- sum(is.finite(prediction_df$layer))
n_grid_cells <- nrow(prediction_df)
cat("Prediction points:", n_grid_cells, "| domain (finite covariate):", n_domain_cells, "\n\n")

#==============================================================================
# ADAPTIVE PARAMETERS
#==============================================================================

correlation <- cor(training_df$Data, training_df$layer, use = "complete.obs")
if (abs(correlation) >= 0.50) {
  nugget_value <- 0.01; nmax_value <- nmax_base; category <- "Strong"
} else if (abs(correlation) >= 0.30) {
  nugget_value <- 0.02; nmax_value <- max(20, round(nmax_base * 0.75)); category <- "Moderate"
} else {
  nugget_value <- 0.05; nmax_value <- max(15, round(nmax_base * 0.50)); category <- "Weak"
}
cat(sprintf("Adaptive parameters: %s correlation (%.3f) -> nugget=%.3f, nmax=%d (base %d)\n",
            category, abs(correlation), nugget_value, nmax_value, nmax_base))
cat("\n")

#==============================================================================
# Z-SCORE STANDARDIZATION
#==============================================================================

data_mean <- mean(training_df$Data, na.rm = TRUE)
data_sd <- sd(training_df$Data, na.rm = TRUE)
layer_mean <- mean(training_df$layer, na.rm = TRUE)
layer_sd <- sd(training_df$layer, na.rm = TRUE)

if (data_sd < 1e-10) {
  cat("WARNING: Near-zero data SD, writing constant field\n")
  pred_raster <- target_gmted
  values(pred_raster) <- data_mean
  se_raster <- target_gmted
  values(se_raster) <- 0
  crs(pred_raster) <- CRS("+proj=longlat +datum=WGS84")
  crs(se_raster) <- CRS("+proj=longlat +datum=WGS84")
  names(pred_raster) <- variable_name
  writeRaster(pred_raster, output_file, overwrite = TRUE)
  if (se_file != "NULL" && nchar(se_file) > 0) writeRaster(se_raster, se_file, overwrite = TRUE)
  n_finite <- sum(is.finite(values(pred_raster)))
  n_grid <- ncell(pred_raster)
  const_res <- list(cov_domain = 1, cov_grid = n_finite / n_grid,
                    sel_cov = ifelse(coverage_denominator == "grid", n_finite / n_grid, 1),
                    passed = TRUE)
  emit_coverage_tokens(const_res)
  cat("SUCCESS - Constant field (zero variance)\n")
  quit(status = 0)
}

training_df$Data_std <- scale(training_df$Data)[, 1]
training_df$layer_std <- scale(training_df$layer)[, 1]
prediction_df$layer_std <- (prediction_df$layer - layer_mean) / layer_sd

#==============================================================================
# REPROJECT TO UTM
#==============================================================================

coordinates(training_df) <- ~x + y
proj4string(training_df) <- CRS("+proj=longlat +datum=WGS84")
training_utm <- spTransform(training_df, CRS(utm_info$proj4))
coordinates(prediction_df) <- ~x + y
proj4string(prediction_df) <- CRS("+proj=longlat +datum=WGS84")
prediction_utm <- spTransform(prediction_df, CRS(utm_info$proj4))
cat("Reprojection complete (UTM zone", utm_info$zone, ")\n\n")

#==============================================================================
# KRIGING - CASCADING STRATEGIES (best field kept; driver enforces the gate)
#==============================================================================

tryCatch({

  training_utm$Data <- training_utm$Data_std
  training_utm$layer <- training_utm$layer_std
  prediction_utm$layer <- prediction_utm$layer_std

  pred_coords <- as.data.frame(prediction_df)
  best_result <- NULL

  consider <- function(result) {
    if (is.null(result)) return(FALSE)
    if (is.null(best_result) || result$sel_cov > best_result$sel_cov) best_result <<- result
    isTRUE(result$passed)
  }

  # STRATEGY 1
  cat(strrep("-", 70), "\nSTRATEGY 1: autoKrige adaptive (nmax =", nmax_value, ")\n")
  tryCatch({
    kr <- autoKrige(formula = Data ~ layer, input_data = training_utm, new_data = prediction_utm,
                    nmax = nmax_value, verbose = FALSE,
                    miscFitOptions = list(merge.small.bins = TRUE, min.np.bin = 5),
                    model = c("Sph", "Exp", "Gau", "Ste"),
                    start_vals = c(psill = 0.5, range = 1, nugget = nugget_value))
    predictions <- kr$krige_output$var1.pred * data_sd + data_mean
    variance <- kr$krige_output$var1.var
    result <- build_and_check(predictions, variance, pred_coords, data_sd, data_mean, "Strategy 1", n_domain_cells)
    if (consider(result)) best_result$strategy <- "autoKrige-adaptive"
  }, error = function(e) cat(sprintf("  Strategy 1 ERROR: %s\n", e$message)))

  # STRATEGY 2
  if (is.null(best_result) || !isTRUE(best_result$passed)) {
    cat(strrep("-", 70), "\nSTRATEGY 2: explicit variogram models (nmax =", nmax_value, ")\n")
    tryCatch({
      sample_vgm <- variogram(Data ~ layer, training_utm)
      max_dist <- max(sample_vgm$dist); max_gamma <- max(sample_vgm$gamma)
      for (model_type in c("Sph", "Exp", "Gau", "Mat")) {
        done <- tryCatch({
          vgm_fit <- fit.variogram(sample_vgm,
                       vgm(psill = max_gamma * 0.9, model = model_type,
                           range = max_dist * 0.5, nugget = max_gamma * 0.1))
          kr2 <- krige(formula = Data ~ layer, locations = training_utm,
                       newdata = prediction_utm, model = vgm_fit, nmax = nmax_value)
          predictions <- kr2$var1.pred * data_sd + data_mean
          variance <- kr2$var1.var
          result <- build_and_check(predictions, variance, pred_coords, data_sd, data_mean,
                                    paste("Strategy 2", model_type), n_domain_cells)
          if (consider(result)) best_result$strategy <- paste0("explicit-", model_type)
          isTRUE(result$passed)
        }, error = function(e) { cat(sprintf("  %s ERROR: %s\n", model_type, e$message)); FALSE })
        if (isTRUE(done)) break
      }
    }, error = function(e) cat(sprintf("  Strategy 2 ERROR: %s\n", e$message)))
  }

  # STRATEGY 3
  if (is.null(best_result) || !isTRUE(best_result$passed)) {
    nmax_boosted <- min(nmax_base * 2, 100)
    cat(strrep("-", 70), "\nSTRATEGY 3: autoKrige boosted nmax (", nmax_boosted, ")\n")
    tryCatch({
      kr3 <- autoKrige(formula = Data ~ layer, input_data = training_utm, new_data = prediction_utm,
                       nmax = nmax_boosted, verbose = FALSE,
                       miscFitOptions = list(merge.small.bins = TRUE, min.np.bin = 5),
                       model = c("Sph", "Exp", "Gau", "Ste"),
                       start_vals = c(psill = 0.5, range = 1, nugget = nugget_value))
      predictions <- kr3$krige_output$var1.pred * data_sd + data_mean
      variance <- kr3$krige_output$var1.var
      result <- build_and_check(predictions, variance, pred_coords, data_sd, data_mean, "Strategy 3", n_domain_cells)
      if (consider(result)) best_result$strategy <- paste0("autoKrige-nmax", nmax_boosted)
    }, error = function(e) cat(sprintf("  Strategy 3 ERROR: %s\n", e$message)))
  }

  # STRATEGY 4
  if (is.null(best_result) || !isTRUE(best_result$passed)) {
    nmax_reduced <- 20; nugget_high <- 0.10
    cat(strrep("-", 70), "\nSTRATEGY 4: autoKrige reduced nmax (", nmax_reduced, ") + high nugget (", nugget_high, ")\n")
    tryCatch({
      kr4 <- autoKrige(formula = Data ~ layer, input_data = training_utm, new_data = prediction_utm,
                       nmax = nmax_reduced, verbose = FALSE,
                       miscFitOptions = list(merge.small.bins = TRUE, min.np.bin = 5),
                       model = c("Sph", "Exp", "Gau", "Ste"),
                       start_vals = c(psill = 0.5, range = 1, nugget = nugget_high))
      predictions <- kr4$krige_output$var1.pred * data_sd + data_mean
      variance <- kr4$krige_output$var1.var
      result <- build_and_check(predictions, variance, pred_coords, data_sd, data_mean, "Strategy 4", n_domain_cells)
      if (consider(result)) best_result$strategy <- paste0("autoKrige-nmax", nmax_reduced, "-nugget", nugget_high)
    }, error = function(e) cat(sprintf("  Strategy 4 ERROR: %s\n", e$message)))
  }

  # WRITE BEST FIELD + EMIT TOKENS
  if (is.null(best_result)) {
    cat(strrep("!", 70), "\nALL STRATEGIES FAILED to produce any field\n")
    emit_coverage_fail_none()
    quit(status = 1)
  }
  if (is.null(best_result$strategy)) best_result$strategy <- "unknown"

  names(best_result$pred) <- variable_name
  writeRaster(best_result$pred, output_file, overwrite = TRUE)
  cat("Saved prediction:", basename(output_file), "\n")
  if (se_file != "NULL" && nchar(se_file) > 0) {
    names(best_result$se) <- paste0("se_", variable_name)
    writeRaster(best_result$se, se_file, overwrite = TRUE)
    cat("Saved standard error:", basename(se_file), "\n")
  }
  if (!file.exists(output_file)) {
    cat("ERROR: Output file was not created:", output_file, "\n")
    emit_coverage_fail_none()
    quit(status = 1)
  }

  emit_coverage_tokens(best_result)
  cat(strrep("=", 70), "\n")
  cat(sprintf("Strategy used: %s\n", best_result$strategy))
  cat(sprintf("Coverage [%s]: %.2f%% (domain %.2f%%, grid %.2f%%); threshold %.2f%% -> %s\n",
              coverage_denominator, 100 * best_result$sel_cov,
              100 * best_result$cov_domain, 100 * best_result$cov_grid,
              100 * min_coverage, ifelse(isTRUE(best_result$passed), "PASS", "FAIL")))
  cat(strrep("=", 70), "\n")

}, error = function(e) {
  cat(strrep("!", 70), "\nFATAL ERROR in kriging:", e$message, "\n")
  emit_coverage_fail_none()
  quit(status = 1)
})
