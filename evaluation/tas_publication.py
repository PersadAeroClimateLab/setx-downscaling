#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
#
# Copyright (c) 2026 Dr Ifeanyichukwu Nduka
#
# This file is part of the SETx-1km-CMIP6 evaluation code release.
#   Code   : MIT License                (see LICENSE)
#   Figures: CC-BY-4.0                  (see LICENSE-CC-BY-4.0)
#   Data   : deposited separately on ESGF; see README.md for identifiers.
"""SETx-1km-CMIP6 evaluation — Daily Mean Temperature (tas) figure set.

VERSION
-------
1.0

PURPOSE
-------
Reproduces the Daily Mean Temperature (tas) evaluation figures for the companion manuscripts
describing 1 km dual-kriging (KrigR) statistical downscaling of NEX-GDDP-CMIP6
over Southeast Texas (28.50-33.75 N, 97.50-93.50 W) for the historical period
1980-2014. Downscaled products are evaluated against gridded reference data
(AORC, Daymet, ERA5) and, where the variable permits, against GHCN-D station
observations.

FIGURES PRODUCED (12)
---------------------
  Fig1_climatology_maps              Time-averaged climatology, multi-panel
  Fig2a_station_metrics              Per-station Bias / RMSE / r boxplots
  Fig3_station_validation_summary    Scatter, monthly bias, per-station bias maps
  Fig4_ensemble_diag                 Bootstrap CI, extreme indices, spread, sign agreement
  Fig5_added_value                   |Bias|/RMSE, KGE, KGE components, KGE skill score
  Fig6_taylor_seasonal               Seasonal Taylor diagrams (DJF/MAM/JJA/SON)
  4a_taylor_Annual                   Annual Taylor diagram
  FigA_kriging_se                    Kriging standard-error maps, distributions, ratio
  FigB_bias_rmse_maps                Spatial bias and RMSE maps vs AORC
  SI_S1.1_daily                      Daily Pearson-r maps
  SI_S1.2b_corr_matrix               Downscaled-only pattern-correlation matrix
  SI_S1.7_station_map                GHCN-D station network over land cover

Figures are written to OUT_DIR as paired PDF (vector) and PNG (raster) files.

TABLES PRODUCED (CSV side-cars)
-------------------------------
  Table_domain_metrics.csv           Per-model domain-averaged metrics
  Table_gridded_summary.csv          Per-model Taylor statistics by season
  Table_station_metrics.csv          Per-station metrics (station variables only)
  Fig3_scatter_regression_details.csv  Regression slope/intercept/N for Fig 3a

INPUT DATA
----------
  NEX-GDDP-CMIP6 downscaled  Zarr stores, one per CMIP6 model (1 km, dual-kriging)
  NEX-GDDP-CMIP6 native      NetCDF, extracted to the SETx domain
  ERA5 / ERA5 downscaled     NetCDF, incl. kriging standard-error fields
  AORC                       NetCDF, 1 km reference forcing
  Daymet                     NetCDF, 1 km reference (where available)
  GHCN-D                     CSV export of daily station observations
  ESA WorldCover             GeoTIFF, land-cover backdrop for the station map

  NOTE: all input paths below are SITE-SPECIFIC absolute paths on the TACC
  Lonestar6 filesystem. They MUST be edited before running elsewhere. See the
  "Paths" configuration block immediately after the imports.

REQUIREMENTS
------------
  python >= 3.10, numpy, pandas, xarray, zarr, netCDF4, scipy, matplotlib
  optional: cartopy (high-resolution coastlines), GDAL or rioxarray (WorldCover)
  See requirements.txt for pinned versions.

USAGE
-----
  python tas_publication.py [--outdir PATH]

  The script also runs unmodified inside a Jupyter kernel, in which case each
  figure is displayed inline as it is produced (see INLINE_DISPLAY_MODE).

CITATION
--------
If you use this code or the figures it produces, please cite both the software
archive and the accompanying dataset/manuscript. See CITATION.cff.
"""
import warnings,csv,struct,gc
import numpy as np,pandas as pd
import matplotlib
def _in_jupyter():
    try:
        from IPython import get_ipython
        ip=get_ipython()
        return ip is not None and "IPKernelApp" in ip.config
    except Exception:
        return False
_IN_JUPYTER=_in_jupyter()
if _IN_JUPYTER:
    try:
        get_ipython=__import__("IPython").get_ipython
        get_ipython().run_line_magic("matplotlib","inline")
    except Exception:
        try: matplotlib.use("module://matplotlib_inline.backend_inline")
        except Exception: pass
    from IPython.display import display as _ipy_display
else:
    def _ipy_display(*a,**kw): pass
import matplotlib.pyplot as plt,matplotlib.colors as mcolors
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF
from scipy import stats
from pathlib import Path
warnings.filterwarnings("ignore")
try: import xarray as xr; HAS_XR=True
except: HAS_XR=False
try: import netCDF4; HAS_NC4=True
except: HAS_NC4=False

plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","DejaVu Serif"],
    "font.size":12,"mathtext.fontset":"dejavuserif","axes.linewidth":.8,
    "axes.labelsize":13,"axes.titlesize":15,"axes.titleweight":"bold",
    "xtick.major.width":.8,"ytick.major.width":.8,"xtick.labelsize":11,
    "ytick.labelsize":11,"xtick.direction":"in","ytick.direction":"in",
    "legend.fontsize":11,"legend.framealpha":.92,"figure.dpi":1200,
    "savefig.dpi":1200,"savefig.bbox":"tight","savefig.pad_inches":.05})

FS_SUPTITLE=16   # page heading (fig.suptitle)         (was 13)
FS_TITLE=14      # figure/axes heading (ax.set_title)  (was 11)
FS_PANEL=14      # panel-letter heading ((a),(b),...)  (was 11)
FS_AXLABEL=12    # axis labels                         (was 10)
FS_TICK=10       # tick labels                         (was 8)
FS_LEGEND=11     # legend text                         (was 8)

LONGNAME="Daily Mean Temperature"
UNIT_LABEL="\u00b0C"
# Single-hue sequential colormap used for every Fig A panel, so the kriging-SE
# maps, their shared colorbar and the SE-ratio map all read as one colour
# ramp from light (low) to dark (high).
SE_CMAP="Blues"
DAYMET_ALWAYS_KEEP=False
FIG5_LEGEND_NCOL=1
# Corner used by the Fig 3 legends. Defaults suit temperature; variables whose
# data occupies the default corner override these per-variant.
FIG3A_LEGEND_LOC="upper left"
FIG3B_LEGEND_LOC="upper right"

ZARR_DIR=Path("/PATH/TO/ZARR_ROOT/tas")
ERA5_PATH=Path("/PATH/TO/EVAL_DATA/evaluation/era5/davg_t2m_era5.nc")
ERA5_DS_PATH=Path("/PATH/TO/EVAL_DATA/evaluation/downs_eval/tas_era5.nc")
ERA5_DS_SE_PATH=Path("/PATH/TO/EVAL_DATA/evaluation/downs_eval/SE_tas_era5.nc")
AORC_PATH=Path("/PATH/TO/AORC/TMP_2maboveground.nc")
DAYMET_PATH=Path("/PATH/TO/EVAL_DATA/evaluation/outputs/tas_daymet_daily_1980-2014_ncss.nc")
NEX_NATIVE_GLOB="/PATH/TO/EVAL_DATA/evaluation/outputs/extracted/tas_day_*"
# GHCN-D TAVG is derived from TMAX and TMIN station files: TAVG=(TMAX+TMIN)/2.
# Both files use the standard GHCND export schema with 'value_converted' in degC.
OBS_TMAX_CSV=Path("/PATH/TO/EVAL_DATA/Met/ghcnd_setx_output/GHCND_TMAX_1980_2020_SETx.csv")
OBS_TMIN_CSV=Path("/PATH/TO/EVAL_DATA/Met/ghcnd_setx_output/GHCND_TMIN_1980_2020_SETx.csv")
STN_META=Path("/PATH/TO/EVAL_DATA/evaluation/met/SETX_stations.csv")
GHCND_INVENTORY=Path("/PATH/TO/EVAL_DATA/evaluation/met/ghcnd-stations.txt")
BOUNDARY_SHP=Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/setx_boundary.shp")
LULC_PATH=Path("/PATH/TO/EVAL_DATA/evaluation/outputs/lulc_setx/worldcover/worldcover_setx.tif")
STATE_SHP=Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/us_states.shp")
OUT_DIR=Path("/PATH/TO/EVAL_DATA/evaluation/outputs/tas_publication"); OUT_DIR.mkdir(parents=True,exist_ok=True)

SEASONS={"DJF":{"months":[12,1,2],"label":"DJF"},"MAM":{"months":[3,4,5],"label":"MAM"},
         "JJA":{"months":[6,7,8],"label":"JJA"},"SON":{"months":[9,10,11],"label":"SON"}}
YR0,YR1=1980,2014; LAT0,LAT1,LON0,LON1=28.50,33.75,-97.50,-93.50
# Observation-data year range (set dynamically from the obs CSV in load_obs).
# Defaults to the gridded window; obs-based figures use these instead of YR0/YR1.
OBS_YR0,OBS_YR1=YR0,YR1
ERA5_C="#D55E00"; ERA5N_C="#CC79A7"; MMM_C="#0072B2"; AORC_C="#009E73"
MMMN_C="#F0E442"
DS_C="#7B1FA2"; MCMAP=plt.cm.tab20; CORR_CMAP="cividis"

# === Shapefile ===
def read_shp_polylines(p):
    p=Path(p)
    if not p.exists(): return []
    rings=[]
    with open(p,"rb") as f:
        f.read(100)
        while True:
            rh=f.read(8)
            if len(rh)<8: break
            _,cl=struct.unpack(">ii",rh); c=f.read(cl*2)
            if len(c)<4: break
            st=struct.unpack("<i",c[0:4])[0]
            if st==0: continue
            if st in (3,5,13,15) and len(c)>=44:
                np_=struct.unpack("<i",c[36:40])[0]; npt=struct.unpack("<i",c[40:44])[0]
                off=44; parts=[]
                for _ in range(np_): parts.append(struct.unpack("<i",c[off:off+4])[0]); off+=4
                pts=[]
                for _ in range(npt): x,y=struct.unpack("<dd",c[off:off+16]); pts.append((x,y)); off+=16
                parts.append(npt)
                for k in range(np_):
                    ring=np.array(pts[parts[k]:parts[k+1]])
                    if len(ring)>1: rings.append(ring)
    return rings

_COASTLINE_CACHE={"rings":None,"source":None}

_FALLBACK_COASTLINE=np.array([
    # SW corner ~Corpus Christi to Matagorda Bay
    [-97.50,28.55],[-97.45,28.58],[-97.38,28.58],[-97.28,28.56],[-97.18,28.55],
    [-97.10,28.58],[-97.05,28.65],[-97.00,28.72],[-96.93,28.78],[-96.85,28.80],
    # Matagorda Bay complex (TX) — indent inward then back out
    [-96.78,28.82],[-96.70,28.85],[-96.62,28.83],[-96.55,28.70],[-96.50,28.55],
    [-96.45,28.45],[-96.40,28.50],[-96.35,28.62],[-96.30,28.70],[-96.25,28.68],
    [-96.18,28.65],[-96.10,28.62],[-96.00,28.60],[-95.92,28.62],[-95.85,28.65],
    # Brazos/Galveston approach
    [-95.78,28.70],[-95.70,28.75],[-95.62,28.82],[-95.55,28.90],[-95.48,28.95],
    [-95.42,29.00],[-95.35,29.05],[-95.28,29.10],[-95.22,29.15],[-95.18,29.20],
    # Galveston Bay entrance (Bolivar Roads ~94.7W)
    [-95.12,29.25],[-95.08,29.32],[-95.05,29.42],[-95.02,29.55],[-95.00,29.65],
    [-94.98,29.72],[-94.95,29.78],[-94.93,29.72],[-94.92,29.55],[-94.90,29.42],
    [-94.88,29.32],[-94.85,29.30],[-94.80,29.32],[-94.78,29.38],[-94.75,29.45],
    # Trinity Bay (NE of Galveston Bay) — deeper indent inland
    [-94.78,29.55],[-94.82,29.65],[-94.78,29.75],[-94.72,29.78],[-94.68,29.70],
    [-94.65,29.55],[-94.63,29.45],[-94.62,29.38],[-94.60,29.35],[-94.55,29.38],
    [-94.50,29.40],[-94.45,29.42],[-94.40,29.43],[-94.35,29.43],[-94.30,29.45],
    # High Island / Sabine Pass approach
    [-94.25,29.48],[-94.20,29.52],[-94.15,29.55],[-94.10,29.58],[-94.05,29.62],
    [-94.00,29.65],[-93.95,29.68],[-93.90,29.70],[-93.88,29.72],[-93.85,29.72],
    # Sabine Lake entrance (~93.85W) — sharp inland indent then out
    [-93.83,29.72],[-93.82,29.80],[-93.83,29.90],[-93.84,29.95],[-93.85,29.92],
    [-93.83,29.85],[-93.80,29.78],[-93.78,29.72],[-93.75,29.70],[-93.72,29.72],
    # Calcasieu Lake (LA) entrance ~93.35W (just past SETx eastern edge at 93.5W)
    [-93.68,29.72],[-93.65,29.72],[-93.62,29.74],[-93.58,29.76],[-93.55,29.76],
    [-93.52,29.76],[-93.50,29.77],
])

def _load_coastline_rings():
    """Lazy-load coastline geometry; returns (rings, source_string).
    Priority: cartopy 10m -> cartopy 50m -> local high-res shapefiles ->
    local low-res shapefiles -> hard-coded fallback."""
    if _COASTLINE_CACHE["rings"] is not None:
        return _COASTLINE_CACHE["rings"],_COASTLINE_CACHE["source"]

    # ---------- Priority 1 & 2: cartopy Natural Earth physical coastline ----------
    # Cartopy auto-downloads NE shapefiles on first use and caches them under
    # ~/.local/share/cartopy/shapefiles/natural_earth. Both resolutions are tried
    # (10m first for journal-quality detail, 50m as fallback if 10m fetch fails).
    for res in ("10m","50m"):
        try:
            from cartopy.io import shapereader as _csr
            ne_path=_csr.natural_earth(resolution=res,category="physical",name="coastline")
            rings=read_shp_polylines(Path(ne_path))
            if rings:
                # Clip to SETx domain bounding box with a small margin to speed up plotting
                # (NE coastline is global — we'd otherwise plot the whole world's polylines).
                bbox=(LON0-0.5,LON1+0.5,LAT0-0.5,LAT1+0.5)
                clipped=[]
                for r in rings:
                    if r.size==0: continue
                    m=(r[:,0]>=bbox[0])&(r[:,0]<=bbox[1])&(r[:,1]>=bbox[2])&(r[:,1]<=bbox[3])
                    if m.any():
                        # Find contiguous spans of in-domain points to keep
                        idx=np.where(m)[0]
                        if len(idx)>=2:
                            # Split into runs where consecutive indices differ by 1
                            splits=np.where(np.diff(idx)>1)[0]+1
                            for run in np.split(idx,splits):
                                if len(run)>=2: clipped.append(r[run])
                if clipped:
                    _COASTLINE_CACHE["rings"]=clipped
                    _COASTLINE_CACHE["source"]=f"cartopy NE {res}"
                    print(f"  [coastline] loaded {len(clipped)} polylines from Cartopy "
                          f"Natural Earth {res} (high-resolution; clipped to SETx+margin)")
                    return clipped,_COASTLINE_CACHE["source"]
        except ImportError:
            if res=="10m":
                print("  [coastline] cartopy not installed; trying local shapefiles...")
            break  # no point trying 50m if cartopy itself is missing
        except Exception as e:
            print(f"  [coastline] cartopy NE {res} unavailable ({e}); trying next source")
            continue

    # ---------- Priority 3: local high-resolution shapefile candidates ----------
    hires_candidates=[
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/ne_10m_coastline.shp"),
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/ne_50m_coastline.shp"),
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/gshhs_h.shp"),
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/gshhs_f.shp"),
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/coastline.shp"),
        Path("/PATH/TO/EVAL_DATA/evaluation/shapefiles/us_coastline.shp"),
    ]
    for p in hires_candidates:
        try:
            if p.exists():
                rings=read_shp_polylines(p)
                if rings:
                    _COASTLINE_CACHE["rings"]=rings
                    _COASTLINE_CACHE["source"]=f"local high-res {p.name}"
                    print(f"  [coastline] loaded {len(rings)} polylines from local high-res shapefile {p.name}")
                    return rings,_COASTLINE_CACHE["source"]
        except Exception as e:
            print(f"  [coastline] failed to read {p.name}: {e}")

    # ---------- Priority 4: low-resolution local shapefiles (state borders / domain) ----------
    for p,label in [(STATE_SHP,"STATE_SHP"),(BOUNDARY_SHP,"BOUNDARY_SHP")]:
        try:
            if p.exists():
                rings=read_shp_polylines(p)
                if rings:
                    _COASTLINE_CACHE["rings"]=rings
                    _COASTLINE_CACHE["source"]=label
                    print(f"  [coastline] loaded {len(rings)} polylines from {label} "
                          f"({p.name}) — LOW resolution; install cartopy for higher detail")
                    return rings,label
        except Exception as e:
            print(f"  [coastline] failed to read {label}: {e}")

    # ---------- Priority 5: hard-coded fallback ----------
    _COASTLINE_CACHE["rings"]=[_FALLBACK_COASTLINE]
    _COASTLINE_CACHE["source"]="builtin"
    print(f"  [coastline] no shapefile / cartopy available; using built-in high-density "
          f"Gulf-Coast fallback ({len(_FALLBACK_COASTLINE)} waypoints)")
    return _COASTLINE_CACHE["rings"],"builtin"

def add_coastline(ax,color="black",lw=0.8,zorder=10):
    """Overlay coastline on a map axis. Robust: always renders SOMETHING.
    Uses high-resolution Natural Earth 10m data via cartopy when available.
    Returns True iff drawn."""
    rings,_=_load_coastline_rings()
    drawn=0
    for r in rings:
        try:
            ax.plot(r[:,0],r[:,1],color=color,lw=lw,zorder=zorder,clip_on=True)
            drawn+=1
        except Exception:
            pass
    return drawn>0

# === Data IO ===
def load_nc(path,vname=None):
    p=str(path)
    if HAS_XR:
        ds=xr.open_dataset(p,use_cftime=True)
        if vname is None:
            dv=[v for v in ds.data_vars if ds[v].ndim>=2 and v not in ("time_bnds",)]
            vname=dv[0] if dv else list(ds.data_vars)[0]
        data=ds[vname].values; lat=lon=time=None
        for c in ds.coords:
            cl=c.lower()
            if cl in ("lat","latitude","y"): lat=ds[c].values
            elif cl in ("lon","longitude","x"): lon=ds[c].values
            elif cl=="time":
                try:
                    tidx=ds[c].values
                    if hasattr(tidx[0],"timetuple") and not isinstance(tidx[0],(pd.Timestamp,)):
                        time=pd.DatetimeIndex([pd.Timestamp(t.year,t.month,t.day) for t in tidx])
                    else:
                        time=pd.DatetimeIndex(tidx)
                except: time=None
        ds.close()
    else:
        nc=netCDF4.Dataset(p)
        skip={"lat","latitude","lon","longitude","x","y","time","time_bnds","crs","bnds"}
        if vname is None:
            vname=[v for v in nc.variables if v.lower() not in skip and nc[v].ndim>=2][0]
        data=np.asarray(nc[vname][:]); lat=lon=time=None
        for n in nc.variables:
            nl=n.lower()
            if nl in ("lat","latitude","y"): lat=np.asarray(nc[n][:])
            elif nl in ("lon","longitude","x"): lon=np.asarray(nc[n][:])
            elif nl=="time":
                try: time=pd.DatetimeIndex(netCDF4.num2date(nc[n][:],nc[n].units))
                except: time=None
        nc.close()
    if hasattr(data,"filled"): data=data.filled(np.nan)
    return np.asarray(data,dtype=np.float64),lat,lon,time

SKIP_K_TO_C=False

def maybe_k_to_c(d,l=""):
    if d.size==0: return d
    if SKIP_K_TO_C:
        # Diagnostic so the suppression is visible in the log
        print(f"    {l} K->C check SKIPPED (SKIP_K_TO_C=True; variant not in Kelvin)")
        return d
    mn=np.nanmean(d[0] if d.ndim==3 else d)
    if mn>200: print(f"    {l} K->C (mean={mn:.1f})"); return d-273.15
    return d

ERA5_RAD_DIVIDE_BY_3600=False

def maybe_jhr_to_wm2(d,l=""):
    """Convert ERA5 radiation from J m^-2 hr^-1 to W m^-2 if flag is set.

    Auto-detects: only divides if the data magnitude is well above plausible
    W m^-2 values (mean > 5000), in case the variant flag is True but the
    input file is already in W m^-2. The autodetect keeps the function safe
    if upstream preprocessing was eventually fixed.
    """
    if not ERA5_RAD_DIVIDE_BY_3600: return d
    if d is None or d.size==0: return d
    mn=float(np.nanmean(d[0] if d.ndim==3 else d))
    if mn>5000:
        print(f"    {l} J/m^2/hr -> W/m^2 (mean was {mn:.1f}, dividing by 3600)")
        return d/3600.0
    print(f"    {l} J/m^2/hr conversion SKIPPED (mean={mn:.1f} already in W/m^2 range)")
    return d

QC_BOUNDS_MIN=-60.0
QC_BOUNDS_MAX=60.0

def physical_qc(d,label="",vmin=None,vmax=None):
    """Clip physically-implausible values in a gridded climate array.

    Values below `vmin` are clipped UP to `vmin`; values above `vmax` are
    clipped DOWN to `vmax`. NaN cells (already missing) are left untouched.
    Returns the cleaned array. The diagnostic prints the count and extreme
    values that were clipped so corrupted inputs remain visible.

    Args:
        d:     numpy array, any shape (typically 2D or 3D time x lat x lon)
        label: short string for the diagnostic message
        vmin:  lower bound; defaults to module-level QC_BOUNDS_MIN
        vmax:  upper bound; defaults to module-level QC_BOUNDS_MAX
    """
    if d is None or d.size==0: return d
    if vmin is None: vmin=QC_BOUNDS_MIN
    if vmax is None: vmax=QC_BOUNDS_MAX
    finite=np.isfinite(d)
    too_low =finite & (d<vmin)
    too_high=finite & (d>vmax)
    n_lo=int(np.count_nonzero(too_low))
    n_hi=int(np.count_nonzero(too_high))
    n_fin=int(np.count_nonzero(finite))
    if n_lo>0 or n_hi>0:
        # Build a concise one-line diagnostic with low/high breakdown.
        msgs=[]
        if n_lo>0:
            lo_min=float(np.min(d[too_low]))
            msgs.append(f"{n_lo:,} below {vmin} (min={lo_min:.3f}) -> clipped to {vmin}")
        if n_hi>0:
            hi_max=float(np.max(d[too_high]))
            msgs.append(f"{n_hi:,} above {vmax} (max={hi_max:.3f}) -> clipped to {vmax}")
        pct=100.0*(n_lo+n_hi)/max(n_fin,1)
        print(f"    {label} QC: {' | '.join(msgs)}  ({pct:.4f}% of finite cells)")
        # Work in float so we can write the clipped values back even if input
        # was integer-typed (rare for climate data, but defensive).
        if not np.issubdtype(d.dtype,np.floating):
            d=d.astype(np.float32)
        else:
            d=d.copy() if not d.flags.writeable else d
        # np.clip handles NaN correctly — NaN stays NaN, finite values are
        # bounded into [vmin, vmax].
        d=np.clip(d,vmin,vmax)
    else:
        # Quiet success line — confirms QC ran on each dataset.
        if n_fin>0:
            print(f"    {label} QC: all {n_fin:,} finite cells within [{vmin}, {vmax}]")
    return d

def ensure_lat_descending(d,lat,lon,l=""):
    if lat is not None and len(lat)>1 and lat[0]<lat[-1]:
        print(f"    {l} flip lat"); lat=lat[::-1]
        if d.ndim==3: d=d[:,::-1,:]
        elif d.ndim==2: d=d[::-1,:]
    return d,lat,lon
HAS_STATION_OBS=True

ZARR_VAR="tas"
ZARR_GLOB=f"*_hist_{ZARR_VAR}_1950-2014.zarr"

def discover_zarr_models(zarr_dir):
    zd=Path(zarr_dir); m={}
    if not zd.exists():
        print(f"  ERROR: Zarr dir not found: {zd}"); return m
    for p in sorted(zd.glob(ZARR_GLOB)):
        if not p.is_dir(): continue
        # Robustly strip the "_hist_<var>_<period>.zarr" suffix to get a clean
        # model name (e.g. "BCC-CSM2-MR"), regardless of variable/period tokens.
        name=p.name
        if name.endswith(".zarr"): name=name[:-5]
        for tok in ["_hist_","_historical_","_ssp245_","_ssp585_"]:
            if tok in name: name=name.split(tok)[0]; break
        m[name]=p
    print(f"  Found {len(m)} models (Zarr: {zd}/{ZARR_GLOB})"); return m

def load_zarr_model(zarr_path,vname=None,se_vname=None):
    """Load a CMIP6 Zarr store. vname defaults to module-level ZARR_VAR so each
    variant script reads its own variable name from the store. SE candidates
    include the variable-specific `<var>_SE` form first, then generic fallbacks."""
    if vname is None: vname=ZARR_VAR
    if not HAS_XR:
        raise RuntimeError("xarray required for zarr loading")
    try:
        ds=xr.open_zarr(str(zarr_path),consolidated=True,use_cftime=True)
    except Exception:
        ds=xr.open_zarr(str(zarr_path),consolidated=False,use_cftime=True)
    if vname not in ds.data_vars:
        dv=[v for v in ds.data_vars if ds[v].ndim>=2]
        vname=dv[0] if dv else list(ds.data_vars)[0]
    data=ds[vname].values
    se_data=None
    se_candidates=[se_vname] if se_vname else []
    se_candidates+=[f"{ZARR_VAR}_SE",f"SE_{ZARR_VAR}",f"{ZARR_VAR}_se",
                    "tas_SE","SE","se","SE_tas","tas_se","standard_error","std_error"]
    for sv in se_candidates:
        if sv and sv in ds.data_vars:
            try:
                se_data=ds[sv].values
                if hasattr(se_data,"filled"): se_data=se_data.filled(np.nan)
                se_data=np.asarray(se_data,dtype=np.float32)
                print(f"[SE:{sv}]",end=" ")
            except: pass
            break
    lat=lon=time=None
    for c in ds.coords:
        cl=c.lower()
        if cl in ("lat","latitude","y"): lat=ds[c].values
        elif cl in ("lon","longitude","x"): lon=ds[c].values
        elif cl=="time":
            try:
                tidx=ds[c].values
                if hasattr(tidx[0],"timetuple") and not isinstance(tidx[0],(pd.Timestamp,)):
                    time=pd.DatetimeIndex([pd.Timestamp(t.year,t.month,t.day) for t in tidx])
                else:
                    time=pd.DatetimeIndex(tidx)
            except: time=None
    ds.close()
    if hasattr(data,"filled"): data=data.filled(np.nan)
    return np.asarray(data,dtype=np.float32),lat,lon,time,se_data

def load_worldcover(tif_path,max_dim=2000):
    p=str(tif_path)
    try:
        from osgeo import gdal
        ds=gdal.Open(p)
        if ds is None: raise RuntimeError(f"GDAL cannot open {p}")
        gt=ds.GetGeoTransform(); nrows=ds.RasterYSize; ncols=ds.RasterXSize
        band=ds.GetRasterBand(1)
        skip=max(1, nrows//max_dim, ncols//max_dim)
        out_h=nrows//skip; out_w=ncols//skip
        data=band.ReadAsArray(buf_xsize=out_w,buf_ysize=out_h).astype(np.float32)
        lon_=np.linspace(gt[0]+gt[1]/2, gt[0]+gt[1]*(ncols-0.5), out_w)
        lat_=np.linspace(gt[3]+gt[5]/2, gt[3]+gt[5]*(nrows-0.5), out_h)
        ds=None
        print(f"    WorldCover via GDAL: {nrows}x{ncols} -> {out_h}x{out_w} (skip={skip})")
        return data,lat_,lon_
    except (ImportError,RuntimeError) as e:
        print(f"    GDAL: {e}")
    try:
        import rioxarray
        ds=xr.open_dataarray(p,engine="rasterio")
        nrows,ncols=ds.shape[-2],ds.shape[-1]
        skip=max(1,nrows//max_dim,ncols//max_dim)
        if skip>1:
            ds=ds.coarsen({ds.dims[-2]:skip,ds.dims[-1]:skip},boundary="trim").max()
        data=np.squeeze(ds.values).astype(np.float32)
        if "y" in ds.dims: lat_=ds.y.values
        elif "latitude" in ds.dims: lat_=ds.latitude.values
        else: lat_=np.arange(data.shape[0])
        if "x" in ds.dims: lon_=ds.x.values
        elif "longitude" in ds.dims: lon_=ds.longitude.values
        else: lon_=np.arange(data.shape[1])
        ds.close()
        print(f"    WorldCover via rioxarray: {data.shape} (skip={skip})")
        return data,lat_,lon_
    except ImportError:
        pass
    raise RuntimeError("No GDAL or rioxarray available to read GeoTIFF")

def crop_to_common(rnl,rnc,d):
    _,dl,dc=d.shape
    if dl==rnl and dc==rnc: return d
    r0=max((dl-rnl)//2,0); c0=max((dc-rnc)//2,0)
    out=d[:,r0:r0+rnl,c0:c0+rnc]
    if out.shape[1]<rnl or out.shape[2]<rnc:
        p=np.full((d.shape[0],rnl,rnc),np.nan); p[:,:out.shape[1],:out.shape[2]]=out; out=p
    return out

def load_obs(tmax_path,tmin_path=None,meta_path=None):
    """Load GHCN-D temperature obs and derive daily-mean TAVG = (TMAX+TMIN)/2.

    Reads TWO GHCN-D files (TMAX and TMIN), each in the standard export schema:
        station_id, station_name, latitude, longitude, elevation_m, state,
        date, value_raw, value_converted, unit, mflag, qflag, sflag
    'value_converted' holds degC. TMAX and TMIN are matched on (station, date)
    and averaged; only station-days where BOTH TMAX and TMIN exist are kept.
    Station metadata is taken from the embedded columns. Internal value field is
    'tas_C' (degC). For backward compatibility, if tmin_path is None the single
    file is treated as a ready-made daily-mean series."""
    def _find(cols,cands):
        low={c.lower():c for c in cols}
        for cand in cands:
            if cand in low: return low[cand]
        return None
    def _read_ghcnd(path,label):
        """Read one GHCN-D CSV -> DataFrame[id,date,val(degC),latitude,longitude,name].

        GHCN-D quality-control filtering applied. The CSV schema includes
        three flag columns (mflag, qflag, sflag); the QUALITY flag `qflag` is the
        one that matters for soundness: any non-blank value indicates the obs
        failed at least one of NCEI's quality checks (D=duplicate, G=gap, I=internal
        inconsistency, K=streak, M=mega-consistency, N=naught check, O=outlier,
        R=lagged-range, S=spatial-consistency, T=temporal-consistency, W=temperature
        warm-fail, X=bounds, Z=Datzilla). Such records are dropped before any
        downstream metric calculation. mflag (measurement) and sflag (source) are
        retained for now since they don't generally indicate bad data."""
        o=pd.read_csv(path); o.columns=o.columns.str.strip().str.strip('"')
        for c in o.select_dtypes("object").columns: o[c]=o[c].str.strip().str.strip('"')
        idc=_find(o.columns,["station_id","id","station","stn","ghcnd_id","staid"]) or o.columns[0]
        if idc!="id": o=o.rename(columns={idc:"id"})
        dc=_find(o.columns,["date","time","day","datetime"])
        if dc and dc!="date": o=o.rename(columns={dc:"date"})
        o["date"]=pd.to_datetime(o["date"],errors="coerce")

        # ---- QC filter on qflag ----
        # Drop any row with a non-blank quality flag. The CSV stores blanks
        # variously as empty string, NaN, "nan", or the literal characters " ";
        # any of those should be treated as "passed QC" (keep). Anything else
        # (single-letter quality-failure code) -> drop.
        qfc=_find(o.columns,["qflag","qc_flag","quality_flag","q_flag"])
        n_pre=len(o)
        if qfc is not None:
            qf=o[qfc].astype(str).str.strip()
            bad=qf.notna() & (qf!="") & (qf.str.lower()!="nan") & (qf!="0") & (qf!="0.0")
            o=o.loc[~bad].copy()
            n_drop=int(bad.sum())
            if n_drop>0:
                print(f"    {label}: QC filter dropped {n_drop:,} of {n_pre:,} records "
                      f"with non-blank qflag ({100.0*n_drop/max(n_pre,1):.2f}%)")
        # ---- end QC filter ----

        vc=_find(o.columns,["value_converted","tmax","tmin","tavg","value","value_raw"])
        if vc is None: raise KeyError(f"{label}: no value column in {list(o.columns)}")
        v=pd.to_numeric(o[vc],errors="coerce")
        o=o.assign(_v=v).dropna(subset=["date","_v"])
        # raw/tenths fields (median magnitude >60) -> scale to degC
        if vc.lower() in ("tmax","tmin","tavg","value_raw") and np.nanmedian(np.abs(o["_v"].values))>60:
            o["val"]=o["_v"]/10.0
        else: o["val"]=o["_v"]
        latc=_find(o.columns,["latitude","lat"]); lonc=_find(o.columns,["longitude","lon","long"])
        namec=_find(o.columns,["station_name","name","stnname"])
        if latc and latc!="latitude": o=o.rename(columns={latc:"latitude"})
        if lonc and lonc!="longitude": o=o.rename(columns={lonc:"longitude"})
        if namec and namec!="name": o=o.rename(columns={namec:"name"})
        if "name" not in o.columns: o["name"]=o["id"]
        keep=["id","date","val","name"]
        if "latitude" in o.columns and "longitude" in o.columns:
            o["latitude"]=pd.to_numeric(o["latitude"],errors="coerce")
            o["longitude"]=pd.to_numeric(o["longitude"],errors="coerce")
            keep+=["latitude","longitude"]
        print(f"    {label}: {len(o):,} rec (QC-passed), value column '{vc}'")
        return o[keep]

    tmx=_read_ghcnd(tmax_path,"TMAX")
    if tmin_path is not None:
        tmn=_read_ghcnd(tmin_path,"TMIN")
        # Match TMAX & TMIN on (id, date); TAVG = (TMAX+TMIN)/2 where both present
        m=tmx.merge(tmn[["id","date","val"]],on=["id","date"],how="inner",suffixes=("_max","_min"))
        m["tas_C"]=(m["val_max"]+m["val_min"])/2.0
        n_max,n_min,n_both=len(tmx),len(tmn),len(m)
        print(f"    TAVG=(TMAX+TMIN)/2: {n_both:,} station-days with both "
              f"(TMAX={n_max:,}, TMIN={n_min:,})")
        merged=m.drop(columns=["val_max","val_min"])
    else:
        merged=tmx.rename(columns={"val":"tas_C"})
    # station metadata: prefer embedded; else external meta file / inventory
    if "latitude" not in merged.columns or "longitude" not in merged.columns:
        meta=pd.read_csv(meta_path); meta.columns=meta.columns.str.strip().str.strip('"')
        for c in meta.select_dtypes("object").columns: meta[c]=meta[c].str.strip().str.strip('"')
        mid=_find(meta.columns,["station_id","id","station","stn","staid"])
        if mid and mid!="id": meta=meta.rename(columns={mid:"id"})
        mlat=_find(meta.columns,["latitude","lat"]); mlon=_find(meta.columns,["longitude","lon","long"])
        if mlat and mlat!="latitude": meta=meta.rename(columns={mlat:"latitude"})
        if mlon and mlon!="longitude": meta=meta.rename(columns={mlon:"longitude"})
        if "name" not in meta.columns:
            mname=_find(meta.columns,["station_name","name","stnname"])
            meta["name"]=meta[mname] if mname else meta.get("id","")
        meta["latitude"]=pd.to_numeric(meta["latitude"],errors="coerce")
        meta["longitude"]=pd.to_numeric(meta["longitude"],errors="coerce")
        meta=meta[["id","latitude","longitude","name"]].drop_duplicates("id").dropna(subset=["latitude","longitude"])
        missing=set(merged["id"].unique())-set(meta["id"])
        if missing and GHCND_INVENTORY.exists():
            inv=[]
            with open(GHCND_INVENTORY) as f:
                for l in f:
                    if len(l)<31: continue
                    sid=l[0:11].strip()
                    if sid in missing:
                        inv.append({"id":sid,"latitude":float(l[12:20]),"longitude":float(l[21:30]),
                                    "name":l[41:71].strip() if len(l)>41 else sid})
            if inv: meta=pd.concat([meta,pd.DataFrame(inv)],ignore_index=True).drop_duplicates("id"); print(f"    +{len(inv)} from inventory")
        merged=merged.drop(columns=[c for c in ("name",) if c in merged.columns]).merge(meta,on="id",how="left")
    merged=merged.dropna(subset=["latitude","longitude"])
    merged=merged[(merged.latitude>=LAT0)&(merged.latitude<=LAT1)&(merged.longitude>=LON0)&(merged.longitude<=LON1)&(merged.date.dt.year>=YR0)&(merged.date.dt.year<=YR1)]
    global OBS_YR0,OBS_YR1
    if len(merged):
        _yr=merged["date"].dt.year
        OBS_YR0,OBS_YR1=int(_yr.min()),int(_yr.max())
        print(f"  Obs: {len(merged):,} rec, {merged['id'].nunique()} stns, "
              f"years {OBS_YR0}-{OBS_YR1}")
        if OBS_YR0>YR0:
            print(f"  [obs-note] obs begin {OBS_YR0}; no station data for "
                  f"{YR0}-{OBS_YR0-1}. Obs-based figures now labelled {OBS_YR0}-{OBS_YR1}.")
    else: print("  Obs: 0 records after domain/year filter")
    return merged

def extract_at_stations(data_3d,lat_g,lon_g,obs_df,dates_g):
    T,nlat,nlon=data_3d.shape; la=np.asarray(lat_g); lo=np.asarray(lon_g)
    di=pd.DatetimeIndex(dates_g); parts=[]
    for sid,grp in obs_df.groupby("id"):
        ri=np.argmin(np.abs(la-grp["latitude"].iloc[0])); ci=np.argmin(np.abs(lo-grp["longitude"].iloc[0]))
        r0,r1=max(ri-1,0),min(ri+2,nlat); c0,c1=max(ci-1,0),min(ci+2,nlon)
        ts=np.nanmedian(data_3d[:,r0:r1,c0:c1],axis=(1,2))
        lookup=pd.Series(ts,index=di,dtype=np.float64); g=grp.copy()
        g["grid_val"]=lookup.reindex(g["date"].values).values; parts.append(g)
    return pd.concat(parts,ignore_index=True)

def pix_r(sim,obs,mask,mn=30):
    T=sim.shape[0]; vi=np.where(mask.ravel())[0]
    sv=sim.reshape(T,-1)[:,vi]; ov=obs.reshape(T,-1)[:,vi]
    fin=np.isfinite(sv)&np.isfinite(ov); nv=fin.sum(0)
    sf=np.where(fin,sv,0.); of_=np.where(fin,ov,0.)
    sm=np.where(nv>0,sf.sum(0)/nv,0.); om=np.where(nv>0,of_.sum(0)/nv,0.)
    sa=np.where(fin,sv-sm,0.); oa=np.where(fin,ov-om,0.)
    num=(sa*oa).sum(0); den=np.sqrt((sa**2).sum(0)*(oa**2).sum(0))
    den[den==0]=np.nan; rv=num/den; rv[nv<mn]=np.nan
    out=np.full(mask.size,np.nan); out[vi]=rv; return out.reshape(mask.shape)
def pix_r_pval(sim,obs,mask,mn=30):
    rmap=pix_r(sim,obs,mask,mn); T=sim.shape[0]; vi=np.where(mask.ravel())[0]
    sv=sim.reshape(T,-1)[:,vi]; ov=obs.reshape(T,-1)[:,vi]
    fin=np.isfinite(sv)&np.isfinite(ov); nv=fin.sum(0).astype(float)
    pout=np.full(mask.size,np.nan); rflat=rmap.ravel()[vi]
    r2=rflat**2; t_=np.abs(rflat)*np.sqrt(np.clip(nv-2,1,None))/np.sqrt(np.clip(1-r2,1e-15,None))
    p_=2*stats.t.sf(t_,np.clip(nv-2,1,None)); p_[nv<mn]=np.nan; pout[vi]=p_
    return rmap,pout.reshape(mask.shape)
def pattern_cor(x,y,w=None):
    v=np.isfinite(x)&np.isfinite(y)
    if v.sum()<10: return np.nan
    x_,y_=x[v],y[v]; w_=w[v] if w is not None else np.ones_like(x_)
    sw=w_.sum(); xw=(w_*x_).sum()/sw; yw=(w_*y_).sum()/sw
    dx=x_-xw; dy=y_-yw; num=(w_*dx*dy).sum(); den=np.sqrt((w_*dx**2).sum()*(w_*dy**2).sum())
    return num/den if den>0 else np.nan
def sp_corr(fields,mask):
    labs=list(fields); n=len(labs); C=np.ones((n,n))
    for i in range(n):
        for j in range(i+1,n): C[i,j]=C[j,i]=pattern_cor(fields[labs[i]][mask],fields[labs[j]][mask])
    return C,labs
def fld_stats(s2d,o2d,mask):
    s=s2d[mask]; o=o2d[mask]; v=np.isfinite(s)&np.isfinite(o); s,o=s[v],o[v]
    if len(s)<10: return np.nan,np.nan,np.nan
    ss,so=np.std(s,ddof=1),np.std(o,ddof=1); r=stats.pearsonr(s,o)[0]
    return(ss/so if so>0 else np.nan,r,np.sqrt(ss**2+so**2-2*ss*so*r)/so if so>0 else np.nan)
def domain_metrics(sim_3d,obs_3d,mask):
    sv=sim_3d[:,mask]; ov=obs_3d[:,mask]; fin=np.isfinite(sv)&np.isfinite(ov)
    diff=np.where(fin,sv-ov,np.nan); mbe=np.nanmean(diff)
    osum=np.nansum(np.where(fin,ov,0)); pbias=100*np.nansum(diff)/osum if osum!=0 else np.nan
    rmse=np.sqrt(np.nanmean(diff**2)); sf=sv[fin]; of_=ov[fin]
    if len(sf)<10: return dict(MBE=np.nan,PBIAS=np.nan,RMSE=np.nan,NSE=np.nan,KGE=np.nan,SD_ratio=np.nan,TX90p_bias=np.nan,TN10p_bias=np.nan)
    nse=1-np.sum((sf-of_)**2)/np.sum((of_-np.mean(of_))**2)
    N=float(len(sf)); ms=np.mean(sf); mo=np.mean(of_)
    ss=np.sqrt(np.sum((sf-ms)**2)/N); so=np.sqrt(np.sum((of_-mo)**2)/N)
    if so>0 and ss>0:
        cov=np.sum((sf-ms)*(of_-mo))/N; r=cov/(ss*so)
    else: r=np.nan
    alpha=ss/so if so>0 else np.nan
    beta=ms/mo if mo!=0 else np.nan
    kge=1-np.sqrt((r-1)**2+(alpha-1)**2+(beta-1)**2) if all(np.isfinite(x) for x in [r,alpha,beta]) else np.nan
    p90_obs=np.nanpercentile(ov,90,axis=0); p10_obs=np.nanpercentile(ov,10,axis=0)
    tx90_obs=np.nanmean(ov>p90_obs,axis=0); tx90_sim=np.nanmean(np.where(fin,sv>p90_obs,np.nan),axis=0)
    tn10_obs=np.nanmean(ov<p10_obs,axis=0); tn10_sim=np.nanmean(np.where(fin,sv<p10_obs,np.nan),axis=0)
    tx90_bias=100*(np.nanmean(tx90_sim)-np.nanmean(tx90_obs))
    tn10_bias=100*(np.nanmean(tn10_sim)-np.nanmean(tn10_obs))
    return dict(MBE=mbe,PBIAS=pbias,RMSE=rmse,NSE=nse,KGE=kge,SD_ratio=alpha,TX90p_bias=tx90_bias,TN10p_bias=tn10_bias)
def bootstrap_ci(sim_3d,obs_3d,mask,n_boot=1000):
    diff=np.where(np.isfinite(sim_3d[:,mask])&np.isfinite(obs_3d[:,mask]),sim_3d[:,mask]-obs_3d[:,mask],np.nan)
    db=np.nanmean(diff,axis=1); db=db[np.isfinite(db)]
    if len(db)<30: return np.nan,np.nan
    rng=np.random.default_rng(42); boots=np.array([np.mean(rng.choice(db,len(db))) for _ in range(n_boot)])
    return np.percentile(boots,2.5),np.percentile(boots,97.5)
def ensemble_spread(mc,obs_clim,mask):
    biases=np.stack([mc[m]-obs_clim for m in mc])
    es=np.nanstd(biases[:,mask],axis=0)
    obs_flat=obs_clim[mask]; obs_flat=obs_flat[np.isfinite(obs_flat)]
    spread_abs=np.nanmean(es)
    spread_rel=spread_abs/np.abs(np.nanmean(obs_flat)) if np.nanmean(obs_flat)!=0 else np.nan
    return spread_abs, spread_rel

# === Taylor ===
class Taylor:
    def __init__(self,fig=None,rect=111,label="AORC",sr=(0,1.65)):
        tr=PolarAxes.PolarTransform(); rl=np.array([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99,1.])
        tmax=np.pi/2; tl=np.arccos(rl)
        gh=FA.GridHelperCurveLinear(tr,extremes=(0,tmax,sr[0],sr[1]),grid_locator1=GF.FixedLocator(tl),tick_formatter1=GF.DictFormatter(dict(zip(tl,map(str,rl)))),grid_locator2=GF.MaxNLocator(6))
        if fig is None: fig=plt.figure(figsize=(5.5,5.5))
        ax=FA.FloatingSubplot(fig,rect,grid_helper=gh); fig.add_subplot(ax)
        ax.axis["top"].set_axis_direction("bottom"); ax.axis["top"].toggle(ticklabels=True,label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")
        ax.axis["top"].label.set_fontsize(12)
        ax.axis["left"].set_axis_direction("bottom"); ax.axis["left"].label.set_text("Normalized Std. Dev.")
        ax.axis["left"].label.set_fontsize(12)
        ax.axis["right"].set_axis_direction("top"); ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left"); ax.axis["bottom"].set_visible(False)
        for _axk in ("top","left","right"):
            ax.axis[_axk].major_ticklabels.set_fontsize(9)
        self.ax=ax.get_aux_axes(tr); self.fig=fig; self.pts=[]; self._ax_host=ax
        SLATE="#5C6B7A"; TAUPE="#9A8C7A"
        tt=np.linspace(0,tmax,100)
        for r_corr in (0.30,0.60,0.90):
            th=np.arccos(r_corr)
            self.ax.plot([th,th],[sr[0],sr[1]],color=SLATE,lw=0.7,ls=(0,(5,4)),
                         alpha=0.28,zorder=1)
        # --- Reference std=1 arc (heavier slate dashed, as in template) ---
        self.ax.plot(tt,np.ones_like(tt),color=SLATE,lw=1.4,ls=(0,(6,3)),alpha=0.55,zorder=2)
        # --- Centered-RMSE arcs (taupe dashed, transparent) about the reference point ---
        rs,ts_=np.meshgrid(np.linspace(sr[0],sr[1],200),np.linspace(0,tmax,200))
        crmse_field=np.sqrt(1+rs**2-2*rs*np.cos(ts_))
        ct=self.ax.contour(ts_,rs,crmse_field,levels=[0.2,0.4,0.6,0.8,1.0],
                           colors=[TAUPE],linewidths=1.0,linestyles="dashed",
                           alpha=0.55,zorder=2)
        self.ax.clabel(ct,inline=True,fontsize=9,fmt="%.1f",inline_spacing=4,colors=TAUPE)
        # Reference point marker (star) at (r=1, std=1).
        self.ax.plot([0],1.,"k*",ms=13,label=label,zorder=12,clip_on=False)
        # Axis-family colour cues: match each label to its arc family colour.
        ax.axis["top"].label.set_color(SLATE)
        ax.axis["left"].label.set_color("k")
        # Small taupe CRMSE key on an inner arc (lower-left, clear of data).
        self.ax.text(np.deg2rad(34),0.60,"Centered RMSE",ha="center",va="center",
                     fontsize=8.5,color=TAUPE,style="italic",fontweight="bold",zorder=6)
        self._SLATE=SLATE; self._TAUPE=TAUPE
    def add(self,sd,r,*a,**k):
        self.pts.append(self.ax.plot(np.arccos(np.clip(r,-1,1)),sd,*a,**k)[0])
    def legend(self,**k): self.fig.legend(self.pts,[p.get_label() for p in self.pts],**k)

INLINE_DISPLAY_MODE="figure"
_SAVED_FIGURES=[]

def _save(fig,name):
    """Write a figure to OUT_DIR as both PDF (vector) and PNG (raster), echo the
    resulting file sizes, optionally display it inline in Jupyter, then close it.

    Every figure in this script is emitted through this single function, so the
    console log doubles as a manifest of the run's outputs.
    """
    pdf_path=OUT_DIR/f"{name}.pdf"; png_path=OUT_DIR/f"{name}.png"
    for e,p in (("pdf",pdf_path),("png",png_path)): fig.savefig(p)
    try:
        pdf_sz=pdf_path.stat().st_size; png_sz=png_path.stat().st_size
        print(f"  -> {name}  [pdf {pdf_sz/1024:.0f} KB, png {png_sz/1024:.0f} KB]")
    except Exception:
        print(f"  -> {name}")
    _SAVED_FIGURES.append((name,str(png_path),str(pdf_path)))
    if _IN_JUPYTER and INLINE_DISPLAY_MODE!="none":
        try:
            if INLINE_DISPLAY_MODE=="file":
                from IPython.display import Image as _IPImage
                _ipy_display(_IPImage(filename=str(png_path),width=900,embed=True))
            else:
                _ipy_display(fig)
        except Exception:
            try: plt.show()
            except Exception: pass
    plt.close(fig)
def _nice_ticks(vmin,vmax,target=9):
    """Return ~target evenly-spaced 'nice' colorbar ticks over [vmin,vmax].
    Picks a step from a 1/2/2.5/5 x10^k ladder so labels never crowd/overlap
    regardless of the data range (fixes cramped colorbars on wide-range vars)."""
    import numpy as _np
    rng=float(vmax-vmin)
    if rng<=0: return _np.array([vmin])
    raw=rng/max(target,1); mag=10**_np.floor(_np.log10(raw))
    for m in (1,2,2.5,5,10):
        if raw<=m*mag: step=m*mag; break
    else: step=10*mag
    lo=_np.ceil(vmin/step)*step; hi=_np.floor(vmax/step)*step
    t=_np.arange(lo,hi+step*0.5,step)
    return t[(t>=vmin-1e-9)&(t<=vmax+1e-9)]
def _pl(ax,t,x=-.08,y=1.04,fs=FS_PANEL):
    ax.text(x,y,t,transform=ax.transAxes,fontsize=fs,fontweight="bold",va="top")

STANDARD_LON_TICKS=np.arange(-97.0,-93.99,1.0)
OCEAN_GREY="#D0D0D0"
def set_standard_lon_ticks(ax,labelsize=None):
    ax.set_xticks(STANDARD_LON_TICKS)
    ax.tick_params(axis="x",labelsize=labelsize if labelsize is not None else 8)
def compute_monthly_domain_mean(d3,mask,dates):
    T=d3.shape[0]; n=min(T,len(dates)); flat=d3[:n].reshape(n,-1)[:,mask.ravel()]
    daily=np.nanmean(flat,axis=1); df=pd.DataFrame({"date":dates[:n],"val":daily})
    df["ym"]=df["date"].dt.to_period("M"); gm=df.groupby("ym")["val"].mean()
    return gm.index.to_timestamp(),gm.values

# ===================== FIGURE FUNCTIONS =====================


def fig_corr_maps_combined(rmap_pairs,lat,lon,mask,period,fn):
    """Combined 1x2 Pearson r maps for two products on a
    shared color scale. rmap_pairs is a list of (source_label, rmap) tuples;
    typically 2 entries (e.g. ERA5 Downscaled and MMM-DS).

    figsize tuned so the SETx portrait-aspect maps fill the panel cells
    without large empty horizontal gaps (was 5.2x5.2 per panel which left a
    wide blank stripe between maps; now ~3.6x6.5 per panel + small wspace).
    """
    import matplotlib.gridspec as gridspec
    n=len(rmap_pairs)
    # Per-map size tuned to the SETx aspect (~4 deg lon x 5.25 deg lat).
    # Width-per-map=4.0", height=7.2", and only a thin gap between maps.
    fig=plt.figure(figsize=(4.5*n+0.6,8.4))
    gs=gridspec.GridSpec(2,n,figure=fig,height_ratios=[1.0,0.05],hspace=0.10,wspace=0.06,
                         left=0.08,right=0.97,top=0.88,bottom=0.07)
    # Shared color range across both panels for fair visual comparison
    all_r=np.concatenate([r[mask][np.isfinite(r[mask])] for _,r in rmap_pairs]) if rmap_pairs else np.array([])
    if len(all_r):
        vmin_=round(max(-1,float(np.nanmin(all_r))-.05),2)
    else: vmin_=-1
    nm=mcolors.Normalize(vmin_,1.0)
    im_last=None
    for i,(src,rmap) in enumerate(rmap_pairs):
        ax=fig.add_subplot(gs[0,i])
        ax.set_facecolor(OCEAN_GREY)
        im_last=ax.pcolormesh(lon,lat,np.where(mask,rmap,np.nan),cmap="RdYlGn",norm=nm,shading="auto",rasterized=True)
        ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1); ax.set_aspect("equal")
        ax.set_xlabel("Lon (°W)",fontsize=11); ax.tick_params(labelsize=10)
        set_standard_lon_ticks(ax)
        if i==0: ax.set_ylabel("Lat (°N)",fontsize=11)
        else: ax.set_yticklabels([])
        ax.set_title(f"({'ab'[i]}) {src}",fontsize=FS_TITLE,fontweight="bold")
    cax=fig.add_subplot(gs[1,:])
    cb=fig.colorbar(im_last,cax=cax,orientation="horizontal",extend="min")
    cb.set_label("Pearson r",fontsize=11); cb.ax.tick_params(labelsize=10)
    if "annual" in str(period).lower() or "interannual" in str(period).lower():
        _ttl=(f"Interannual Pearson Correlation Between Downscaled Climate Products "
              f"and AORC Reference Data\nfor {LONGNAME} Over the Historical Period ({YR0}-{YR1})")
    else:
        _ttl=(f"Pearson Correlation Between Downscaled Climate Products "
              f"and AORC Reference Data\nfor {LONGNAME} Over the Historical Period ({YR0}-{YR1})")
    fig.suptitle(_ttl,fontsize=FS_SUPTITLE,fontweight="bold",y=0.985)
    _save(fig,f"{fn}")


def fig_corrmat(C,labels,fn="SI_S1.2"):
    n=len(labels); fig,ax=plt.subplots(figsize=(max(7,n*.55),max(6,n*.45)),constrained_layout=True)
    vmin=max(.7,np.nanmin(C[C<1])-.02) if np.any(C[C<1]) else .7
    im=ax.imshow(C,vmin=vmin,vmax=1.,cmap=CORR_CMAP)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels,rotation=45,ha="right",fontsize=11); ax.set_yticklabels(labels,fontsize=11)
    for i in range(n):
        for j in range(n):
            v=C[i,j]; c_="white" if v<(vmin+.05) else "black"
            ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=7.5,color=c_)
    fig.colorbar(im,ax=ax,shrink=.82,label="Pattern correlation (centered)")
    ax.set_title(f"Pattern Correlation Matrix\nfor {LONGNAME} ({YR0}-{YR1})",fontsize=FS_TITLE,fontweight="bold")
    _save(fig,f"{fn}_corr_matrix")

def fig_taylor(ms,mmm_s,era5_s,season,fn):
    n=len(ms); cols=[MCMAP(i/max(n-1,1)) for i in range(n)]
    fig=plt.figure(figsize=(7,8.5)); td=Taylor(fig=fig)
    for i,(m,(ns,r)) in enumerate(ms.items()):
        if np.isfinite(ns) and np.isfinite(r):
            td.add(ns,r,marker="o",ms=7,color=cols[i],markeredgecolor="none",label=m,zorder=5)
    if mmm_s and all(np.isfinite(x) for x in mmm_s):
        td.add(mmm_s[0],mmm_s[1],marker="D",ms=10,color=MMM_C,markeredgecolor="gold",markeredgewidth=1.5,label="MMM",zorder=8)
    if era5_s and all(np.isfinite(x) for x in era5_s):
        td.add(era5_s[0],era5_s[1],marker="^",ms=10,color=ERA5_C,markeredgecolor="k",markeredgewidth=.8,label="ERA5 Downscaled",zorder=9)
    # Legend BELOW the Taylor diagram so it never overlaps data quadrant.
    # AORC reference (black star) added manually since it's a static reference.
    handles=list(td.pts)
    labels=[h.get_label() for h in handles]
    handles.append(plt.Line2D([0],[0],marker="*",color="k",ms=12,lw=0)); labels.append("AORC")
    fig.legend(handles,labels,loc="lower center",ncol=4,fontsize=11.5,
               bbox_to_anchor=(0.5,0.01),framealpha=0.5)
    fig.subplots_adjust(top=0.88,bottom=0.22)
    fig.suptitle(f"Taylor diagram ({season})\n{LONGNAME} vs AORC ({YR0}-{YR1})",fontsize=FS_SUPTITLE,fontweight="bold",y=0.97)
    _save(fig,f"{fn}_taylor_{season}")

def fig_taylor_seasonal(all_ms,all_mmm,all_era5,fn="Fig6"):
    skeys=list(SEASONS)
    fig=plt.figure(figsize=(11,12))
    fig.subplots_adjust(bottom=0.14,top=0.90,hspace=0.30)
    rects=[221,222,223,224]
    for sk in skeys:
        ms=all_ms.get(sk,{})
        skipped=[m for m,(ns,r) in ms.items() if not (np.isfinite(ns) and np.isfinite(r))]
        if skipped:
            print(f"  [Fig 6 Taylor {sk}] {len(skipped)} model(s) with NaN stats "
                  f"(skipped on diagram): {skipped}")
            for m in skipped:
                ns,r=ms[m]
                print(f"    {m}: norm_sd={ns}, r={r} -- check seasonal variance / mask coverage")
    for pi,(sk,rect) in enumerate(zip(skeys,rects)):
        ms=all_ms.get(sk,{}); mmm_s=all_mmm.get(sk,(np.nan,np.nan)); era5_s=all_era5.get(sk,(np.nan,np.nan))
        n=len(ms); cols=[MCMAP(i/max(n-1,1)) for i in range(n)]
        td=Taylor(fig=fig,rect=rect,label="AORC")
        for i,(m,(ns,r)) in enumerate(ms.items()):
            if np.isfinite(ns) and np.isfinite(r):
                td.add(ns,r,marker="o",ms=6,color=cols[i],markeredgecolor="none",label=m,zorder=5)
        if mmm_s and all(np.isfinite(x) for x in mmm_s):
            td.add(mmm_s[0],mmm_s[1],marker="D",ms=9,color=MMM_C,markeredgecolor="gold",markeredgewidth=1.2,label="MMM",zorder=8)
        if era5_s and all(np.isfinite(x) for x in era5_s):
            td.add(era5_s[0],era5_s[1],marker="^",ms=9,color=ERA5_C,markeredgecolor="k",markeredgewidth=.7,label="ERA5 Downscaled",zorder=9)
        row,col=divmod(pi,2)
        fx=0.31+col*0.5; fy=0.91-row*0.47
        fig.text(fx,fy,sk,fontsize=13,fontweight="bold",ha="center",va="bottom")
    handles=[]; labels=[]
    ms0=all_ms.get(skeys[0],{})
    for i,(m,_) in enumerate(ms0.items()):
        handles.append(plt.Line2D([0],[0],marker="o",color="w",markerfacecolor=MCMAP(i/max(len(ms0)-1,1)),ms=6)); labels.append(m)
    handles.append(plt.Line2D([0],[0],marker="D",color="w",markerfacecolor=MMM_C,ms=8,markeredgecolor="gold",markeredgewidth=1.2)); labels.append("MMM")
    handles.append(plt.Line2D([0],[0],marker="^",color="w",markerfacecolor=ERA5_C,ms=8,markeredgecolor="k",markeredgewidth=.7)); labels.append("ERA5 Downscaled")
    handles.append(plt.Line2D([0],[0],marker="*",color="k",ms=10)); labels.append("AORC")
    fig.legend(handles,labels,loc="lower center",ncol=4,fontsize=11,bbox_to_anchor=(.5,0.005))
    fig.suptitle(f"Seasonal Taylor diagrams\n{LONGNAME} vs AORC ({YR0}-{YR1})",fontsize=FS_SUPTITLE,fontweight="bold",y=0.98)
    _save(fig,f"{fn}_taylor_seasonal")




def fig_stn_metrics(sm,fn="Fig2"):
    methods=[("AORC",AORC_C),("Downscaled",DS_C),("ERA5nat",ERA5N_C),("MMM",MMM_C),("MMM_native","#888888")]
    method_labels={"AORC":"AORC","Downscaled":"ERA5 Downscaled","ERA5nat":"ERA5 Native",
                   "MMM":"MMM-Downscaled","MMM_native":"MMM Native"}
    mets=[("bias","Bias (°C)"),("rmse","RMSE (°C)"),("r","Pearson r")]
    stns=list(sm)
    fig,axes=plt.subplots(3,1,figsize=(6,10),constrained_layout=True)
    for ai,(mk,ylabel) in enumerate(mets):
        ax=axes[ai]
        data_by_method=[]; colors_plot=[]; labels_plot=[]
        for mkey,mcol in methods:
            vals=[sm[s].get(f"{mk}_{mkey}",np.nan) for s in stns]
            vals=[v for v in vals if np.isfinite(v)]
            if vals:
                data_by_method.append(vals); colors_plot.append(mcol)
                labels_plot.append(method_labels[mkey])
        if not data_by_method: continue
        bp=ax.boxplot(data_by_method,patch_artist=True,widths=.55,
                      medianprops=dict(color="k",lw=1.2),
                      flierprops=dict(marker="d",ms=4,markerfacecolor=".3",markeredgecolor=".3"))
        for patch,col in zip(bp["boxes"],colors_plot):
            patch.set_facecolor(col); patch.set_alpha(.85)
        ax.set_xticklabels(labels_plot,fontsize=11,rotation=25,ha="right")
        ax.set_ylabel(ylabel)
        if mk=="bias": ax.axhline(0,color=".4",lw=.5,ls="--")
        ax.set_title(ylabel,fontsize=FS_TITLE,fontweight="bold")
        _pl(ax,f"({'abc'[ai]})",x=-.10,y=1.05,fs=FS_PANEL)
    fig.suptitle(f"Validation for Gridded {LONGNAME}\nUsing GHCN-D Station Data ({OBS_YR0}-{OBS_YR1})",fontsize=FS_SUPTITLE,fontweight="bold",y=1.06)
    _save(fig,f"{fn}a_station_metrics")


SCATTER_LO_OVERRIDE=-30.0
SCATTER_HI_OVERRIDE=60.0
SCATTER_PAD_FRAC=0.05
SCATTER_USE_PERCENTILES=True
SCATTER_ABS_PAD=None

def _scatter_extents(obs_vals,gridded_arrays_list,pad_frac=None):
    """Return (lo, hi) suitable for a 1:1 scatter plot.

    Range source is controlled by SCATTER_USE_PERCENTILES:
      True  -> robust 1st/99th-percentile bounds (default; used for tas).
      False -> raw data min/max (used by non-tas variants).
    Padding:
      SCATTER_ABS_PAD not None -> absolute offset (raw min-pad, raw max+pad).
      else                     -> fractional pad (SCATTER_PAD_FRAC * range).
    """
    if pad_frac is None: pad_frac=SCATTER_PAD_FRAC
    pool=[np.asarray(obs_vals).ravel()]
    for g in gridded_arrays_list:
        if g is None: continue
        pool.append(np.asarray(g).ravel())
    all_v=np.concatenate(pool) if pool else np.array([])
    all_v=all_v[np.isfinite(all_v)]
    if all_v.size<10:
        lo,hi=-30.0,60.0
    elif SCATTER_ABS_PAD is not None:
        # Absolute offset around the TRUE data extrema (Task 4).
        p_lo=float(np.nanmin(all_v)); p_hi=float(np.nanmax(all_v))
        lo=p_lo-float(SCATTER_ABS_PAD); hi=p_hi+float(SCATTER_ABS_PAD)
    else:
        if SCATTER_USE_PERCENTILES:
            p_lo=float(np.percentile(all_v,1)); p_hi=float(np.percentile(all_v,99))
        else:
            p_lo=float(np.nanmin(all_v)); p_hi=float(np.nanmax(all_v))
        rng=p_hi-p_lo
        if rng<=0: rng=max(abs(p_hi),1.0)
        lo=p_lo-pad_frac*rng; hi=p_hi+pad_frac*rng
        if float(np.nanmin(all_v))>=0.0: lo=0.0
    if SCATTER_LO_OVERRIDE is not None: lo=float(SCATTER_LO_OVERRIDE)
    if SCATTER_HI_OVERRIDE is not None: hi=float(SCATTER_HI_OVERRIDE)
    return lo,hi

_LINE_STYLE_BY_LABEL={
    # 1:1 reference gets its own black dashed pattern in the scatter plots.
    "AORC":            (0,(6,2)),        # long dashes
    "ERA5 Downscaled": (0,(4,1,1,1)),    # dash-dot-dot
    "ERA5 Native":     (0,(3,2)),        # medium dashes
    "MMM-Downscaled":  (0,(2,1)),        # short dashes
    "MMM-DS":          (0,(2,1)),        # (alias)
    "MMM Downscaled":  (0,(2,1)),        # (alias with space, no hyphen)
    "MMM Native":      (0,(1,1)),        # dotted
}

def _dash_for(label):
    """Return the dash tuple for a source label; falls back to solid if unknown."""
    return _LINE_STYLE_BY_LABEL.get(label,"-")



def fig_climatology_maps(panels,lat_1km,lon_1km,mask_1km,fn="Fig1"):
    import matplotlib.gridspec as gridspec
    from scipy.ndimage import distance_transform_edt
    npan=len(panels)
    nrows_maps=int(np.ceil(npan/2))
    fig=plt.figure(figsize=(9.5,4.5*nrows_maps+0.9))
    gs=gridspec.GridSpec(nrows_maps+1,2,figure=fig,
                         height_ratios=[1]*nrows_maps+[0.07],hspace=0.30,wspace=0.22,
                         left=0.09,right=0.95,top=0.93,bottom=0.075)
    axes=[[fig.add_subplot(gs[r,c]) for c in range(2)] for r in range(nrows_maps)]
    cax=fig.add_subplot(gs[nrows_maps,:])
    # Hide any trailing unused cells (e.g. the 6th cell when only 5 panels).
    for _k in range(npan,nrows_maps*2):
        _r,_c=divmod(_k,2); axes[_r][_c].set_visible(False)

    def _fill_nan_nearest(arr2d):
        """Fill NaN pixels in a 2D array using value of nearest non-NaN pixel."""
        a=np.asarray(arr2d,dtype=np.float64)
        nan_mask=~np.isfinite(a)
        if not nan_mask.any(): return a
        if nan_mask.all(): return a   # nothing to fill from
        # distance_transform_edt returns the index of nearest valid (non-NaN) cell
        _,(yy,xx)=distance_transform_edt(nan_mask,return_distances=True,return_indices=True)
        out=a.copy(); out[nan_mask]=a[yy[nan_mask],xx[nan_mask]]
        return out

    vals_all=[]
    for p in panels:
        label,clim,plat,plon,is_coarse=p
        v=clim[mask_1km][np.isfinite(clim[mask_1km])] if not is_coarse and clim.shape==mask_1km.shape else clim[np.isfinite(clim)]
        vals_all.append(v)
    vcat=np.concatenate(vals_all) if vals_all else np.array([0])
    vmin_raw=np.nanpercentile(vcat,0.5); vmax_raw=np.nanpercentile(vcat,99.5)
    _dr=max(float(vmax_raw-vmin_raw),1e-12)
    # Want approximately 8-12 colorbar ticks; pick a step that gives that.
    # Step candidates span 6 orders of magnitude: 5e-4, 1e-3, 5e-3, 1e-2,
    # 5e-2, 1e-1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0.
    _candidates=[1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0]
    _step=_candidates[0]
    for s in _candidates:
        if _dr/s<=12:
            _step=s; break
    vmin=np.floor(vmin_raw/_step)*_step
    vmax=np.ceil(vmax_raw/_step)*_step
    # Defensive: if data is degenerate (all equal), pad to avoid Normalize error
    if vmax<=vmin: vmax=vmin+_step
    print(f"  [Fig 1] data range {vmin_raw:.4g}-{vmax_raw:.4g}, "
          f"step={_step:g}, colorbar bounds {vmin:.4g}-{vmax:.4g}")
    norm=mcolors.Normalize(vmin,vmax)
    plabs=["(a)","(b)","(c)","(d)","(e)","(f)"]
    im_last=None
    for idx,(label,clim,plat,plon,is_coarse) in enumerate(panels):
        row,col=divmod(idx,2); ax=axes[row][col]
        if is_coarse:
            from scipy.interpolate import RegularGridInterpolator as _RGI
            try:
                # Sort source axes ascending; align data accordingly
                sl=np.asarray(plat,dtype=np.float64).copy()
                dd=np.asarray(clim,dtype=np.float64).copy()
                if sl[0]>sl[-1]: sl=sl[::-1]; dd=dd[::-1,:]
                so=np.asarray(plon,dtype=np.float64).copy()
                if so[0]>so[-1]: so=so[::-1]; dd=dd[:,::-1]
                dd=_fill_nan_nearest(dd)
                interp_=_RGI((sl,so),dd,method="nearest",bounds_error=False,fill_value=np.nan)
                tgt_lat_=np.asarray(lat_1km,dtype=np.float64)
                tgt_lon_=np.asarray(lon_1km,dtype=np.float64)
                tgt_lat_asc=tgt_lat_ if tgt_lat_[0]<tgt_lat_[-1] else tgt_lat_[::-1]
                tgt_lon_asc=tgt_lon_ if tgt_lon_[0]<tgt_lon_[-1] else tgt_lon_[::-1]
                gla_,glo_=np.meshgrid(tgt_lat_asc,tgt_lon_asc,indexing="ij")
                # Clip target coordinates to source bounds for nearest-edge extrapolation
                gla_clip=np.clip(gla_,sl.min(),sl.max())
                glo_clip=np.clip(glo_,so.min(),so.max())
                disp=interp_((gla_clip,glo_clip)).astype(np.float32)
                # Flip back to original orientation
                if tgt_lat_[0]>tgt_lat_[-1]: disp=disp[::-1,:]
                if tgt_lon_[0]>tgt_lon_[-1]: disp=disp[:,::-1]
            except Exception as e:
                print(f"  [Fig 11] coarse-grid extrapolation failed for '{label}': {e}")
                disp=clim
        else:
            disp=np.where(mask_1km,clim,np.nan) if clim.shape==mask_1km.shape else clim
        ax.set_facecolor(OCEAN_GREY)
        im_last=ax.imshow(disp,extent=[LON0,LON1,LAT0,LAT1],
                          origin="lower" if lat_1km[0]<lat_1km[-1] else "upper",
                          cmap="RdYlBu_r",norm=norm,aspect="auto",interpolation="nearest")
        ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1)
        ax.set_title(f"{plabs[idx]} {label}",fontsize=FS_TITLE,fontweight="bold",pad=6)
        ax.tick_params(labelsize=10)
        set_standard_lon_ticks(ax,labelsize=9)
        _is_bottom = (idx+2>=npan)  # no panel two slots further down in same col
        if not _is_bottom: ax.set_xticklabels([])
        else: ax.set_xlabel("Lon (°W)",fontsize=11.5)
        if col>0: ax.set_yticklabels([])
        else: ax.set_ylabel("Lat (°N)",fontsize=11.5)
    # Colorbar: bounds (vmin, vmax) now equal tick endpoints, so extend='neither'
    # is safe and there are no white gaps between ticks and the color bar. Any
    # extreme outlier pixels beyond [vmin, vmax] simply saturate to the end colors.
    cb=fig.colorbar(im_last,cax=cax,orientation="horizontal",extend="neither")
    cb.set_label(f"{LONGNAME} (°C)",fontsize=11)
    cb.set_ticks(_nice_ticks(vmin,vmax)); cb.ax.tick_params(labelsize=11)
    cb.ax.set_xlim(vmin,vmax)   # belt-and-braces: also tighten the axis limits
    fig.suptitle(f"Time Averaged {LONGNAME} ({YR0}-{YR1})",fontsize=FS_SUPTITLE,fontweight="bold")
    _save(fig,f"{fn}_climatology_maps")

def fig_station_map(obs_stns,lat,lon,mask,
                    lulc_data=None,lulc_lat=None,lulc_lon=None,fn="SI_S1.7"):
    fig,ax=plt.subplots(figsize=(7,6.5),constrained_layout=True)
    used_bg=False
    if lulc_data is not None and lulc_lat is not None and lulc_lon is not None:
        wc_colors={10:"#006400",20:"#FFBB22",30:"#FFFF4C",40:"#F096FF",50:"#FA0000",
                   60:"#B4B4B4",70:"#F0F0F0",80:"#0064C8",90:"#0096A0",95:"#00CF75",100:"#FAE6A0"}
        wc_labels={10:"Tree cover",20:"Shrubland",30:"Grassland",40:"Cropland",
                   50:"Built-up",60:"Bare",70:"Snow/Ice",80:"Water",90:"Wetland",95:"Mangroves",100:"Moss/lichen"}
        try:
            lat_ok=(lulc_lat>=LAT0-.05)&(lulc_lat<=LAT1+.05)
            lon_ok=(lulc_lon>=LON0-.05)&(lulc_lon<=LON1+.05)
            if lat_ok.any() and lon_ok.any():
                li0,li1=np.where(lat_ok)[0][[0,-1]]; lj0,lj1=np.where(lon_ok)[0][[0,-1]]
                lcd=lulc_data[li0:li1+1,lj0:lj1+1]; lla=lulc_lat[li0:li1+1]; llo=lulc_lon[lj0:lj1+1]
            else: lcd,lla,llo=lulc_data,lulc_lat,lulc_lon
            max_px=2000
            if lcd.shape[0]>max_px or lcd.shape[1]>max_px:
                rstep=max(1,lcd.shape[0]//max_px); cstep=max(1,lcd.shape[1]//max_px)
                lcd=lcd[::rstep,::cstep]; lla=lla[::rstep]; llo=llo[::cstep]
            uvals=sorted(set(int(v) for v in np.unique(lcd[np.isfinite(lcd)]) if int(v) in wc_colors))
            if len(uvals)>=2:
                colors=[wc_colors[v] for v in uvals]; cmap_lc=mcolors.ListedColormap(colors)
                bounds=uvals+[uvals[-1]+1]; norm_lc=mcolors.BoundaryNorm(bounds,cmap_lc.N)
                ax.pcolormesh(llo,lla,lcd,cmap=cmap_lc,norm=norm_lc,shading="nearest",rasterized=True,zorder=1)
                import matplotlib.patches as mpatches
                patches=[mpatches.Patch(color=wc_colors[v],label=wc_labels.get(v,str(v))) for v in uvals]
                ax.legend(handles=patches,fontsize=6,loc="lower right",title="WorldCover",
                          title_fontsize=8,ncol=2,framealpha=0.5,
                          handlelength=1.0,handletextpad=0.4,columnspacing=0.8,
                          labelspacing=0.3,borderpad=0.4)
                used_bg=True
        except Exception as e: print(f"    WorldCover rendering failed: {e}")
    if not used_bg:
        ax.pcolormesh(lon,lat,np.where(mask,1.,np.nan),cmap=mcolors.ListedColormap(["#f0efe9"]),shading="auto",rasterized=True,zorder=1)
        drawn=False
        if STATE_SHP.exists():
            rings=read_shp_polylines(STATE_SHP)
            for r in rings: ax.plot(r[:,0],r[:,1],color=".35",lw=.6,zorder=2)
            if rings: drawn=True
        if not drawn and BOUNDARY_SHP.exists():
            rings=read_shp_polylines(BOUNDARY_SHP)
            for r in rings: ax.plot(r[:,0],r[:,1],color=".35",lw=.6,zorder=2)
            if rings: drawn=True
        if not drawn: ax.plot([LON0,LON1,LON1,LON0,LON0],[LAT0,LAT0,LAT1,LAT1,LAT0],"k-",lw=.5,zorder=2)
    ax.scatter(obs_stns["longitude"],obs_stns["latitude"],s=25,marker="o",
               c="k",edgecolors="none",label=f"GHCN-D stations ({len(obs_stns)})",zorder=10)
    ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1); ax.set_aspect("equal")
    ax.set_xlabel("Longitude (°W)"); ax.set_ylabel("Latitude (°N)")
    if not used_bg: ax.legend(fontsize=11,loc="lower right")
    ax.set_title(f"Positions of the GHCN-D stations for {LONGNAME}",fontsize=FS_TITLE,fontweight="bold")
    _save(fig,f"{fn}_station_map")

# ===================== NEW FIGURES =====================

def fig_bias_rmse_maps(era5ds,mmm,aorc,lat,mask,fn="FigB",obs_df=None,sm=None):
    """Fig B — 2x2 spatial bias/RMSE maps with per-row shared horizontal colorbars.

    Layout:
        (a) ERA5 Downscaled Bias        (b) MMM-Downscaled Bias    <- shared "Bias" colorbar below
        (c) ERA5 Downscaled RMSE        (d) MMM-Downscaled RMSE    <- shared "RMSE" colorbar below

    When obs_df + sm are provided, bias panels (a,b) also overlay per-station
    bias markers using the same diverging colour scale as the underlying map.
    """
    import matplotlib.gridspec as gridspec
    sources=[("ERA5 Downscaled",era5ds),("MMM-Downscaled",mmm)]
    bias_maps=[]; rmse_maps=[]
    for _,src in sources:
        bias_maps.append(np.nanmean(src-aorc,0))
        rmse_maps.append(np.sqrt(np.nanmean((src-aorc)**2,0)))
    bvals=np.concatenate([b[mask][np.isfinite(b[mask])] for b in bias_maps])
    bmax=max(abs(np.nanpercentile(bvals,1)),abs(np.nanpercentile(bvals,99)))
    rvals=np.concatenate([r[mask][np.isfinite(r[mask])] for r in rmse_maps])
    rmax=np.nanpercentile(rvals,99)

    stn_lons=stn_lats=None; b_era5ds_stn=b_mmm_stn=None
    if obs_df is not None and sm is not None:
        stn_coords=(obs_df[["id","latitude","longitude"]].dropna()
                    .drop_duplicates("id").set_index("id"))
        sids=[s for s in sm if s in stn_coords.index]
        if sids:
            stn_lons=np.array([float(stn_coords.loc[s,"longitude"]) for s in sids])
            stn_lats=np.array([float(stn_coords.loc[s,"latitude"]) for s in sids])
            b_era5ds_stn=np.array([sm[s].get("bias_Downscaled",np.nan) for s in sids])
            b_mmm_stn   =np.array([sm[s].get("bias_MMM",np.nan) for s in sids])
            all_stn=np.concatenate([b_era5ds_stn,b_mmm_stn])
            all_stn=all_stn[np.isfinite(all_stn)]
            if all_stn.size:
                bmax=max(bmax,float(np.nanpercentile(np.abs(all_stn),95)))

    fig=plt.figure(figsize=(11.5,11.2))
    gs=gridspec.GridSpec(2,3,figure=fig,
                         width_ratios=[1.0,1.0,0.05],
                         height_ratios=[1.0,1.0],
                         hspace=0.28,wspace=0.10,
                         left=0.07,right=0.92,top=0.88,bottom=0.06)

    stn_marker_kw=dict(cmap="coolwarm",vmin=-bmax,vmax=bmax,
                       s=22,edgecolor="k",lw=0.4,zorder=5)

    im_b_last=None; bias_axes=[]
    for col,(label,_) in enumerate(sources):
        ax=fig.add_subplot(gs[0,col]); bias_axes.append(ax)
        ax.set_facecolor(OCEAN_GREY)
        im_b_last=ax.imshow(np.where(mask,bias_maps[col],np.nan),extent=[LON0,LON1,LAT0,LAT1],
                            origin="lower" if lat[0]<lat[-1] else "upper",
                            cmap="coolwarm",vmin=-bmax,vmax=bmax,aspect="auto",interpolation="nearest")
        if stn_lons is not None:
            b_stn=b_era5ds_stn if col==0 else b_mmm_stn
            ax.scatter(stn_lons,stn_lats,c=b_stn,**stn_marker_kw)
        add_coastline(ax,color="black",lw=0.5,zorder=10)
        letter="ab"[col]
        ax.set_title(f"({letter}) {label}",fontsize=FS_TITLE,fontweight="bold")
        ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1); ax.set_aspect("equal"); ax.tick_params(labelsize=9)
        set_standard_lon_ticks(ax,labelsize=8)
        ax.set_xlabel("Lon (\u00b0W)",fontsize=10)
        if col==0: ax.set_ylabel("Lat (\u00b0N)",fontsize=10)
        else: ax.set_yticklabels([])

    # Bias colorbar — VERTICAL, right of the top row.
    cax_b=fig.add_subplot(gs[0,2])
    cb_b=fig.colorbar(im_b_last,cax=cax_b,orientation="vertical",extend="both")
    cb_b.set_label("Bias (\u00b0C)",fontsize=11); cb_b.ax.tick_params(labelsize=9)

    im_r_last=None; rmse_axes=[]
    for col,(label,_) in enumerate(sources):
        ax=fig.add_subplot(gs[1,col]); rmse_axes.append(ax)
        ax.set_facecolor(OCEAN_GREY)
        im_r_last=ax.imshow(np.where(mask,rmse_maps[col],np.nan),extent=[LON0,LON1,LAT0,LAT1],
                            origin="lower" if lat[0]<lat[-1] else "upper",
                            cmap="viridis",vmin=0,vmax=rmax,aspect="auto",interpolation="nearest")
        add_coastline(ax,color="white",lw=0.5,zorder=10)
        letter="cd"[col]
        ax.set_title(f"({letter}) {label}",fontsize=FS_TITLE,fontweight="bold")
        ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1); ax.set_aspect("equal"); ax.tick_params(labelsize=9)
        set_standard_lon_ticks(ax,labelsize=8)
        ax.set_xlabel("Lon (\u00b0W)",fontsize=10)
        if col==0: ax.set_ylabel("Lat (\u00b0N)",fontsize=10)
        else: ax.set_yticklabels([])

    # RMSE colorbar — VERTICAL, right of the bottom row. NOTE: no stray "RMSE"
    # text label above the bar anymore; the row super-title carries that word.
    cax_r=fig.add_subplot(gs[1,2])
    cb_r=fig.colorbar(im_r_last,cax=cax_r,orientation="vertical",extend="max")
    cb_r.set_label("RMSE (\u00b0C)",fontsize=11); cb_r.ax.tick_params(labelsize=9)

    fig.canvas.draw()
    def _row_title(axes_row,text):
        tops=[a.get_position().y1 for a in axes_row]
        xs  =[a.get_position().x0 for a in axes_row]+[a.get_position().x1 for a in axes_row]
        y=max(tops)+0.035
        xmid=(min(xs)+max(xs))/2
        fig.text(xmid,min(y,0.99),text,fontsize=FS_SUPTITLE+1,fontweight="bold",ha="center",va="bottom")
    _row_title(bias_axes,"Bias")
    _row_title(rmse_axes,"RMSE")

    fig.suptitle(f"Mean Bias and RMSE Between Downscaled Products and AORC "
                 f"for {LONGNAME} ({YR0}-{YR1})",
                 fontsize=FS_SUPTITLE,fontweight="bold",y=0.995)
    _save(fig,f"{fn}_bias_rmse_maps")


def _monthly_mean_se(d3,mask,dates,se_map=None):
    T=d3.shape[0]; n=min(T,len(dates)); flat=d3[:n].reshape(n,-1)[:,mask.ravel()]
    daily_mean=np.nanmean(flat,axis=1)
    if se_map is not None and se_map.shape==mask.shape:
        se_val=np.nanmean(se_map[mask]); daily_se=np.full(n,se_val,dtype=np.float32)
    else:
        daily_se=np.nanstd(flat,axis=1)/np.sqrt(np.sum(np.isfinite(flat),axis=1).clip(1))
    df=pd.DataFrame({"date":dates[:n],"val":daily_mean,"se":daily_se})
    df["ym"]=df["date"].dt.to_period("M")
    gm=df.groupby("ym").agg(val=("val","mean"),se=("se","mean"))
    return gm.index.to_timestamp(),gm["val"].values,gm["se"].values







def fig_se_spatial_hist(se_panels,lat,lon,mask,fn="FigA"):
    """Fig A — Kriging standard-error 2x2 layout per tas_results_v3.docx S5.2.4:
       (a) ERA5 Downscaled time-mean SE map
       (b) MMM-Downscaled time-mean SE map
       (c) Shared histogram of grid-cell SE for both products
       (d) Ratio map MMM-SE / ERA5 Downscaled-SE.
       panels (a) and (b) now SHARE a single horizontal
       colorbar below the top row (both use the same SE_CMAP + norm_se
       colormap, so a shared bar is honest and doubles as more real estate).
       Ratio panel (d) keeps its own colorbar because it uses a different
       cmap and TwoSlopeNorm centred on 1.0. All colorbars enlarged.
    """
    import matplotlib.gridspec as gridspec
    fig=plt.figure(figsize=(12,10.5))
    # Layout: 3 rows (maps, shared bar, ratio+hist), 2 cols. Height ratios
    # give the maps and bottom panels equal room and the shared colorbar
    # ~10% (thick enough to read after enlargement).
    gs=gridspec.GridSpec(3,2,figure=fig,
                         height_ratios=[1.0,0.10,1.0],
                         hspace=0.35,wspace=0.15,
                         left=0.07,right=0.97,top=0.92,bottom=0.06)
    se_by={lab:(smap,scol) for (lab,smap,scol) in se_panels}
    era5ds_label=next((k for k in se_by if "ERA5" in k and "Down" in k),None)
    mmm_label   =next((k for k in se_by if "MMM" in k),None)
    se_era5ds=se_by.get(era5ds_label,(None,ERA5_C))[0]
    se_mmm   =se_by.get(mmm_label,(None,MMM_C))[0]
    def _disp(arr):
        if arr is None: return None
        return np.where(mask,arr,np.nan) if arr.shape==mask.shape else arr
    se_e_d=_disp(se_era5ds); se_m_d=_disp(se_mmm)
    vals=[]
    for v in (se_e_d,se_m_d):
        if v is not None: vals.append(v[np.isfinite(v)])
    cat=np.concatenate(vals) if vals else np.array([0.0,0.5])
    def _nice_ceil(x):
        if not np.isfinite(x) or x<=0: return 1.0
        k=np.floor(np.log10(x)); base=x/(10**k)
        nice=1.0 if base<=1 else (2.0 if base<=2 else (5.0 if base<=5 else 10.0))
        return nice*(10**k)
    _p99=float(np.nanpercentile(cat,99)) if cat.size else 0.5
    vmax_se=_nice_ceil(_p99) if _p99>0 else 0.5
    norm_se=mcolors.Normalize(0.0,vmax_se)

    # (a) ERA5 Downscaled SE map — no per-panel colorbar (shared bar below).
    ax_a=fig.add_subplot(gs[0,0])
    ax_a.set_facecolor(OCEAN_GREY)
    im_ref=None
    if se_e_d is not None:
        im_ref=ax_a.pcolormesh(lon,lat,se_e_d,cmap=SE_CMAP,norm=norm_se,shading="auto",rasterized=True)
    ax_a.set_xlim(LON0,LON1); ax_a.set_ylim(LAT0,LAT1); ax_a.set_aspect("equal"); ax_a.tick_params(labelsize=11)
    set_standard_lon_ticks(ax_a)
    ax_a.set_xlabel("Lon (\u00b0W)",fontsize=10); ax_a.set_ylabel("Lat (\u00b0N)",fontsize=10)
    ax_a.set_title("(a) ERA5 Downscaled",fontsize=FS_TITLE,fontweight="bold")

    # (b) MMM-DS SE map — no per-panel colorbar (shared bar below).
    ax_b=fig.add_subplot(gs[0,1])
    ax_b.set_facecolor(OCEAN_GREY)
    if se_m_d is not None:
        im_b=ax_b.pcolormesh(lon,lat,se_m_d,cmap=SE_CMAP,norm=norm_se,shading="auto",rasterized=True)
        if im_ref is None: im_ref=im_b
    ax_b.set_xlim(LON0,LON1); ax_b.set_ylim(LAT0,LAT1); ax_b.set_aspect("equal"); ax_b.tick_params(labelsize=11)
    set_standard_lon_ticks(ax_b)
    ax_b.set_xlabel("Lon (\u00b0W)",fontsize=10)
    ax_b.set_yticklabels([])
    ax_b.set_title("(b) MMM-Downscaled",fontsize=FS_TITLE,fontweight="bold")

    # Shared horizontal colorbar for (a) + (b).
    if im_ref is not None:
        cax_ab=fig.add_subplot(gs[1,:])
        cb_ab=fig.colorbar(im_ref,cax=cax_ab,orientation="horizontal",extend="max")
        cb_ab.set_label(f"SE ({UNIT_LABEL})",fontsize=12); cb_ab.ax.tick_params(labelsize=11)

    def _fmt_mean(x):
        ax_=abs(x)
        if ax_==0 or not np.isfinite(x): return f"{x:.2f}"
        if ax_>=1: return f"{x:.2f}"
        if ax_>=0.01: return f"{x:.3f}"
        if ax_>=1e-4: return f"{x:.5f}"
        return f"{x:.2e}"

    # (c) Shared histogram
    ax_c=fig.add_subplot(gs[2,0])
    if se_e_d is not None:
        v_e=se_e_d[np.isfinite(se_e_d)]
        ax_c.hist(v_e,bins=60,color=ERA5_C,alpha=0.65,edgecolor="w",lw=0.3,
                  label=f"ERA5 Downscaled (mean={_fmt_mean(np.nanmean(v_e))} {UNIT_LABEL})")
    if se_m_d is not None:
        v_m=se_m_d[np.isfinite(se_m_d)]
        ax_c.hist(v_m,bins=60,color=MMM_C,alpha=0.65,edgecolor="w",lw=0.3,
                  label=f"MMM-Downscaled (mean={_fmt_mean(np.nanmean(v_m))} {UNIT_LABEL})")
    ax_c.set_xlabel(f"SE ({UNIT_LABEL})"); ax_c.set_ylabel("Grid count")
    ax_c.legend(fontsize=10,loc="upper right",framealpha=0.5)
    ax_c.set_title("(c) Grid-cell SE distributions",fontsize=FS_TITLE,fontweight="bold")
    ax_c.grid(True,ls=":",lw=0.4,alpha=0.5)

    # (d) Ratio map MMM-SE / ERA5 Downscaled-SE — enlarged colorbar.
    ax_d=fig.add_subplot(gs[2,1])
    ax_d.set_facecolor(OCEAN_GREY)
    if se_e_d is not None and se_m_d is not None:
        with np.errstate(divide="ignore",invalid="ignore"):
            # Guard the denominator against near-zero SE. An absolute floor
            # (1e-9) is meaningless across variables whose SE magnitudes span
            # orders of magnitude: for hurs the ERA5 SE is O(1e-4), so cells
            # with a merely small denominator survived and produced ratios up
            # to ~1e7, which dominated the percentiles and saturated the map.
            # Threshold relative to the field's own median instead.
            _den_med=float(np.nanmedian(se_e_d[np.isfinite(se_e_d)])) if np.any(np.isfinite(se_e_d)) else 0.0
            _den_min=max(_den_med*0.05,1e-12)
            ratio=np.where(se_e_d>_den_min,se_m_d/se_e_d,np.nan)
        _rv=ratio[np.isfinite(ratio)]
        if _rv.size:
            # Sequential single-hue scale spanning the robust 2nd-98th
            # percentile range of the ratio actually present in the field.
            _vmin=float(np.nanpercentile(_rv,2)); _vmax=float(np.nanpercentile(_rv,98))
            if not np.isfinite(_vmin) or not np.isfinite(_vmax) or _vmax<=_vmin:
                _vmin,_vmax=0.0,2.0
            if _vmax-_vmin < 1e-6: _vmax=_vmin+1e-6
        else:
            _vmin,_vmax=0.0,2.0
        norm_r=mcolors.Normalize(vmin=_vmin,vmax=_vmax)
        print(f"  [FigA(d)] SE-ratio colour range {_vmin:.3g}-{_vmax:.3g} "
              f"(denominator floor {_den_min:.3g}; {_rv.size:,} valid cells)")
        im_d=ax_d.pcolormesh(lon,lat,ratio,cmap=SE_CMAP,norm=norm_r,shading="auto",rasterized=True)
        cb_d=fig.colorbar(im_d,ax=ax_d,orientation="horizontal",shrink=0.95,pad=0.18,
                          label="MMM-Downscaled SE / ERA5 Downscaled SE",extend="both")
        cb_d.ax.tick_params(labelsize=10); cb_d.set_label("MMM-Downscaled SE / ERA5 Downscaled SE",fontsize=11)
        mean_ratio=float(np.nanmean(ratio[mask])) if ratio.shape==mask.shape else float(np.nanmean(ratio))
        print(f"  [FigA(d)] SE ratio domain mean = {mean_ratio:.3f}")
    ax_d.set_xlim(LON0,LON1); ax_d.set_ylim(LAT0,LAT1); ax_d.set_aspect("equal")
    ax_d.tick_params(labelsize=10)
    set_standard_lon_ticks(ax_d,labelsize=10)
    ax_d.set_xlabel("Lon (\u00b0W)",fontsize=11); ax_d.set_ylabel("Lat (\u00b0N)",fontsize=11)
    ax_d.set_title("(d) SE ratio Between Downscaled MMM and ERA5",
                   fontsize=FS_TITLE,fontweight="bold")

    fig.suptitle(f"Kriging Standard Error Results for {LONGNAME} ({YR0}-{YR1})",
                 fontsize=FS_SUPTITLE,fontweight="bold",y=0.985)
    _save(fig,f"{fn}_kriging_se")

def fig_combined_stn_panel(ov,av,dv,e5v,monthly_mse,sm,obs_df,
                           mmm_v=None,mmmn_v=None,fn="Fig3"):
    """Station-validation summary (1-row layout, three panel cells).
    (a) Daily scatter: per-point alpha scatter for AORC/ERA5 Downscaled/ERA5 Native (and
        MMM-Downscaled if supplied), 1:1 line, per-product OLS regression,
        annotated r/slope/intercept/N.
    (b) Monthly mean bias time series with GHCN-D as baseline, +-0.5 deg C
        practical-equivalence band (grey).
    (c) Per-station mean bias spatial maps, 1x4 layout for AORC, ERA5 Downscaled,
        ERA5 Native, and MMM-Downscaled, shared diverging RdBu_r scale.
    """
    import matplotlib.gridspec as gridspec
    fig=plt.figure(figsize=(16,11))
    gs=gridspec.GridSpec(2,1,figure=fig,height_ratios=[1.0,0.95],hspace=0.42,
                         left=0.06,right=0.97,top=0.90,bottom=0.055)
    gs_top=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[0],wspace=0.22,width_ratios=[1.0,1.4])

    # ----- (a) Daily scatter -----
    ax_scat=fig.add_subplot(gs_top[0,0])
    sources=[(av,"AORC",AORC_C),(dv,"ERA5 Downscaled",ERA5_C)]
    if e5v is not None: sources.append((e5v,"ERA5 Native",ERA5N_C))
    if mmm_v is not None: sources.append((mmm_v,"MMM-Downscaled",MMM_C))
    if mmmn_v is not None: sources.append((mmmn_v,"MMM Native",MMMN_C))
    lo,hi=_scatter_extents(ov,[s[0] for s in sources])
    scat_rows=[("Source","r","slope","intercept","N")]
    for si,(gv,lab,col) in enumerate(sources):
        v=np.isfinite(ov)&np.isfinite(gv); o,g=ov[v],gv[v]
        if len(o)<10: continue
        ax_scat.scatter(o,g,s=2,c=col,alpha=0.08,lw=0,rasterized=True,zorder=2+si)
        sl,ic,r,_,_=stats.linregress(o,g)
        ax_scat.plot([lo,hi],[sl*lo+ic,sl*hi+ic],color=col,lw=2.0,ls=_dash_for(lab),
                     alpha=0.9,label=f"{lab}: r={r:.2f}",zorder=10+si)
        scat_rows.append((lab,f"{r:.4f}",f"{sl:.4f}",f"{ic:+.4f}",f"{len(o)}"))
    ax_scat.plot([lo,hi],[lo,hi],color="k",lw=.8,ls=(0,(5,3)),alpha=.55,label="1:1",zorder=1)
    ax_scat.set_xlim(lo,hi); ax_scat.set_ylim(lo,hi); ax_scat.set_aspect("equal")
    ax_scat.set_xlabel("Station Data (GHCN-D) (\u00b0C)"); ax_scat.set_ylabel("Gridded Datasets (\u00b0C)")
    ax_scat.legend(fontsize=11,loc=FIG3A_LEGEND_LOC,framealpha=0.55)
    # Panel letter is carried by the title itself (loc="left"), which keeps
    # "(a)" locked to the start of its label at any figure size instead of
    # floating on a separate transform.
    # The variable name is already carried by the figure suptitle, so the panel
    # titles stay short. Long names (e.g. "Daily Average Relative Humidity")
    # previously pushed the (a) title into the (b) title.
    ax_scat.set_title("(a) Gridded Data vs GHCN-D Stations",
                      fontsize=FS_TITLE-1,fontweight="bold",loc="left",pad=10)
    # Write the regression-details CSV (slope, intercept, N) next to the figure.
    try:
        with open(OUT_DIR/f"{fn}_scatter_regression_details.csv","w",newline="") as _fdet:
            _w=csv.writer(_fdet)
            for row in scat_rows: _w.writerow(row)
    except Exception as _e:
        print(f"  [Fig3] could not write regression-details CSV: {_e}")

    # ----- (b) Monthly mean bias with +-0.5 deg C equivalence band -----
    ax_bias=fig.add_subplot(gs_top[0,1])
    if "GHCN-D obs" in monthly_mse:
        od,ov_,_=monthly_mse["GHCN-D obs"]
        obs_s=pd.Series(np.asarray(ov_,dtype=float),
                        index=pd.to_datetime(od).to_period("M").to_timestamp())
        obs_s=obs_s[~obs_s.index.duplicated(keep="first")].sort_index()
        style={"AORC":            (AORC_C, _LINE_STYLE_BY_LABEL["AORC"]),
               "ERA5 downscaled": (ERA5_C, _LINE_STYLE_BY_LABEL["ERA5 Downscaled"]),
               "ERA5 Native":     (ERA5N_C,_LINE_STYLE_BY_LABEL["ERA5 Native"]),
               "MMM-Downscaled":  (MMM_C,  _LINE_STYLE_BY_LABEL["MMM-Downscaled"]),
               "MMM Native":      (MMMN_C, _LINE_STYLE_BY_LABEL["MMM Native"])}
        for label,(dates,vals,se) in monthly_mse.items():
            if label=="GHCN-D obs": continue
            col,ls=style.get(label,(".5",(0,(4,2))))
            prod_s=pd.Series(np.asarray(vals,dtype=float),
                             index=pd.to_datetime(dates).to_period("M").to_timestamp())
            prod_s=prod_s[~prod_s.index.duplicated(keep="first")].sort_index()
            bias=prod_s-obs_s.reindex(prod_s.index)
            ax_bias.plot(prod_s.index,bias.values,color=col,ls=ls,lw=1.3,
                         alpha=0.85,label=label)
        ax_bias.axhline(0,color=".4",lw=.5,ls="--")
    ax_bias.set_xlim(pd.Timestamp(f"{OBS_YR0}-01-01"),pd.Timestamp(f"{OBS_YR1}-12-31"))
    ax_bias.set_ylabel("Bias (\u00b0C)"); ax_bias.set_xlabel("Date")
    ax_bias.legend(fontsize=10,ncol=2,framealpha=0.3,loc=FIG3B_LEGEND_LOC)
    ax_bias.set_title("(b) Domain-Averaged Monthly Bias vs GHCN-D Stations",
                      fontsize=FS_TITLE-1,fontweight="bold",loc="left",pad=10)

    gs_c=gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=gs[1],
                                          height_ratios=[1,0.08],hspace=0.42,wspace=0.10)
    stn_coords=(obs_df[["id","latitude","longitude"]].dropna()
                .drop_duplicates("id").set_index("id"))
    sids_c=[s for s in sm if s in stn_coords.index]
    lons_c=np.array([float(stn_coords.loc[s,"longitude"]) for s in sids_c])
    lats_c=np.array([float(stn_coords.loc[s,"latitude"]) for s in sids_c])
    b_aorc_c=np.array([sm[s].get("bias_AORC",np.nan) for s in sids_c])
    b_ds_c  =np.array([sm[s].get("bias_Downscaled",np.nan) for s in sids_c])
    b_e5_c  =np.array([sm[s].get("bias_ERA5nat",np.nan) for s in sids_c])
    b_mmm_c =np.array([sm[s].get("bias_MMM",np.nan) for s in sids_c])
    all_b_c=np.concatenate([b_aorc_c,b_ds_c,b_e5_c,b_mmm_c])
    vmax_c=max(0.5,float(np.nanpercentile(np.abs(all_b_c[np.isfinite(all_b_c)]),95)))
    # Cap at +-2 deg C per the manuscript spec ("shared diverging colour scale, -2 to +2 deg C")
    vmax_c=min(vmax_c,2.0)
    norm_c=mcolors.TwoSlopeNorm(vmin=-vmax_c,vcenter=0,vmax=vmax_c)
    sc_last=None
    panel_products=[("AORC",b_aorc_c),("ERA5 Downscaled",b_ds_c),("ERA5 Native",b_e5_c),("MMM-Downscaled",b_mmm_c)]
    aorc_ax_ref=None
    for i_p,(pname,pvals) in enumerate(panel_products):
        axp=fig.add_subplot(gs_c[0,i_p])
        axp.set_facecolor(OCEAN_GREY)
        if i_p==0: aorc_ax_ref=axp
        sc_last=axp.scatter(lons_c,lats_c,c=pvals,cmap="RdBu_r",norm=norm_c,
                            s=14,edgecolor="k",lw=0.25)
        add_coastline(axp,color="black",lw=0.7,zorder=10)
        axp.set_xlim(LON0,LON1); axp.set_ylim(LAT0,LAT1); axp.set_aspect("equal")
        axp.set_title(pname,fontsize=11,fontweight="bold")
        axp.tick_params(labelsize=9)
        set_standard_lon_ticks(axp,labelsize=8)
        if i_p==0: axp.set_ylabel("Lat (\u00b0N)",fontsize=10)
        else: axp.set_yticklabels([])
        axp.set_xlabel("Lon (\u00b0W)",fontsize=10)
    if aorc_ax_ref is not None:
        aorc_bbox=aorc_ax_ref.get_position()
        # Header is centred on the figure, which is correct for the 4-panel row.
        fig.text(0.51,aorc_bbox.y1+0.030,
                 f"(c) Spatial Distribution of Per-Station Mean Bias vs GHCN-D Stations",
                 fontsize=FS_TITLE,fontweight="bold",ha="center")
    cax_c=fig.add_subplot(gs_c[1,:])
    cb_c=fig.colorbar(sc_last,cax=cax_c,orientation="horizontal")
    cb_c.set_label("Mean bias (\u00b0C)",fontsize=11); cb_c.ax.tick_params(labelsize=10)

    fig.suptitle(f"Station Validation Summary for {LONGNAME} ({OBS_YR0}-{OBS_YR1})",
                 fontsize=FS_SUPTITLE,fontweight="bold",y=0.985)
    _save(fig,f"{fn}_station_validation_summary")

def fig_combined_improvement_kge(era5ds,era5n,aorc,mask,_season_idx,mask_e5=None,fn="Fig5"):
    """Fig 5: Downscaling Added-Value 2x2.
       (a) Absolute |Bias|/RMSE bars (Native vs DS) with Delta annotations between pairs.
       (b) KGE component decomposition (r, beta, gamma): Native vs DS per season.
       (c) KGE bars (Native vs DS) with skill thresholds at 0.41 and 0.75.
       (d) KGE skill score with improvement-zone shading above 0.
    """
    from matplotlib.patches import Patch
    skeys=["Annual"]+list(SEASONS); nseas=len(skeys)
    if mask_e5 is None: mask_e5=np.ones(era5n.shape[1:],bool)
    Te5=era5n.shape[0]

    # --- Compute absolute |Bias|/RMSE per season ---
    abs_bias_nat=[]; abs_bias_ds=[]; rmse_nat=[]; rmse_ds=[]
    dbias=[]; drmse=[]
    for sk in skeys:
        if sk=="Annual":
            idx=np.ones(era5ds.shape[0],bool); idx_e5=np.ones(Te5,bool)
        else:
            idx=_season_idx[sk]; idx_e5=idx[:Te5]
        nat_dm=np.nanmean(era5n[idx_e5].reshape(idx_e5.sum(),-1)[:,mask_e5.ravel()],axis=1)
        aorc_dm_e5=np.nanmean(aorc[idx_e5].reshape(idx_e5.sum(),-1)[:,mask.ravel()],axis=1)
        ds_dm=np.nanmean(era5ds[idx].reshape(idx.sum(),-1)[:,mask.ravel()],axis=1)
        aorc_dm_ds=np.nanmean(aorc[idx].reshape(idx.sum(),-1)[:,mask.ravel()],axis=1)
        bn=np.abs(np.nanmean(nat_dm-aorc_dm_e5)); rn=np.sqrt(np.nanmean((nat_dm-aorc_dm_e5)**2))
        bd=np.abs(np.nanmean(ds_dm-aorc_dm_ds));  rd=np.sqrt(np.nanmean((ds_dm-aorc_dm_ds)**2))
        abs_bias_nat.append(bn); abs_bias_ds.append(bd)
        rmse_nat.append(rn);     rmse_ds.append(rd)
        dbias.append(bn-bd); drmse.append(rn-rd)

    # --- Compute KGE + components per season (Kling 2012 form) ---
    def _ts_kge_components(sim_ts,obs_ts):
        fin=np.isfinite(sim_ts)&np.isfinite(obs_ts); s=sim_ts[fin]; o=obs_ts[fin]
        if len(s)<100: return (np.nan,np.nan,np.nan,np.nan)
        N=float(len(s)); ms=np.mean(s); mo=np.mean(o)
        ss=np.sqrt(np.sum((s-ms)**2)/N); so=np.sqrt(np.sum((o-mo)**2)/N)
        if so==0 or ss==0 or mo==0: return (np.nan,np.nan,np.nan,np.nan)
        cov=np.sum((s-ms)*(o-mo))/N; r=cov/(ss*so)
        beta=ms/mo
        gamma=(ss/ms)/(so/mo) if ms!=0 else np.nan
        kge=1-np.sqrt((r-1)**2+(beta-1)**2+(gamma-1)**2)
        return (kge,r,beta,gamma)
    kge_nat=[]; kge_ds=[]; kge_skill=[]
    rN=[]; rD=[]; bN=[]; bD=[]; gN=[]; gD=[]
    for pk in skeys:
        if pk=="Annual":
            idx_ds=np.ones(era5ds.shape[0],bool); idx_n=np.ones(Te5,bool)
        else:
            idx_full=_season_idx[pk]; idx_ds=idx_full; idx_n=idx_full[:Te5]
        nat_dm=np.nanmean(era5n[idx_n].reshape(idx_n.sum(),-1)[:,mask_e5.ravel()],axis=1)
        aorc_dm_n=np.nanmean(aorc[:Te5][idx_n].reshape(idx_n.sum(),-1)[:,mask.ravel()],axis=1)
        ds_dm=np.nanmean(era5ds[idx_ds].reshape(idx_ds.sum(),-1)[:,mask.ravel()],axis=1)
        aorc_dm_d=np.nanmean(aorc[idx_ds].reshape(idx_ds.sum(),-1)[:,mask.ravel()],axis=1)
        kn,rn_,bn_,gn_=_ts_kge_components(nat_dm,aorc_dm_n)
        kd,rd_,bd_,gd_=_ts_kge_components(ds_dm,aorc_dm_d)
        kge_nat.append(kn); kge_ds.append(kd)
        rN.append(rn_); rD.append(rd_); bN.append(bn_); bD.append(bd_); gN.append(gn_); gD.append(gd_)
        skill=1-(1-kd)/(1-kn) if np.isfinite(kn) and np.isfinite(kd) and kn!=1 else np.nan
        kge_skill.append(skill)

    fig=plt.figure(figsize=(14,11.5))
    import matplotlib.gridspec as gridspec
    gs=gridspec.GridSpec(2,2,figure=fig,hspace=0.46,wspace=0.26,
                         left=0.07,right=0.97,top=0.87,bottom=0.06)
    x=np.arange(nseas)

    # ===== Panel (a): Absolute |Bias| / RMSE with Delta annotations =====
    ax=fig.add_subplot(gs[0,0]); w=0.20
    ax.bar(x-1.5*w,abs_bias_nat,w,color=ERA5N_C,label="|Bias| Native",alpha=.85,edgecolor="white",lw=0.4)
    ax.bar(x-0.5*w,abs_bias_ds, w,color=ERA5_C, label="|Bias| Downscaled",alpha=.85,edgecolor="white",lw=0.4)
    ax.bar(x+0.5*w,rmse_nat,    w,color=ERA5N_C,label="RMSE Native",alpha=.85,hatch="//",edgecolor="white",lw=0.4)
    ax.bar(x+1.5*w,rmse_ds,     w,color=ERA5_C, label="RMSE Downscaled",alpha=.85,hatch="//",edgecolor="white",lw=0.4)
    _all_bars=[v for v in (abs_bias_nat+abs_bias_ds+rmse_nat+rmse_ds) if np.isfinite(v)]
    ymax=(max(_all_bars) if _all_bars else 1.0)*1.26
    ax.set_ylim(0,ymax)
    for i in range(nseas):
        if np.isfinite(dbias[i]):
            y_b=max(abs_bias_nat[i],abs_bias_ds[i])+ymax*0.025
            ax.text(i-w, y_b, f"\u0394{dbias[i]:+.2f}",ha="center",va="bottom",fontsize=10.5,
                    color="#222",fontweight="bold")
        if np.isfinite(drmse[i]):
            y_r=max(rmse_nat[i],rmse_ds[i])+ymax*0.025
            ax.text(i+w, y_r, f"\u0394{drmse[i]:+.2f}",ha="center",va="bottom",fontsize=10.5,
                    color="#222",fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(skeys,fontsize=11)
    ax.set_ylabel("Error magnitude (\u00b0C)")
    ax.legend(fontsize=9,loc="upper right",framealpha=0.5,ncol=FIG5_LEGEND_NCOL,
              bbox_to_anchor=(0.995,0.995),handletextpad=0.5,labelspacing=0.35,
              columnspacing=1.0)
    ax.set_title("(a) Absolute |Bias| and RMSE results",
                 fontsize=FS_TITLE,fontweight="bold")
    ax.grid(True,axis="y",ls=":",lw=0.4,alpha=0.5)

    ax=fig.add_subplot(gs[0,1]); w2=.30
    ax.bar(x-w2/2,kge_nat,w2,color=ERA5N_C,label="ERA5 Native",alpha=.85,edgecolor="white",lw=0.4)
    ax.bar(x+w2/2,kge_ds, w2,color=ERA5_C, label="ERA5 Downscaled",alpha=.85,edgecolor="white",lw=0.4)
    for i in range(nseas):
        if np.isfinite(kge_nat[i]): ax.text(i-w2/2,kge_nat[i]+.012,f"{kge_nat[i]:.2f}",ha="center",va="bottom",fontsize=10.5)
        if np.isfinite(kge_ds[i]):  ax.text(i+w2/2,kge_ds[i]+.012, f"{kge_ds[i]:.2f}", ha="center",va="bottom",fontsize=10.5)
    ax.axhline(0.41,color="#666666",lw=1.0,ls="--",alpha=0.85,zorder=2)
    ax.axhline(0.75,color="#222222",lw=1.0,ls="-.",alpha=0.85,zorder=2)
    finite_kge=[v for v in (kge_nat+kge_ds) if np.isfinite(v)]
    if finite_kge:
        ymin_v=min(0.0,min(finite_kge)-0.1); ymax_v=max(0.95,max(finite_kge)+0.18)
    else:
        ymin_v,ymax_v=0,1.0
    ax.set_ylim(ymin_v,ymax_v)
    ax.text(0.02,0.97,
            "KGE = 0.41  inertial-mean benchmark\nKGE = 0.75  high-skill threshold",
            transform=ax.transAxes,fontsize=11,va="top",ha="left",
            bbox=dict(boxstyle="round,pad=0.35",fc="white",ec="0.6",alpha=0.55))
    ax.set_xticks(x); ax.set_xticklabels(skeys,fontsize=11); ax.set_ylabel("KGE")
    ax.legend(fontsize=10,loc="upper right",framealpha=0.35)
    ax.grid(True,axis="y",ls=":",lw=0.4,alpha=0.4)
    ax.set_title("(b) Kling-Gupta Efficiency results",
                 fontsize=FS_TITLE,fontweight="bold")

    ax=fig.add_subplot(gs[1,0])
    # 6 bars per season, all packed inside a group width of ~0.90 units
    # (x-tick spacing is 1 by construction) so groups sit clear of each
    # other. bar_w = group_w / 6, but shrink each drawn bar to 92% of a
    # slot so there is a hairline gap between adjacent bars.
    group_w=0.90
    slot_w=group_w/6.0        # width of one slot
    bar_w=slot_w*0.92         # actual drawn bar width
    # Slot centre positions relative to season centre x (six slots, half-widths
    # from -2.5 to +2.5 slot_w). Order matches the requested column layout.
    slot_off=np.array([-2.5,-1.5,-0.5,0.5,1.5,2.5])*slot_w
    comps=[("r","#0072B2"),("\u03b2","#E69F00"),("\u03b3","#009E73")]
    nat_vals=[rN,bN,gN]   # native values per component (indexed r=0, β=1, γ=2)
    ds_vals =[rD,bD,gD]   # downscaled values per component
    for ci,(clab,col) in enumerate(comps):
        # Native slot: slot 2*ci; DS slot: slot 2*ci+1.
        x_n=x+slot_off[2*ci]     # r Native (ci=0), β Native (ci=1), γ Native (ci=2)
        x_d=x+slot_off[2*ci+1]   # r DS, β DS, γ DS
        ax.bar(x_n,nat_vals[ci],bar_w,color=col,alpha=0.55,edgecolor="white",lw=0.3,zorder=3)
        ax.bar(x_d,ds_vals[ci], bar_w,color=col,alpha=0.95,edgecolor="white",lw=0.3,zorder=3)
    ax.axhline(1.0,color="k",lw=0.9,ls="--",alpha=0.75,zorder=2)
    ax.set_xticks(x); ax.set_xticklabels(skeys,fontsize=FS_TICK+1)
    ax.set_ylabel("Component value",fontsize=FS_AXLABEL)
    ax.tick_params(labelsize=FS_TICK)
    leg_handles=[]
    for clab,col in comps:
        leg_handles.append(Patch(facecolor=col,alpha=0.55,label=f"{clab} Native"))
        leg_handles.append(Patch(facecolor=col,alpha=0.95,label=f"{clab} DS"))
    ax.legend(handles=leg_handles,fontsize=FS_LEGEND-1,ncol=3,loc="lower right",framealpha=0.5)
    ax.set_title("(c) KGE component decomposition (r, \u03b2 = bias ratio, \u03b3 = variability ratio)",
                 fontsize=FS_TITLE,fontweight="bold")
    ax.grid(True,axis="y",ls=":",lw=0.4,alpha=0.4)

    ax=fig.add_subplot(gs[1,1])
    finite_skill=[v for v in kge_skill if np.isfinite(v)]
    if finite_skill:
        ymax_d=max(finite_skill)*1.20 if max(finite_skill)>0 else 0.10
        ymin_d=min(finite_skill)*1.20 if min(finite_skill)<0 else -0.05
    else:
        ymin_d,ymax_d=-0.1,0.5
    ax.axhspan(0,ymax_d,color="#009E73",alpha=0.10,zorder=0)
    if ymin_d<0:
        ax.axhspan(ymin_d,0,color="#D55E00",alpha=0.07,zorder=0)
    bar_colors=["#009E73" if (np.isfinite(s) and s>0) else "#D55E00" for s in kge_skill]
    ax.bar(x,kge_skill,0.55,color=bar_colors,alpha=.85,edgecolor="white",lw=0.4,zorder=3)
    for i,v in enumerate(kge_skill):
        if np.isfinite(v):
            offset=(ymax_d-ymin_d)*0.015
            ax.text(i,v+(offset if v>=0 else -offset),f"{v:+.2f}",
                    ha="center",va="bottom" if v>=0 else "top",fontsize=11.5,fontweight="bold")
    ax.axhline(0,color="k",lw=0.9,zorder=4)
    ax.set_ylim(ymin_d,ymax_d)
    ax.set_xticks(x); ax.set_xticklabels(skeys,fontsize=11)
    ax.set_ylabel(r"$KGE_{ss}$")
    # Corrected legend: two explicit Patches matching the bar/shading colors.
    leg_d=[Patch(facecolor="#009E73",alpha=0.85,label="DS improves skill ($KGE_{ss}>0$)"),
           Patch(facecolor="#D55E00",alpha=0.85,label="DS degrades skill ($KGE_{ss}<0$)")]
    ax.legend(handles=leg_d,fontsize=11,loc="lower right",framealpha=0.5)
    ax.grid(True,axis="y",ls=":",lw=0.4,alpha=0.4)
    ax.set_title(r"(d) KGE Skill Score ($KGE_{ss}$)",fontsize=FS_TITLE,fontweight="bold")

    fig.suptitle(f"Measuring Downscaling Added Value for {LONGNAME} ({YR0}-{YR1})\n"
                 f"by Comparing ERA5 (Native and Downscaled) Against AORC",
                 fontsize=FS_SUPTITLE,fontweight="bold",y=0.985)
    _save(fig,f"{fn}_added_value")

# ===================== METRICS TABLE & ENSEMBLE =====================


def per_pixel_sign_agreement(mc,obs_clim,mask):
    """Per-grid sign-agreement: at each pixel, fraction of ensemble members
    whose bias-against-AORC has the same sign as the majority. Returns 2D array
    in percent (50-100); NaN outside mask.
    """
    biases=[(mc[m]-obs_clim) for m in mc]
    if not biases: return np.full_like(obs_clim,np.nan,dtype=float)
    stack=np.stack(biases,axis=0)  # (N,lat,lon)
    pos=np.sum(stack>0,axis=0); neg=np.sum(stack<0,axis=0)
    agree=100.0*np.maximum(pos,neg)/max(1,len(biases))
    out=np.where(mask,agree,np.nan)
    return out

def fig_ensemble_diag(spread_abs,boot_cis,mnames,
                      model_metrics=None,sign_agree_map=None,grid_lat=None,
                      fn="Fig4"):
    """Fig 4 (manuscript order): Ensemble Diagnostics 2x2.
    Bottom row uses unequal column widths: (c) is compact on the left; (d) is
    wider on the right so the sign-agreement map fills available space without
    empty margins."""
    import matplotlib.gridspec as gridspec
    fig=plt.figure(figsize=(14,10))
    gs=gridspec.GridSpec(2,1,figure=fig,height_ratios=[1.0,1.10],hspace=0.42,
                         left=0.07,right=0.97,top=0.93,bottom=0.07)
    gs_top=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[0],wspace=0.28)
    gs_bot=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1],wspace=0.22,width_ratios=[1.0,1.65])

    n=len(mnames)
    # Pre-compute model bias means + colors so they're shared by (a) and (c).
    model_biases=[]
    for m in mnames:
        lo,hi=boot_cis.get(m,(np.nan,np.nan))
        model_biases.append((lo+hi)/2)
    model_biases=np.array(model_biases)
    model_colors=[MCMAP(i/max(n-1,1)) for i in range(n)]

    # ----- Panel (a): Bootstrap CI on MBE [top-left] -----
    ax=fig.add_subplot(gs_top[0,0])
    for i,m in enumerate(mnames):
        lo,hi=boot_cis.get(m,(np.nan,np.nan)); mid=(lo+hi)/2
        ax.errorbar(i,mid,yerr=[[mid-lo],[hi-mid]],fmt="o",ms=5,color=model_colors[i],capsize=3)
    ax.axhline(0,color=".4",lw=.5,ls="--"); ax.set_xticks(range(n))
    ax.set_xticklabels(mnames,rotation=45,ha="right",fontsize=10.5); ax.set_ylabel("MBE (\u00b0C)")
    ax.set_title("(a) 95% Bootstrap CI on MBE",fontsize=FS_TITLE,fontweight="bold")

    # ----- Panel (b): Extreme-index bias bars [top-right] (formerly c) -----
    ax=fig.add_subplot(gs_top[0,1])
    if model_metrics:
        tx=[model_metrics.get(m,{}).get("TX90p_bias",np.nan) for m in mnames]
        tn=[model_metrics.get(m,{}).get("TN10p_bias",np.nan) for m in mnames]
        x=np.arange(n); w=0.4
        ax.barh(x-w/2,tx,w,color="#D55E00",label="TX90p bias (hot-day freq)",alpha=0.85)
        ax.barh(x+w/2,tn,w,color="#0072B2",label="TN10p bias (cold-night freq)",alpha=0.85)
        ax.axvline(0,color="k",lw=0.6)
        ax.set_yticks(x); ax.set_yticklabels(mnames,fontsize=10.5)
        ax.set_xlabel("Bias (percentage points)")
        ax.legend(fontsize=11,loc="lower left",framealpha=0.5)
        ax.grid(True,axis="x",ls=":",lw=.4,alpha=.5)
    else:
        ax.text(0.5,0.5,"TX90p/TN10p metrics not provided",ha="center",va="center",
                transform=ax.transAxes,fontsize=12,color="gray")
    ax.set_title("(b) Extreme-index bias by model",fontsize=FS_TITLE,fontweight="bold")

    ax=fig.add_subplot(gs_bot[0,0])
    x_jitter=np.random.default_rng(42).uniform(-0.15,0.15,n)
    dot_handles=[]
    for i,m in enumerate(mnames):
        h,=ax.plot(x_jitter[i],model_biases[i],"o",ms=8,color=model_colors[i],
                   markeredgecolor="white",markeredgewidth=0.5,zorder=5,label=m)
        dot_handles.append(h)
    ens_mean=np.nanmean(model_biases); ens_std=np.nanstd(model_biases)
    ax.axhline(ens_mean,color=MMM_C,lw=1.5,ls="-",zorder=3,label=f"Ens. mean = {ens_mean:.2f} \u00b0C")
    ax.axhspan(ens_mean-ens_std,ens_mean+ens_std,color=MMM_C,alpha=0.15,zorder=2,
               label=f"\u00b11 SD = {spread_abs:.2f} \u00b0C")
    ax.axhline(0,color=".4",lw=.5,ls="--",zorder=1)
    ax.set_xlim(-0.5,0.5); ax.set_xticks([]); ax.set_ylabel("MBE (\u00b0C)")
    # Model legend OUTSIDE the plot on the right so dots aren't obscured
    ax.legend(handles=dot_handles,loc="center left",bbox_to_anchor=(1.02,0.5),
              fontsize=10.5,framealpha=0.5,title="Model",title_fontsize=11)
    ax.set_title(f"(c) Inter-model spread ({spread_abs:.2f} \u00b0C)",fontsize=FS_TITLE,fontweight="bold")

    # ----- Panel (d): Per-grid sign-agreement map (gets the wider right column) -----
    ax=fig.add_subplot(gs_bot[0,1])
    ax.set_facecolor(OCEAN_GREY)
    if sign_agree_map is not None and grid_lat is not None:
        im=ax.imshow(sign_agree_map,extent=[LON0,LON1,LAT0,LAT1],
                     origin="lower" if grid_lat[0]<grid_lat[-1] else "upper",
                     cmap="YlOrRd",vmin=50,vmax=100,aspect="auto",interpolation="nearest")
        cb=fig.colorbar(im,ax=ax,shrink=0.90,pad=0.02,label="Sign agreement (%)")
        cb.ax.tick_params(labelsize=10)
    ax.set_xlim(LON0,LON1); ax.set_ylim(LAT0,LAT1); ax.set_aspect("equal")
    set_standard_lon_ticks(ax)
    ax.tick_params(labelsize=10)
    ax.set_xlabel("Lon (\u00b0W)",fontsize=11); ax.set_ylabel("Lat (\u00b0N)",fontsize=11)
    ax.set_title("(d) Per-grid ensemble sign agreement\n(% of members with same-sign bias)",
                 fontsize=FS_TITLE,fontweight="bold")

    fig.suptitle(f"Ensemble Diagnostics for {LONGNAME} ({YR0}-{YR1})",
                 fontsize=FS_SUPTITLE,fontweight="bold",y=0.985)
    _save(fig,f"{fn}_ensemble_diag")

# ===================== MAIN =====================
def main():
    import time as _time; _t0=_time.time()
    def _elapsed(): return f"[{(_time.time()-_t0)/60:.1f} min]"
    print("="*72); print("  tas Spatial Performance v1.0 (manuscript-ordered)"); print("="*72)

    print(f"\n{_elapsed()} [1] AORC ..."); aorc,alat,alon,aorc_time=load_nc(AORC_PATH)
    aorc=maybe_k_to_c(aorc,"AORC"); aorc=physical_qc(aorc,"AORC"); aorc,alat,alon=ensure_lat_descending(aorc,alat,alon,"AORC")
    aorc=aorc.astype(np.float32); print(f"    shape={aorc.shape}")
    print(f"{_elapsed()} [2] ERA5 native ..."); era5n,e5lat,e5lon,era5n_time=load_nc(ERA5_PATH)
    era5n=maybe_jhr_to_wm2(era5n,"ERA5-nat")
    era5n=maybe_k_to_c(era5n,"ERA5-nat"); era5n=physical_qc(era5n,"ERA5-nat"); era5n,e5lat,e5lon=ensure_lat_descending(era5n,e5lat,e5lon,"ERA5-nat")
    era5n=era5n.astype(np.float32); print(f"    shape={era5n.shape}")
    print(f"{_elapsed()} [3] ERA5 downscaled ..."); era5ds,_,_,era5ds_time=load_nc(ERA5_DS_PATH)
    era5ds=maybe_jhr_to_wm2(era5ds,"ERA5 Downscaled")
    era5ds=maybe_k_to_c(era5ds,"ERA5 Downscaled"); era5ds=physical_qc(era5ds,"ERA5 Downscaled"); era5ds=era5ds.astype(np.float32); print(f"    shape={era5ds.shape}")
    for _nm,_tc,_dat in [("AORC",aorc_time,aorc),("ERA5-nat",era5n_time,era5n),("ERA5 Downscaled",era5ds_time,era5ds)]:
        if _tc is not None and len(_tc):
            print(f"    [time-axis] {_nm}: real time {pd.Timestamp(_tc[0]):%Y-%m-%d} -> "
                  f"{pd.Timestamp(_tc[-1]):%Y-%m-%d} ({len(_tc)} steps, data has {_dat.shape[0]})")
        else:
            print(f"    [time-axis] {_nm}: NO embedded time coord (will synthesize from {YR0})")

    daymet_clim=None; daymet_lat=None; daymet_lon=None
    if DAYMET_PATH.exists():
        print(f"{_elapsed()} [3b] Daymet ...")
        try:
            daymet_raw,daymet_lat_orig,daymet_lon_orig,_=load_nc(DAYMET_PATH)
            daymet_raw=maybe_k_to_c(daymet_raw,"Daymet")
            daymet_raw=physical_qc(daymet_raw,"Daymet")
            daymet_clim_orig=np.nanmean(daymet_raw,0).astype(np.float32)
            print(f"    Daymet raw shape={daymet_raw.shape}, clim={daymet_clim_orig.shape}")
            del daymet_raw; gc.collect()
            daymet_lat=daymet_lat_orig; daymet_lon=daymet_lon_orig; daymet_clim=daymet_clim_orig
        except Exception as e: print(f"    Daymet load failed: {e}")
    else: print(f"  Daymet file not found: {DAYMET_PATH}")

    nex_native_clim={}; nex_native_monthly={}; nex_native_lat=None; nex_native_lon=None
    import glob as _glob
    nex_files=sorted(_glob.glob(NEX_NATIVE_GLOB))
    if nex_files:
        print(f"{_elapsed()} [3c] Native NEX-GDDP ({len(nex_files)} files) ...")
        for nf_path in nex_files:
            nf=Path(nf_path)
            try:
                parts=nf.stem.split("_"); nm="_".join(parts[2:]) if len(parts)>2 else nf.stem
                for tok in ["_historical","_ssp245","_ssp585","_1950","_1980","_r1i1p1"]:
                    if tok in nm: nm=nm[:nm.index(tok)]
                # Skip ERA5 / reanalysis files — not CMIP6 native models and often on a
                # different native grid (e.g. 22x16 vs 21x16), which breaks the stack.
                if "ERA5" in nm.upper():
                    print(f"    {nm}: skipped (ERA5 is not a CMIP6 native model)"); continue
                nd,nlat_,nlon_,nt_=load_nc(str(nf)); nd=maybe_k_to_c(nd,f"Native-{nm}"); nd=physical_qc(nd,f"Native-{nm}")
                nex_native_clim[nm]=np.nanmean(nd,0).astype(np.float32)
                if nex_native_lat is None: nex_native_lat=nlat_; nex_native_lon=nlon_
                nmask=np.ones(nd.shape[1:],bool); nT=nd.shape[0]
                ndates=nt_ if nt_ is not None else pd.date_range(f"{YR0}-01-01",periods=nT,freq="D")
                nex_native_monthly[nm]=compute_monthly_domain_mean(nd,nmask,ndates)
                del nd; gc.collect(); print(f"    {nm}: loaded ({nex_native_clim[nm].shape})")
            except Exception as e: print(f"    {nf.name}: ERR {e}")
        print(f"  {len(nex_native_clim)} native models loaded")
    else: print(f"  No native NEX-GDDP files found matching: {NEX_NATIVE_GLOB}")

    print(f"{_elapsed()} [4] CMIP6 models (Zarr) ...")
    mf=discover_zarr_models(ZARR_DIR)
    if not mf: print("  ERROR: No models found in Zarr dir. Exiting."); return

    print(f"{_elapsed()}   Peeking grid shape ...")
    mn0=list(mf.keys())[0]; d0,_,_,_,_=load_zarr_model(mf[mn0])
    NLAT,NLON=d0.shape[1:]; T_model=d0.shape[0]; del d0; gc.collect()
    lat=np.linspace(LAT1,LAT0,NLAT); lon=np.linspace(LON0,LON1,NLON)
    mask=np.ones((NLAT,NLON),bool); print(f"  Grid: {NLAT}x{NLON}")

    if daymet_clim is not None and daymet_clim.shape!=(NLAT,NLON):
        print(f"  Regridding Daymet {daymet_clim.shape} -> ({NLAT},{NLON}) ...")
        if daymet_lat is not None and daymet_lon is not None:
            try:
                from scipy.interpolate import RegularGridInterpolator
                dl=daymet_lat; dd=daymet_clim
                if dl[0]>dl[-1]: dl=dl[::-1]; dd=dd[::-1,:]
                interp=RegularGridInterpolator((dl,daymet_lon),dd,method="nearest",bounds_error=False,fill_value=np.nan)
                tgt_lat=np.linspace(LAT1,LAT0,NLAT); tgt_lon=np.linspace(LON0,LON1,NLON)
                grid_lat,grid_lon=np.meshgrid(tgt_lat,tgt_lon,indexing="ij")
                daymet_clim=interp((grid_lat,grid_lon)).astype(np.float32)
                daymet_lat=tgt_lat; daymet_lon=tgt_lon
            except Exception as e:
                daymet_clim=crop_to_common(NLAT,NLON,daymet_clim)
                daymet_lat=np.linspace(LAT1,LAT0,NLAT); daymet_lon=np.linspace(LON0,LON1,NLON)
        else:
            daymet_clim=crop_to_common(NLAT,NLON,daymet_clim)
            daymet_lat=np.linspace(LAT1,LAT0,NLAT); daymet_lon=np.linspace(LON0,LON1,NLON)

    aorc=crop_to_common(NLAT,NLON,aorc); era5ds=crop_to_common(NLAT,NLON,era5ds)
    T=min(aorc.shape[0],era5ds.shape[0],T_model)
    print(f"  Common T={T}"); aorc=aorc[:T]; era5ds=era5ds[:T]
    Te5=min(era5n.shape[0],T); era5n=era5n[:Te5]
    def _dates_for(real_time,n,name):
        if real_time is not None and len(real_time)>=n:
            d=pd.DatetimeIndex(pd.to_datetime(real_time[:n]))
            print(f"    [dates] {name}: using embedded time {d[0]:%Y-%m} -> {d[-1]:%Y-%m} (n={n})")
            return d
        d=pd.date_range(f"{YR0}-01-01",periods=n,freq="D")
        print(f"    [dates] {name}: SYNTHESIZED {d[0]:%Y-%m} -> {d[-1]:%Y-%m} (n={n})"
              + ("" if real_time is None else f"  [WARN embedded time too short: {len(real_time)}<{n}]"))
        return d
    gdates=_dates_for(aorc_time,T,"gdates(AORC/DS/MMM)")
    e5dates=_dates_for(era5n_time,Te5,"e5dates(ERA5-nat)")
    # ==================== DETAILED TIME-AXIS LOG ====================
    print("\n" + "="*72)
    print("  TIME-AXIS DIAGNOSTIC LOG (v1.0)")
    print("="*72)
    print(f"  Config: YR0={YR0}, YR1={YR1}  (expected span {YR0}-01 .. {YR1}-12 = {(YR1-YR0+1)*12} months)")
    for _nm,_tc,_dat in [("AORC",aorc_time,aorc),("ERA5-nat",era5n_time,era5n),("ERA5 Downscaled",era5ds_time,era5ds)]:
        if _tc is not None and len(_tc):
            _t=pd.DatetimeIndex(pd.to_datetime(_tc))
            _yrs=np.unique(_t.year)
            print(f"  [{_nm}] embedded time: {_t[0]:%Y-%m-%d} .. {_t[-1]:%Y-%m-%d} | "
                  f"{len(_t)} steps | data rows={_dat.shape[0]} | "
                  f"years {_yrs.min()}-{_yrs.max()} | step0->1 = {(_t[1]-_t[0]).days if len(_t)>1 else '?'}d")
        else:
            print(f"  [{_nm}] NO embedded time -> SYNTHESIZED from {YR0}; data rows={_dat.shape[0]}")
    print(f"  [gdates chosen]  {gdates[0]:%Y-%m-%d} .. {gdates[-1]:%Y-%m-%d} ({len(gdates)} days, "
          f"{gdates.year.min()}-{gdates.year.max()})")
    print(f"  [e5dates chosen] {e5dates[0]:%Y-%m-%d} .. {e5dates[-1]:%Y-%m-%d} ({len(e5dates)} days)")
    print("="*72 + "\n")
    # ======================================================================
    _season_idx={}
    _gmonths=pd.DatetimeIndex(gdates).month.values
    for sk,si in SEASONS.items(): _season_idx[sk]=np.isin(_gmonths,si["months"])

    MF=list(mf.keys()); n_mod=len(MF)
    print(f"\n{_elapsed()}   Processing {n_mod} models ...")
    mmm_sum=np.zeros((T,NLAT,NLON),dtype=np.float64)
    mmm_cnt=np.zeros((T,NLAT,NLON),dtype=np.int16)
    m_clim={}; m_taylor={}; m_metrics={}; m_boot={}; m_se_clim={}
    ac_ann=np.nanmean(aorc,0); csv1=[]
    # AORC seasonal climatologies are identical for every model, so compute
    # them ONCE here rather than re-reducing the full AORC cube inside the
    # per-model loop (was 4 seasonal passes x N models).
    ac_seas={sk:np.nanmean(aorc[_season_idx[sk]],0) for sk in SEASONS}
    for mi,mn in enumerate(MF):
        print(f"    {_elapsed()} [{mi+1}/{n_mod}] {mn} ...",end=" ",flush=True)
        try:
            d,_,_,_,se_d=load_zarr_model(mf[mn]); d=maybe_k_to_c(d,mn); d=physical_qc(d,mn); d=crop_to_common(NLAT,NLON,d)[:T]; d=d.astype(np.float32)
            if se_d is not None:
                if se_d.ndim==3: se_d=crop_to_common(NLAT,NLON,se_d)[:T]
                m_se_clim[mn]=np.nanmean(se_d,0).astype(np.float32) if se_d.ndim==3 else se_d.astype(np.float32)
        except Exception as e: print(f"ERR {e}"); continue
        if np.all(np.isnan(d)) or d.size==0: print(f"SKIP (no valid data)"); continue
        v=np.isfinite(d); mmm_sum=np.where(v,mmm_sum+d,mmm_sum); mmm_cnt+=v.astype(np.int16)
        m_clim[mn]=np.nanmean(d,0)
        m_taylor[mn]={}; ns,r,_=fld_stats(m_clim[mn],ac_ann,mask); m_taylor[mn]["Annual"]=(ns,r)
        csv1.append([mn,"Annual",f"{ns:.2f}",f"{r:.2f}","",""])
        for sk in SEASONS:
            ns,r,_=fld_stats(np.nanmean(d[_season_idx[sk]],0),ac_seas[sk],mask); m_taylor[mn][sk]=(ns,r)
        m_metrics[mn]=domain_metrics(d,aorc,mask)
        lo,hi=bootstrap_ci(d,aorc,mask,1000); m_boot[mn]=(lo,hi)
        del d; gc.collect(); print("done")

    MF=[m for m in MF if m in m_clim]
    print(f"  {len(MF)} models successfully loaded")

    MMM_DROP_OUTLIERS=True
    MMM_MAX_DROP_FRAC=0.5
    outliers=[]
    if len(MF)>=4:
        r_per_model={mn:m_taylor[mn]["Annual"][1] for mn in MF if np.isfinite(m_taylor[mn]["Annual"][1])}
        if r_per_model:
            r_vals=np.array(list(r_per_model.values()))
            med=float(np.median(r_vals)); q1,q3=np.percentile(r_vals,[25,75]); iqr=q3-q1
            if med<0.5:
                lo_thr=med-2.0*iqr
                print(f"  [MMM anomaly check] median r ({med:.3f}) is below 0.5; "
                      f"using pure IQR-based threshold (no 0.5 floor).")
            else:
                lo_thr=max(0.5,med-2.0*iqr)
            outliers=[mn for mn,r in r_per_model.items() if r<lo_thr]
            print(f"  [MMM anomaly check] per-model r vs AORC: median={med:.3f}, "
                  f"IQR=[{q1:.3f},{q3:.3f}], threshold for outlier flag={lo_thr:.3f}")
            for mn in sorted(r_per_model,key=lambda k:r_per_model[k]):
                marker=" <-- OUTLIER (excluded from MMM)" if mn in outliers else ""
                print(f"    {mn:20s}  r={r_per_model[mn]:+.3f}{marker}")
            # Safeguard: refuse to drop more than MMM_MAX_DROP_FRAC of the ensemble
            if len(outliers)>MMM_MAX_DROP_FRAC*len(r_per_model):
                print(f"  *** WARNING: outlier rule flagged {len(outliers)} of "
                      f"{len(r_per_model)} models (>{MMM_MAX_DROP_FRAC*100:.0f}%). "
                      f"Refusing to drop to preserve a meaningful MMM. "
                      f"Keeping ALL models for the ensemble mean.")
                outliers=[]
            if outliers and MMM_DROP_OUTLIERS:
                print(f"  *** {len(outliers)} model(s) excluded from MMM: {outliers}")
                print(f"  *** Reloading these to subtract their contribution from mmm_sum / mmm_cnt ...")
                # Subtract each outlier from the running sum so MMM reflects only
                # well-behaved models. Reload each outlier's data, then zero its
                # contribution. This is more memory-efficient than rebuilding
                # mmm_sum from scratch.
                for mn_out in outliers:
                    try:
                        d_out,_,_,_,_=load_zarr_model(mf[mn_out])
                        d_out=maybe_k_to_c(d_out,mn_out)
                        d_out=physical_qc(d_out,mn_out)
                        d_out=crop_to_common(NLAT,NLON,d_out)[:T].astype(np.float32)
                        v_out=np.isfinite(d_out)
                        mmm_sum=np.where(v_out,mmm_sum-d_out,mmm_sum)
                        mmm_cnt-=v_out.astype(np.int16)
                        print(f"    {mn_out} contribution removed from MMM")
                        del d_out,v_out; gc.collect()
                    except Exception as e:
                        print(f"    {mn_out} reload for MMM removal failed: {e}")

    print(f"\n{_elapsed()}   Computing MMM ..."); no_data=(mmm_cnt==0)
    mmm_cnt[no_data]=1; mmm=(mmm_sum/mmm_cnt.astype(np.float64)).astype(np.float32); mmm[no_data]=np.nan
    del mmm_sum,mmm_cnt; gc.collect()
    _mmm_finite=np.count_nonzero(np.isfinite(mmm)); _mmm_total=mmm.size
    _mmm_pct=100.0*_mmm_finite/max(_mmm_total,1)
    print(f"    [MMM] {_mmm_finite:,} of {_mmm_total:,} cells finite ({_mmm_pct:.2f}%); "
          f"mean={np.nanmean(mmm):.3f}, contributing models in pool={len(MF)-len(outliers)}")
    if _mmm_finite==0:
        print(f"    *** WARNING: MMM is entirely NaN! Every downstream MMM figure "
              f"will be blank. Likely causes: (a) all models dropped; (b) QC clipped "
              f"every model; (c) ZARR_VAR='{ZARR_VAR}' does not match Zarr var name.")

    e5n_lat_=e5lat if e5lat is not None else np.linspace(LAT1,LAT0,era5n.shape[1])
    e5n_lon_=e5lon if e5lon is not None else np.linspace(LON0,LON1,era5n.shape[2])
    e5n_shape=(len(e5n_lat_),len(e5n_lon_))

    mmm_native_clim=None
    if nex_native_clim and nex_native_lat is not None:
        from scipy.interpolate import RegularGridInterpolator
        shapes=set(v.shape for v in nex_native_clim.values())
        if len(shapes)>1 or list(shapes)[0]!=e5n_shape:
            print(f"    [native] grid shapes present: {sorted(shapes)} -> regridding all to ERA5-native {e5n_shape}")
        for nm in list(nex_native_clim.keys()):
            nc_=nex_native_clim[nm]
            if nc_.shape==e5n_shape: continue
            try:
                nl_=np.asarray(nex_native_lat,dtype=float); dd=nc_
                if nl_[0]>nl_[-1]: nl_=nl_[::-1]; dd=dd[::-1,:]
                no_=np.asarray(nex_native_lon,dtype=float)
                if no_[0]>no_[-1]: no_=no_[::-1]; dd=dd[:,::-1]
                el_=np.asarray(e5n_lat_,dtype=float).copy()
                flip_back=el_[0]>el_[-1]
                if flip_back: el_=el_[::-1]
                eo_=np.asarray(e5n_lon_,dtype=float).copy()
                flip_lon=eo_[0]>eo_[-1]
                if flip_lon: eo_=eo_[::-1]
                interp=RegularGridInterpolator((nl_,no_),dd,method="nearest",bounds_error=False,fill_value=np.nan)
                grid_la,grid_lo=np.meshgrid(el_,eo_,indexing="ij")
                regridded=interp((grid_la,grid_lo)).astype(np.float32)
                if flip_back: regridded=regridded[::-1,:]
                if flip_lon: regridded=regridded[:,::-1]
                nex_native_clim[nm]=regridded
            except Exception as e:
                print(f"    {nm} regrid failed ({e}); dropping from native MMM")
                del nex_native_clim[nm]
        nex_native_lat=e5n_lat_; nex_native_lon=e5n_lon_
        good=[v for v in nex_native_clim.values() if v.shape==e5n_shape]
        if good:
            mmm_native_clim=np.nanmean(np.stack(good),axis=0)
            print(f"    [native] MMM native climatology from {len(good)} models, shape={mmm_native_clim.shape}")
        else:
            print("    [native] WARNING: no native models share the ERA5 grid; MMM native unavailable")

    mmm_native_monthly=None
    if nex_native_monthly:
        try:
            _frames=[]
            for nm,(d,v) in nex_native_monthly.items():
                _idx=pd.to_datetime(d).to_period("M").to_timestamp()
                _s=pd.Series(np.asarray(v,dtype=float),index=_idx,name=nm)
                _s=_s[~_s.index.duplicated(keep="first")].sort_index()
                _frames.append(_s)
            if _frames:
                _df=pd.concat(_frames,axis=1)  # outer join on month index
                _shapes=set(len(f) for f in _frames)
                if len(_shapes)>1:
                    print(f"    [native-monthly] aligning {len(_frames)} models with "
                          f"different month-counts {sorted(_shapes)} via date-keyed reindex")
                _mmm_vals=np.nanmean(_df.values,axis=1)
                mmm_native_monthly=(_df.index,_mmm_vals)
                print(f"    [native-monthly] MMM native monthly series: "
                      f"{len(_df)} months, {_df.index[0]:%Y-%m} -> {_df.index[-1]:%Y-%m}")
        except Exception as _e:
            print(f"    [native-monthly] alignment failed ({_e}); mmm_native_monthly = None")
            mmm_native_monthly=None

    print(f"{_elapsed()}   Loading SE fields ...")
    era5ds_se=None
    if ERA5_DS_SE_PATH.exists():
        try:
            se_raw,_,_,_=load_nc(ERA5_DS_SE_PATH)
            se_raw=maybe_jhr_to_wm2(se_raw,"ERA5 Downscaled SE")
            se_raw=crop_to_common(NLAT,NLON,se_raw) if se_raw.ndim==3 else se_raw
            era5ds_se=np.nanmean(se_raw,0).astype(np.float32) if se_raw.ndim==3 else se_raw.astype(np.float32)
        except Exception as e: print(f"    ERA5 Downscaled SE load failed: {e}")
    if m_se_clim:
        se_stack=np.stack([m_se_clim[m] for m in m_se_clim])
        mmm_se=np.sqrt(np.nanmean(se_stack**2,axis=0))
    else:
        mmm_se=np.nanstd(np.stack([m_clim[m] for m in MF]),axis=0)/np.sqrt(len(MF))
    aorc_se=np.nanstd(aorc,axis=0)/np.sqrt(aorc.shape[0])

    print(f"\n{_elapsed()} [SI S1.1 — Daily Pearson r maps]")
    products=[(era5ds,"ERA5 downscaled","era5ds"),(mmm,"MMM-Downscaled","mmm")]
    daily_pairs=[]
    for sd,sl,tag in products:
        rd,_=pix_r_pval(sd,aorc,mask); daily_pairs.append((sl,rd))
    fig_corr_maps_combined(daily_pairs,lat,lon,mask,"Daily","SI_S1.1_daily")


    print(f"{_elapsed()} [SI S1.2b — DS-only pattern matrix]")
    # Time-mean climatologies are reused by the pattern-correlation matrix,
    # the Taylor statistics and the Fig 1 panels; compute each once.
    mmm_ann=np.nanmean(mmm,0); e5ds_ann=np.nanmean(era5ds,0); e5n_ann=np.nanmean(era5n,0)
    cl_ds=dict(m_clim); cl_ds["MMM-Downscaled"]=mmm_ann; cl_ds["ERA5 Downscaled"]=e5ds_ann; cl_ds["AORC"]=ac_ann
    C_ds,lb_ds=sp_corr(cl_ds,mask); fig_corrmat(C_ds,lb_ds,fn="SI_S1.2b"); del cl_ds

    print(f"{_elapsed()} [Fig 6 — Seasonal Taylor diagrams]")
    mmm_t={"Annual":fld_stats(mmm_ann,ac_ann,mask)[:2]}
    era5_t={"Annual":fld_stats(e5ds_ann,ac_ann,mask)[:2]}
    for sk in SEASONS:
        idx=_season_idx[sk]
        mmm_t[sk]=fld_stats(np.nanmean(mmm[idx],0),ac_seas[sk],mask)[:2]
        era5_t[sk]=fld_stats(np.nanmean(era5ds[idx],0),ac_seas[sk],mask)[:2]
    mt_ann={m:m_taylor[m]["Annual"] for m in MF if "Annual" in m_taylor.get(m,{})}
    fig_taylor(mt_ann,mmm_t["Annual"],era5_t["Annual"],"Annual","4a")
    all_ms={sk:{m:m_taylor[m][sk] for m in MF if sk in m_taylor.get(m,{})} for sk in SEASONS}
    fig_taylor_seasonal(all_ms,mmm_t,era5_t)

    print(f"{_elapsed()}   Domain-mean monthly series (input to Fig 3b) ...")
    monthly={"AORC":compute_monthly_domain_mean(aorc,mask,gdates),
             "ERA5 downscaled":compute_monthly_domain_mean(era5ds,mask,gdates),
             "MMM-Downscaled":compute_monthly_domain_mean(mmm,mask,gdates),
             "ERA5 Native":compute_monthly_domain_mean(era5n,np.ones(era5n.shape[1:],bool),e5dates)}
    if mmm_native_monthly is not None: monthly["MMM Native"]=mmm_native_monthly

    print(f"{_elapsed()} [Fig 1 — Climatology maps]")
    e5p_lat=e5lat if e5lat is not None else np.linspace(LAT1,LAT0,era5n.shape[1])
    e5p_lon=e5lon if e5lon is not None else np.linspace(LON0,LON1,era5n.shape[2])
    nex_lat_=nex_native_lat if nex_native_lat is not None else e5p_lat
    nex_lon_=nex_native_lon if nex_native_lon is not None else e5p_lon
    panels_11=[("ERA5 Native",e5n_ann,e5p_lat,e5p_lon,True),
               ("ERA5 Downscaled",e5ds_ann,lat,lon,False),
               ("MMM-Native",mmm_native_clim if mmm_native_clim is not None else e5n_ann,nex_lat_,nex_lon_,True),
               ("MMM-Downscaled",mmm_ann,lat,lon,False),
               ("AORC",ac_ann,lat,lon,False)]
    _daymet_has_data = (daymet_clim is not None) and np.any(np.isfinite(daymet_clim))
    if _daymet_has_data or DAYMET_ALWAYS_KEEP:
        panels_11.append(("Daymet",
                          daymet_clim if daymet_clim is not None else np.full((NLAT,NLON),np.nan),
                          daymet_lat if daymet_lat is not None else lat,
                          daymet_lon if daymet_lon is not None else lon,False))
    else:
        print("  [Fig 1] Daymet panel OMITTED (no Daymet data found for this variable).")
    fig_climatology_maps(panels_11,lat,lon,mask)

    print(f"{_elapsed()}   Per-model domain metrics (input to Fig 4b + CSV) ...")
    mmm_metrics=domain_metrics(mmm,aorc,mask)
    m_metrics_all=dict(m_metrics); m_metrics_all["MMM-Downscaled"]=mmm_metrics
    with open(OUT_DIR/"Table_domain_metrics.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["Model","MBE","PBIAS","RMSE","NSE","KGE","SD_ratio","TX90p_bias","TN10p_bias"])
        for m in m_metrics_all: d=m_metrics_all[m]; w.writerow([m]+[f"{d.get(c,np.nan):.2f}" for c in ["MBE","PBIAS","RMSE","NSE","KGE","SD_ratio","TX90p_bias","TN10p_bias"]])

    print(f"{_elapsed()} [Fig 4 — Ensemble diagnostics (2x2)]")
    spread_abs,_=ensemble_spread(m_clim,ac_ann,mask)
    lo_mmm,hi_mmm=bootstrap_ci(mmm,aorc,mask,1000); m_boot["MMM-Downscaled"]=(lo_mmm,hi_mmm)
    sign_map=per_pixel_sign_agreement(m_clim,ac_ann,mask)
    fig_ensemble_diag(spread_abs,m_boot,MF+["MMM-Downscaled"],
                      model_metrics=m_metrics_all,sign_agree_map=sign_map,grid_lat=lat)

    print(f"{_elapsed()} [Fig 5 — Added-value 2x2]"); fig_combined_improvement_kge(era5ds,era5n,aorc,mask,_season_idx)
    print(f"{_elapsed()} [Fig A — Kriging SE]")
    se_panels=[]
    if era5ds_se is not None: se_panels.append(("ERA5 Downscaled",era5ds_se,ERA5_C))
    se_panels.append(("MMM-Downscaled",mmm_se,MMM_C)); fig_se_spatial_hist(se_panels,lat,lon,mask)

    if HAS_STATION_OBS:
        print(f"{_elapsed()} [SI S1.7 — Station network map]"); obs_df=load_obs(OBS_TMAX_CSV,OBS_TMIN_CSV,STN_META)
        used=obs_df[["id","latitude","longitude","name"]].drop_duplicates("id").reset_index(drop=True)
        stn_num={row["id"]:i+1 for i,row in used.iterrows()}
        stn_csv=OUT_DIR/"station_metadata.csv"
        used_out=used.copy(); used_out.insert(0,"number",[stn_num[x] for x in used["id"]])
        used_out.to_csv(stn_csv,index=False)
        lulc_d=lulc_la=lulc_lo=None
        if LULC_PATH.exists():
            try:
                lulc_d,lulc_la,lulc_lo=load_worldcover(LULC_PATH)
                if lulc_la is not None and len(lulc_la)>1 and lulc_la[0]<lulc_la[-1]: lulc_la=lulc_la[::-1]; lulc_d=lulc_d[::-1,:]
            except Exception as e: print(f"  WorldCover load failed ({e})"); lulc_d=lulc_la=lulc_lo=None
        fig_station_map(used,lat,lon,mask,lulc_data=lulc_d,lulc_lat=lulc_la,lulc_lon=lulc_lo)

        print(f"{_elapsed()} [Fig 2 — Station validation boxplots; SI station scatter + seasonal cycle]")
        obs_a=extract_at_stations(aorc,lat,lon,obs_df,gdates)
        obs_d=extract_at_stations(era5ds,lat,lon,obs_df,gdates)
        e5n_la=e5lat if e5lat is not None else np.linspace(LAT1,LAT0,era5n.shape[1])
        e5n_lo=e5lon if e5lon is not None else np.linspace(LON0,LON1,era5n.shape[2])
        obs_e5=extract_at_stations(era5n,e5n_la,e5n_lo,obs_df,e5dates)
        obs_mmm=extract_at_stations(mmm,lat,lon,obs_df,gdates)
        obs_df["aorc_val"]=obs_a["grid_val"].values; obs_df["ds_val"]=obs_d["grid_val"].values
        obs_df["e5nat_val"]=obs_e5["grid_val"].values; obs_df["mmm_val"]=obs_mmm["grid_val"].values

        # MMM native station extraction
        if nex_native_monthly:
            nex_files_reload=sorted(__import__("glob").glob(NEX_NATIVE_GLOB))
            if nex_files_reload:
                try:
                    time_lengths=[]
                    for nf_path in nex_files_reload:
                        if "ERA5" in Path(nf_path).stem.upper(): continue  # not a CMIP6 native model
                        try: nd_peek,_,_,_=load_nc(str(nf_path)); time_lengths.append(nd_peek.shape[0]); del nd_peek; gc.collect()
                        except: pass
                    if not time_lengths: raise ValueError("No native files readable")
                    T_common=min(time_lengths); nat_sum=None; nat_cnt=0; nat_lat_r=None; nat_lon_r=None
                    nat_shape=None
                    for nf_path in nex_files_reload:
                        if "ERA5" in Path(nf_path).stem.upper(): continue  # skip ERA5 reanalysis
                        try:
                            nd_,nlat_r,nlon_r,nt_r=load_nc(str(nf_path)); nd_=maybe_k_to_c(nd_,f"NatReload"); nd_=physical_qc(nd_,f"NatReload"); nd_=nd_[:T_common]
                            if nat_sum is None:
                                nat_shape=nd_.shape[1:]
                                nat_sum=np.where(np.isfinite(nd_),nd_,0.0).astype(np.float64)
                                nat_fin=np.isfinite(nd_).astype(np.int16); nat_lat_r=nlat_r; nat_lon_r=nlon_r
                            else:
                                if nd_.shape[1:]!=nat_shape:
                                    print(f"    NatReload skip {Path(nf_path).name}: shape {nd_.shape[1:]} != {nat_shape}")
                                    del nd_; gc.collect(); continue
                                fin=np.isfinite(nd_); nat_sum+=np.where(fin,nd_,0.0); nat_fin+=fin.astype(np.int16)
                            nat_cnt+=1
                        except: pass
                        finally:
                            try: del nd_; gc.collect()
                            except: pass
                    if nat_sum is None: raise ValueError("All reloads failed")
                    nat_fin[nat_fin==0]=1; mmm_nat_3d=(nat_sum/nat_fin.astype(np.float64)).astype(np.float32)
                    nat_dates=pd.date_range(f"{YR0}-01-01",periods=mmm_nat_3d.shape[0],freq="D")
                    obs_mmm_nat=extract_at_stations(mmm_nat_3d,nat_lat_r,nat_lon_r,obs_df,nat_dates)
                    obs_df["mmm_native_val"]=obs_mmm_nat["grid_val"].values
                    del mmm_nat_3d,nat_sum,nat_fin; gc.collect()
                except Exception as e:
                    obs_df["mmm_native_val"]=np.nan; print(f"  MMM native FAILED: {e}")
            else: obs_df["mmm_native_val"]=np.nan
        else: obs_df["mmm_native_val"]=np.nan
        # Fig B needs the full 3D MMM array to compute the time-mean bias and
        # RMSE against AORC, and it runs after the station metrics are built, so
        # keep a reference alive here and free it immediately after that call.
        mmm_for_figB=mmm

        print(f"{_elapsed()} [SI S1.4a — Monthly bias time series]")
        obs_monthly=obs_df.copy(); obs_monthly["ym"]=obs_monthly["date"].dt.to_period("M")
        _od=pd.to_datetime(obs_df["date"])
        print(f"  [obs-raw] GHCN-D daily records: {len(obs_df):,} rows, "
              f"{_od.min():%Y-%m-%d} -> {_od.max():%Y-%m-%d}")
        _yrc=obs_df.assign(_y=_od.dt.year).groupby("_y").size()
        print(f"  [obs-raw] years present: {int(_yrc.index.min())}-{int(_yrc.index.max())} "
              f"({_yrc.index.nunique()} distinct years)")
        _missing=[y for y in range(YR0,YR1+1) if y not in set(_yrc.index)]
        if _missing: print(f"  [obs-raw] *** YEARS WITH NO OBS in {YR0}-{YR1}: {_missing}")
        else: print(f"  [obs-raw] all years {YR0}-{YR1} have at least some obs")
        obs_gm=obs_monthly.groupby("ym")["tas_C"].agg(["mean","sem"])
        obs_dates=obs_gm.index.to_timestamp(); obs_vals=obs_gm["mean"].values; obs_se=obs_gm["sem"].values
        if len(obs_dates):
            print(f"    obs monthly coverage: {pd.Timestamp(obs_dates.min()):%Y-%m} -> "
                  f"{pd.Timestamp(obs_dates.max()):%Y-%m} ({len(obs_dates)} months)")
        mse_aorc=_monthly_mean_se(aorc,mask,gdates,se_map=aorc_se)
        mse_ds=_monthly_mean_se(era5ds,mask,gdates,se_map=era5ds_se)
        e5n_mask=np.ones(era5n.shape[1:],bool); mse_e5n=_monthly_mean_se(era5n,e5n_mask,e5dates)
        monthly_mse={"AORC":mse_aorc,"ERA5 downscaled":mse_ds,"ERA5 Native":mse_e5n,"GHCN-D obs":(obs_dates,obs_vals,obs_se)}
        if "MMM-Downscaled" in monthly:
            d_mmm,v_mmm=monthly["MMM-Downscaled"]; monthly_mse["MMM-Downscaled"]=(d_mmm,v_mmm,np.zeros_like(v_mmm))
        if mmm_native_monthly is not None:
            d_,v_=mmm_native_monthly; monthly_mse["MMM Native"]=(d_,v_,np.zeros_like(v_))

        sm={}; sok=[]
        for sid,grp in obs_df.groupby("id"):
            nm=grp["name"].iloc[0] if "name" in grp.columns else sid; d={"name":nm}
            for src,col in [("AORC","aorc_val"),("Downscaled","ds_val"),("ERA5nat","e5nat_val"),("MMM","mmm_val"),("MMM_native","mmm_native_val")]:
                v=np.isfinite(grp.tas_C)&np.isfinite(grp[col])
                if v.sum()<100: d[f"bias_{src}"]=np.nan; continue
                o,g=grp.loc[v,"tas_C"].values,grp.loc[v,col].values
                d[f"bias_{src}"]=np.mean(g-o); d[f"rmse_{src}"]=np.sqrt(np.mean((g-o)**2)); d[f"r_{src}"]=stats.pearsonr(o,g)[0]
            if any(np.isfinite(d.get(f"bias_{s}",np.nan)) for s in ["AORC","Downscaled","ERA5nat","MMM","MMM_native"]):
                sm[sid]=d; sok.append(sid)
        if sm:
            fig_stn_metrics(sm)
        if sok:
            ov=obs_df.loc[obs_df.id.isin(sok),"tas_C"].values
            av=obs_df.loc[obs_df.id.isin(sok),"aorc_val"].values
            dv=obs_df.loc[obs_df.id.isin(sok),"ds_val"].values
            e5v=obs_df.loc[obs_df.id.isin(sok),"e5nat_val"].values
            mv=obs_df.loc[obs_df.id.isin(sok),"mmm_val"].values if "mmm_val" in obs_df.columns else None
            mnv=obs_df.loc[obs_df.id.isin(sok),"mmm_native_val"].values if "mmm_native_val" in obs_df.columns else None
            print(f"{_elapsed()} [Fig 3 — Station validation summary]")
            fig_combined_stn_panel(ov,av,dv,e5v,monthly_mse,sm,obs_df,
                                   mmm_v=mv,mmmn_v=mnv)
            print(f"{_elapsed()} [Fig B — Spatial bias/RMSE maps + GHCN-D overlay]")
            fig_bias_rmse_maps(era5ds,mmm_for_figB,aorc,lat,mask,obs_df=obs_df,sm=sm)
            # No further figure needs the full 3D MMM array; release it.
            try: del mmm_for_figB, mmm; gc.collect()
            except NameError: pass
        with open(OUT_DIR/"Table_station_metrics.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["Number","ID","Name","Bias_AORC","RMSE_AORC","r_AORC","Bias_DS","RMSE_DS","r_DS","Bias_E5nat","RMSE_E5nat","r_E5nat"])
            for s in sm:
                d=sm[s]; num=stn_num.get(s,"?")
                w.writerow([num,s,d.get("name","")]+[f"{d.get(f'{mk}_{src}',np.nan):.2f}" for src in ["AORC","Downscaled","ERA5nat"] for mk in ["bias","rmse","r"]])
    else:
        print(f"{_elapsed()} [Station-based analysis SKIPPED \u2014 no GHCN-D obs for {LONGNAME}]")
        print(f"  Skipped: Fig 2 (station boxplots), Fig 3 (station validation),")
        print(f"           Fig C (bias attribution vs GHCN-D), SI S1.7 (station map),")
        print(f"           SI scatter helpers, SI ERA5 station validation,")
        print(f"           Table_station_metrics.csv.")
        print(f"  AORC is retained as the gridded reference for the spatial figures.")
        print(f"{_elapsed()} [Fig B — Spatial bias/RMSE maps (no station overlay)]")
        fig_bias_rmse_maps(era5ds,mmm,aorc,lat,mask)
        # mmm is large; free it as the obs block would have done at line `del mmm`.
        # Fig C is obs-dependent so we don't need to stash mmm_clim_for_figC.
        try: del mmm; gc.collect()
        except NameError: pass

    with open(OUT_DIR/"Table_gridded_summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["Model","Season","NormSD","PearsonR","CRMSE","ACC"])
        for row in csv1: w.writerow(row)

    print(f"\n{_elapsed()}  All -> {OUT_DIR.resolve()}")
    print(f"  Total runtime: {(_time.time()-_t0)/60:.1f} minutes")
    try:
        on_disk=sorted(p for p in OUT_DIR.iterdir() if p.suffix.lower() in (".pdf",".png"))
        print(f"\n  ====== Saved figures ({len(on_disk)} files in {OUT_DIR.resolve()}) ======")
        for p in on_disk:
            print(f"    {p.name:48s}  {p.stat().st_size/1024:>8.0f} KB")
        print(f"  =======================================================")
        if len(_SAVED_FIGURES)*2!=len(on_disk):
            print(f"  NOTE: Python recorded {len(_SAVED_FIGURES)} figure calls "
                  f"(expect {len(_SAVED_FIGURES)*2} files), but {len(on_disk)} found on disk.")
    except Exception as _e:
        print(f"  (could not list OUT_DIR contents: {_e})")

# ===================== RUNNER =====================
RUN_MODE="production"

_in_j=False
try:
    cfg=get_ipython().config
    if "IPKernelApp" in cfg or hasattr(get_ipython(),"kernel"): _in_j=True
except NameError: pass

def _run(mode="production",outdir=None):
    global OUT_DIR
    if outdir: OUT_DIR=Path(outdir); OUT_DIR.mkdir(parents=True,exist_ok=True)
    main()

if _in_j: _run(RUN_MODE)
else:
    import argparse; p=argparse.ArgumentParser()
    p.add_argument("--outdir",default=None)
    a=p.parse_args(); _run("production",a.outdir)
