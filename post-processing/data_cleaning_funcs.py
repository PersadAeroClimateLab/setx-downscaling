import xarray as xr
import numpy as np
import dask
import dask.array as da
from os.path import isdir, getsize
from os import listdir
import cftime

def write_ds(zarr_path):
    """
    Writes structure for template zarr store
    """
    input_ds = xr.open_zarr(zarr_path, consolidated=False, chunks=None).rename({"latitude": "lat", "longitude": "lon"}).drop_vars("crs")
    
    model, scenario, variable, start, end = zarr_path.split("/")[-1].split("_")
    if scenario == "hist":
        template_zarr_path = f"{output_dir}/{model}_{scenario}_{variable}_1950-2014.zarr"
    else:
        template_zarr_path = f"{output_dir}/{model}_{scenario}_{variable}_2015-2100.zarr"
    template_zarr = xr.open_zarr(template_zarr_path)

    input_ds = input_ds.assign_coords(dict(lat=template_zarr.lat, lon=template_zarr.lon))
    if template_zarr.sizes != input_ds.sizes:
        print(f"{zarr_path} <-> {template_zarr_path}")
    else:
        for t_start in range(0, len(input_ds.time), 300):
            chunk = input_ds.isel(time=slice(t_start, t_start + 300)).load()
            chunk.to_zarr(template_zarr_path, region="auto")


def build_template_zarrs(variable, output_dir):
    """
    Builds empty template zarr store for variable
    """
    scenarios = ['ssp370', 'ssp245', 'ssp126', 'ssp585', 'hist']
    
    models_calendars = {
        'BCC-CSM2-MR': "365_day",
        'CMCC-ESM2': "365_day",
        'CNRM-ESM2-1': "standard",
        'EC-Earth3': "proleptic_gregorian",
        'FGOALS-g3': "365_day",
        'GFDL-CM4': "365_day",
        'MPI-ESM1-2-HR': "proleptic_gregorian",
        'MRI-ESM2-0': "proleptic_gregorian",
        'NorESM2-MM': "365_day",
        'CESM2': "365_day",
    }
    
    for model in models_calendars:
        calendar = models_calendars[model]
        end_day = 31
        if calendar == "360_day":
            end_day = 30
        
        for scenario in scenarios:
            if scenario == "hist":
                tag = "1950-2014"
                times = xr.date_range(start="1950-01-01", end=f"2014-12-{end_day}", calendar=calendar, use_cftime=True)
            else:
                tag = "2015-2100"
                times = xr.date_range(start="2015-01-01", end=f"2100-12-{end_day}", calendar=calendar, use_cftime=True)
                
            lat_start = 33.745693949461256
            lat_end = 28.504027303764904
            lats = np.linspace(lat_start, lat_end, 630)
            
            lon_start = -97.49597255223884
            lon_end = -93.50430590154502
            lons = np.linspace(lon_start, lon_end, 480)
            
            global_attrs = {
                'activity': 'NEX-GDDP-CMIP6',
                'contact': 'Dr. Geeta Persad: geeta.persad@jsg.utexas.edu, Dr. Ifeanyichukwu Nduka: icnduka@jsg.utexas.edu',
                'Conventions': 'CF-1.7',
                'frequency': 'day',
                'product': 'output',
                'realm': 'atmos',
                'scenario': f"{scenario}",
                'version': '1.0',
                'disclaimer': 'This data is considered provisional and subject to change. This data is provided as is without any warranty of any kind, either express or implied, arising by law or otherwise, including but not limited to warranties of completeness, non-infringement, accuracy, merchantability, or fitness for a particular purpose. The user assumes all risk associated with the use of, or inability to use, this data.',
                'cmip6_source_id': f"{model}",
                'variant_label': 'r1i1p1f1',
                'resolution_id': '0.25 degree',
                'doi': 'https://doi.org/10.7917/OFSG3345',
                'external_variables': 'areacella',
                'institution': 'Jackson School of Geosciences, UT Austin, 23 San Jacinto Blvd, Austin, TX 78712',
                'source': 'BCSD; KrigR',
                'references': 'BCSD method: Thrasher et al., 2012, Hydrol. Earth Syst. Sci.,16, 3309-3314. KrigR method: Erik Kusch and Richard Davy 2022 Environ. Res. Lett. 17 024005. Ref period obs: latest version of the Princeton Global Meteorological Forcings (http://hydrology.princeton.edu/data.php), based on Sheffield et al., 2006, J. Climate, 19 (13), 3088-3111.',
                'title': f'{model}, {scenario}, global downscaled CMIP6 climate projection data'
            }
            
            ds = xr.Dataset(
                {
                    f"{variable}": (["time", "lat", "lon"], da.empty((len(times), len(lats), len(lons)), dtype=np.float32, chunks=(200, 630, 480))),
                    f"{variable}_SE": (["time", "lat", "lon"], da.empty((len(times), len(lats), len(lons)), dtype=np.float32, chunks=(200, 630, 480))),
                    "crs": None
                },
                coords={"time": times, "lat": lats, "lon": lons},
                attrs=global_attrs,
            )
            
            ds["lat"].attrs = {
                'units': 'degrees_north',
                'long_name': 'latitude',
                'standard_name': 'latitude',
                'axis': 'Y'
            }
            
            ds["lon"].attrs = {
                'units': 'degrees_east',
                'long_name': 'longitude',
                'standard_name': 'longitude',
                'axis': 'X'
            }
            
            ds["crs"].attrs = {
                'proj4': '+proj=longlat +datum=WGS84 +no_defs',
                'geotransform': '-97.50013921888883317024 0.008333333299999997126339 0 33.7498606161112562063 0 -0.008333333300000009269404'
            }
            
            ds.to_zarr(f"{output_dir}{model}_{scenario}_{variable}_{tag}.zarr", compute=False, zarr_format=2)


def scan_zarr(path, variable):
    """
    Scans zarr store variable for nans
    """
    nan_count = xr.open_zarr(path)[variable].isnull().sum(dim=["lat", "lon"]).compute()
    counts = np.bincount(nan_count.values.flatten())
    most_common_int = np.argmax(counts)
    return (path, variable, (nan_count.where(nan_count != most_common_int, drop=True) - most_common_int))


def get_slice_orders(head_input_dir, zarr_dir):
    """
    Checks input directory for timestep slice netcdf files
    """
    swap_orders = {}
    
    for filename in listdir(head_input_dir):
        if isdir(f"{head_input_dir}/{filename}") or "_9km_" in filename or "UKESM" in filename:
            continue
        se = None
        if "Kriging_" not in filename and ".nc" in filename:
            if "SE_" in filename:
                se, model, scenario, variable, km, date = filename.split("_")
            else:
                model, scenario, variable, km, date = filename.split("_")
        else:
            continue
        date = date.split(".")[0]

        if scenario == "ssp1":
            scenario = "ssp126"
        elif scenario == "ssp2":
            scenario = "ssp245"
        elif scenario == "ssp3":
            scenario = "ssp370"
        elif scenario == "ssp5":
            scenario = "ssp585"

        if len(date) == 8 and "-" not in date:
            date = f"{date[:4]}-{date[4:6]}-{date[6:]}"
        
        if getsize(f"{head_input_dir}/{filename}") < 10*1024:
            continue
        
        if scenario == "hist":
            zarr_path = f"{zarr_dir}/{model}_{scenario}_{variable}_1950-2014.zarr"
        else:
            zarr_path = f"{zarr_dir}/{model}_{scenario}_{variable}_2015-2100.zarr"
        assert isdir(zarr_path)
    
        if se is not None:
            variable = f"{variable}_SE"
        order_args = (f"{head_input_dir}/{filename}", date, variable)
        
        if zarr_path in swap_orders:
            swap_orders[zarr_path].append(order_args)
        else:
            swap_orders[zarr_path] = [order_args]
    return swap_orders


def swap_slice(template_zarr_path, swap_args):
    """
    Takes output from get_slice_orders and writes slices to zarr store
    """
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    template_zarr = xr.open_zarr(template_zarr_path, decode_times=time_coder)
    calendar = template_zarr.time.values[0].calendar
    
    for swap_path, date, variable in swap_args:
        swap_ds = xr.open_dataset(swap_path).rename({"latitude":"lat", "longitude":"lon"})
        if swap_ds.lat.size != 630 or swap_ds.lon.size != 480:
            continue
        for swap_var in swap_ds.variables:
            if swap_var != "crs" and len(swap_ds[swap_var].dims) == 2:
                swap_ds = swap_ds[swap_var]
                break
        date = cftime.datetime.strptime(date, "%Y-%m-%d", calendar=calendar)
        swap_ds = swap_ds.rename(variable)
        swap_ds = swap_ds.expand_dims("time").assign_coords(dict(time=[date], lat=template_zarr.lat, lon=template_zarr.lon))
        swap_ds.to_zarr(template_zarr_path, region="auto")


def add_full_ds(nc_path, zarr_path):
    """
    Reads post-proc netCDF file and adds it to the target zarr_path
    """
    zarr_ds = xr.open_zarr(zarr_path)
    swap_ds = xr.open_dataset(nc_path).rename({"latitude":"lat", "longitude":"lon"})
    swap_ds = swap_ds.assign_coords(lat=zarr_ds.lat, lon=zarr_ds.lon).drop_vars("crs")
    swap_ds.to_zarr(zarr_path, region="auto")
    

def fill_initial_dataset(zarr_dir):
    """
    Uses zarr template names to find existing netCDF files to fill the store with
    """
    nc_products_dir1 = "/projects/dgs/persad_research/SETx_UIFL_Products"
    nc_products_dir2 = "/projects/dgs/persad_research/SETx_UIFL_Staging"
    
    orders = []
    
    for zarr_name in listdir(zarr_dir):
        model, scenario, variable, yrs = zarr_name.split("_")
        nc_dir = f"{nc_products_dir1}/{model}/{scenario}/day/{variable}"
        if not isdir(nc_dir):
            nc_dir = f"{nc_products_dir2}/{model}/{scenario}/day/{variable}"
        if not isdir(nc_dir):
            continue
    
        for nc_name in listdir(nc_dir):
            orders.append(dask.delayed(add_full_ds)(f"{nc_dir}/{nc_name}", f"{zarr_dir}/{zarr_name}"))
    return orders

def export_missing_slices(file_path, results):
    """
    Exports results from scan_zarr to text file
    """
    with open(file_path, 'w') as f: 
        for path, variable, array in results:
            if array.size != 0:
                model, scenario, var, ext = path.split("/")[-1].split("_")
                for ts in array.time.values:
                    ts = str(ts)
                    if "T" in ts:
                        ts = ts.split("T")[0]
                    else:
                        ts = ts.split(" ")[0]
                    f.write(f"{model}, {variable}, {scenario}, {ts}\n")