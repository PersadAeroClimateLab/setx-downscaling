import xarray as xr
from pathlib import Path
from os.path import isfile
from os import listdir
from uuid import uuid4
import zarr



if __name__ == "__main__":
    head_input_dir = "/local1/SETx_UIFL_Zarr"

    variable_units = {
        "lat": "degrees_north",
        "lon": "degrees_east",
        "tas": "K",
        "tasmax": "K",
        "tasmin": "K",
        "hurs": "%",
        "huss": "1",
        "pr": "mm day-1",
        "rlds": "W m-2",
        "rsds": "W m-2",
        "sfcWind": "m s-1",
    }

    variable_names = {
        "lat": "latitude",
        "lon": "longitude",
        "time": "time",
        "tas": "air_temperature",
        "tasmax": "air_temperature",
        "tasmin": "air_temperature",
        "hurs": "relative_humidity",
        "huss": "specific_humidity",
        "pr": "lwe_precipitation_rate", # could be precipitation_amount depending on units
        "rlds": "surface_downwelling_longwave_radiation",
        "rsds": "surface_downwelling_shortwave_radiation",
        "sfcWind": "wind_speed",
    }

    variable_cellmethods = {
        "tas": "time: mean",
        "tasmax": "time: max",
        "tasmin": "time: min",
        "hurs": "time: mean",
        "huss": "time: mean",
        "pr": "time: mean",
        "rlds": "time: mean",
        "rsds": "time: mean",
        "sfcWind": "time: mean",
    }

    global_variable_attrs = {
        "cell_measures": "area: areacella",
        "grid_mapping": "crs"
    }

    for var_dir in listdir(head_input_dir):
        for zarr_name in listdir(f"{head_input_dir}/{var_dir}"):
            if ".zarr" not in zarr_name:
                continue
            
            ds_path = f"{head_input_dir}/{var_dir}/{zarr_name}"
            print(ds_path)
            ds = xr.open_zarr(ds_path, decode_times=False)
            
            model, scenario, variable, period = ds_path.split("/")[-1].split("_")

            if "variant_label" in ds.attrs:
                ds.attrs.pop("variant_label")
            if "activity" in ds.attrs:
                ds.attrs.pop("activity")
            if "institution" in ds.attrs:
                ds.attrs.pop("institution")
            if "source" in ds.attrs:
                ds.attrs.pop("source")
            if "scenario" in ds.attrs:
                ds.attrs.pop("scenario")
            if "experiment" in ds.attrs:
                ds.attrs.pop("experiment")
            if "resolution_id" in ds.attrs:
                ds.attrs.pop("resolution_id")
            if "doi" in ds.attrs:
                ds.attrs.pop("doi")

            if scenario == "hist":
                scenario = "historical"
                
            ds.attrs["activity_id"] = "SETx-1km-CMIP6"
            ds.attrs["region_id"] = "NAM"
            ds.attrs["institution_id"] = "SETx-UTAustin"
            ds.attrs["source_id"] = "KrigR"
            ds.attrs["driving_mip_era"] = "CMIP6"
            ds.attrs["driving_activity_id"] = "CMIP"
            ds.attrs["driving_source_id"] = model
            ds.attrs["driving_experiment_id"] = scenario
            ds.attrs["driving_variant_label"] = "r1i1p1f1"
            ds.attrs["variable_id"] = variable
            ds.attrs["version"] = "v20260720"
            ds.attrs["Conventions"] = "CF-1.8"
            ds.attrs["nominal_resolution"] = "1 km"
            ds.attrs["geospatial_lat_resolution"] = "0.008333333 degrees"
            ds.attrs["geospatial_lon_resolution"] = "0.008333333 degrees"
            ds.attrs = dict(sorted(ds.attrs.items())) 
            
            for variable_name in ds.variables:
                new_attrs = {}
                standard_err = False
                if "SE" in variable_name:
                    variable_name = variable_name.split("_")[0]
                    standard_err = True
                if variable_name in variable_units:
                    new_attrs |= {
                        "units": variable_units[variable_name],
                    }
                if variable_name in variable_names:
                    new_attrs |= {
                        "standard_name": variable_names[variable_name],
                    }
                if variable_name in variable_cellmethods:
                    new_attrs |= {
                        "cell_methods": variable_cellmethods[variable_name],
                    } | global_variable_attrs
            
                if standard_err:
                    new_attrs["standard_name"] += " standard_error"
                    new_attrs["comment"] = "Standard error introduced by the statistical downscaling method."
                    ds[variable_name + "_SE"].attrs |= new_attrs
                else:
                    ds[variable_name].attrs |= new_attrs
            
            # Open the existing store in read/write mode (metadata only, no data copy)
            group = zarr.open_group(ds_path, mode="r+")
            
            # Write global attrs
            group.attrs.put(ds.attrs)   # or .update(ds.attrs) if you want to merge rather than replace
            
            # Write per-variable attrs
            for name, var in ds.variables.items():
                if name in group:
                    group[name].attrs.put(dict(var.attrs))
            
            # If the store uses consolidated metadata (.zmetadata), refresh it —
            # otherwise readers that use consolidated=True will see stale attrs
            zarr.consolidate_metadata(ds_path)
    
    output_head_dir = "/home/s1/persad/SETx_UIFL"
    head_input_dir = "/local1/SETx_UIFL_Zarr"
    
    for var_dir in listdir(head_input_dir):
        if var_dir == "pr":
            continue
        for zarr_name in listdir(f"{head_input_dir}/{var_dir}"):
            if ".zarr" not in zarr_name:
                continue
            
            ds_path = f"{head_input_dir}/{var_dir}/{zarr_name}"
    
            ds = xr.open_zarr(ds_path)
            a = ds.attrs
            variable_id = a["variable_id"]
            
            # Build the DRS directory path from global attributes
            out_dir = Path(output_head_dir,
                a["activity_id"], a["region_id"], a["institution_id"], a["source_id"],
                a["driving_mip_era"], a["driving_activity_id"], a["driving_experiment_id"],
                a["driving_source_id"], a["driving_variant_label"], a["frequency"],
                variable_id, a["version"],
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            
            # Common filename prefix (everything except the time range)
            fname_prefix = "_".join([
                variable_id, a["region_id"], a["institution_id"], a["source_id"],
                a["driving_mip_era"], a["driving_experiment_id"], a["driving_source_id"],
                a["driving_variant_label"], a["frequency"],
            ])
            # Split into yearly chunks and write one file per year
            for year, ds_year in ds.groupby("time.year"):
                start = ds_year.time.min().dt.strftime("%Y%m%d").item()
                end = ds_year.time.max().dt.strftime("%Y%m%d").item()
                ds_year.attrs["tracking_id"] = f"hdl:21.14100/{uuid4()}"
                
                out_path = out_dir / f"{fname_prefix}_{start}-{end}.nc"
                if not isfile(out_path):
                    ds_year.to_netcdf(out_path, engine="h5netcdf")
                else:
                    print(f"Skipping {out_path}")
            print(ds_path)
