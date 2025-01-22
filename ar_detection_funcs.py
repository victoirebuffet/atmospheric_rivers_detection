#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 19:55:06 2025

@author: vickybuffy
"""

import xarray as xr
import os
import numpy as np
from dask.diagnostics import ProgressBar
from scipy.ndimage import gaussian_filter1d
import glob


def compute_adaptative_threshold(years, scan_extent, source_path='.'):
    
    if os.path.exists(os.path.join(source_path, "adaptative_threshold_IWV.nc")):
        print("Adaptative threshold already computed. Skipping this step.")
    else:
        scan_extent_sorted = validate_latitudes(scan_extent)
        
        list_yearly_mean = []
        
        # Collect yearly means
        for year in years:
            filepath_tcw = os.path.join(source_path, f"total_column_water_{year}_reanaHS.nc")
            tcw = xr.open_dataset(filepath_tcw)['tcw'].sel(latitude=slice(-scan_extent_sorted[0], -scan_extent_sorted[1]))
            yearly_mean = tcw.mean().item()  # Convert to scalar
            list_yearly_mean.append(yearly_mean)
        
        # Compute overall mean and smoothed yearly means
        overall_mean = np.mean(list_yearly_mean)
        smoothed_yearly_mean = gaussian_filter1d(list_yearly_mean, sigma=2)
        adjusted_threshold = overall_mean / smoothed_yearly_mean  # Ratio for threshold adjustment
        
        # Store results in an xarray.Dataset
        ds = xr.Dataset(
            {
                "yearly_mean": (["year"], list_yearly_mean),
                "smoothed_yearly_mean": (["year"], smoothed_yearly_mean),
                "adaptative_threshold": (["year"], adjusted_threshold)
            },
            coords={"year": years}
        )
        save_to_netcdf(ds, os.path.join(source_path, "adaptative_threshold_IWV.nc"), progress_bar=True)
        
def validate_latitudes(scan_extent):
    """Validate latitude extents."""
    scan_extent_sorted = sorted([abs(i) for i in scan_extent])
    if max(scan_extent_sorted) > 90 or min(scan_extent_sorted) < 0:
        raise ValueError("Latitude extents must be within [-90, 90].")
    return scan_extent_sorted


def get_file_path(source_path, long_varname, year):
    """Construct the file path for NetCDF files."""
    pattern = os.path.join(source_path, f"{long_varname}_{year}_reanaHS.nc")
    matching_files = glob.glob(pattern)
    
    if not matching_files:
        raise FileNotFoundError(f"No files match the pattern: {pattern}")
    
    return matching_files[0]

def load_and_process_dataset(filepath, varname, month, hemisphere, scan_extent_sorted):
    """Load and preprocess dataset for a specific month and hemisphere."""
    ds = xr.open_dataset(filepath)[varname]
    ds = ds.sel(time=ds['time.month'] == month)
    if hemisphere == 'ant':
        ds = ds.sel(latitude=slice(-scan_extent_sorted[0], -scan_extent_sorted[1]))
    elif hemisphere == 'arc':
        ds = ds.sel(latitude=slice(scan_extent_sorted[1], scan_extent_sorted[0]))
    else:
        raise ValueError(f"Unknown hemisphere: {hemisphere}")
    return ds


def save_to_netcdf(ds, filename, progress_bar=True):
    """Save xarray Dataset to NetCDF with optional progress bar."""
    delayed_obj = ds.to_netcdf(filename, compute=False, mode='w', format='NETCDF4')
    if progress_bar:
        with ProgressBar():
            delayed_obj.compute()
    else:
        delayed_obj.compute()


def get_percentile(scheme, years, scan_extent, percentile, hemisphere='ant', source_path='.'):
    scan_extent_sorted = validate_latitudes(scan_extent)

    scheme_map = {
        'vIVT': ('vertical_integral_of_northward_water_vapour_flux', 'p72.162'),
        'IWV curved': ('total_column_water', 'tcw')
    }

    if scheme not in scheme_map:
        raise ValueError(f"Unknown scheme: {scheme}")

    long_varname, varname = scheme_map[scheme]

    for month in range(1, 13):
        if os.path.exists(os.path.join(source_path, f"{scheme}_per{percentile}_{month}_{hemisphere}.nc")):
            print(f"Percentile for month {month} already computed. Skipping this one.")
            continue
    
        combined_ds = None
        for year in years:
            try:
                filepath = get_file_path(source_path, long_varname, year)
                ds = load_and_process_dataset(filepath, varname, month, hemisphere, scan_extent_sorted)

                if combined_ds is None:
                    combined_ds = ds
                else:
                    combined_ds = xr.concat([combined_ds, ds], dim='time')
            except FileNotFoundError as e:
                print(e)
                continue

        if combined_ds is None:
            print(f"No data found for scheme {scheme} in month {month}. Skipping...")
            continue

        if scheme == 'vIVT' and hemisphere == 'ant':
            combined_ds = -combined_ds

        ds_per = combined_ds.quantile(percentile / 100, dim='time', skipna=True).rename('per').to_dataset()
        ds_per.attrs['description'] = f"{percentile}th Percentile of {scheme} in ERA-5"

        output_filename = os.path.join(source_path, f"{scheme}_per{percentile}_{month}_{hemisphere}.nc")
        save_to_netcdf(ds_per, output_filename)


def filter_ARs_varying_length(ds, base_len=1000, max_len=2000, lat_min=-85, lat_max=-70):
    """
    Filter atmospheric rivers (ARs) by their varying lengths based on latitude.

    Parameters
    ----------
    ds : xarray.DataArray
        Binary dataset of AR indicators (1s).
    base_len : int
        Minimum length of the AR at `lat_min` in km.
    max_len : int
        Maximum length of the AR at `lat_max` in km.
    lat_min : float
        Latitude corresponding to `base_len`.
    lat_max : float
        Latitude corresponding to `max_len`.

    Returns
    -------
    xarray.DataArray
        Filtered AR binary dataset with the same shape as the input.
    """
    
    from scipy.ndimage import label
    from numpy.polynomial.polynomial import Polynomial
    from geopy.distance import geodesic
    
    def dynamic_length_requirement(lat, base_len, max_len, lat_min, lat_max):
        """Calculate the dynamic length requirement based on latitude."""
        if lat >= lat_max:
            return max_len
        elif lat <= lat_min:
            return base_len
        else:
            return base_len + (max_len - base_len) * (lat - lat_min) / (lat_max - lat_min)

    def polynomial_length_in_km(poly, lat_min, lat_max, num_points=100):
        """Calculate the length of the polynomial curve in kilometers using geodesic distances."""
        lats = np.linspace(lat_min, lat_max, num_points)
        lons = poly(lats)
        distances = [
            geodesic((lats[i], lons[i]), (lats[i + 1], lons[i + 1])).kilometers
            for i in range(len(lats) - 1)
        ]
        return np.sum(distances)

    def calculate_ar_length(poly, lat_min, lat_max, num_points=100):
        """Calculate the length of an AR using geodesic distances."""
        lats = np.linspace(lat_min, lat_max, num_points)
        lons = poly(lats)
        distances = [
            geodesic((lats[i], lons[i]), (lats[i + 1], lons[i + 1])).kilometers
            for i in range(len(lats) - 1)
        ]
        return np.sum(distances)

    structure = np.ones((3, 3), dtype=int)
    out = np.zeros_like(ds, dtype=int)

    for t in range(ds.sizes['time']):
        
        binary_slice = ds.isel(time=t).values
        labeled, num_features = label(binary_slice, structure=structure)

        for label_id in range(1, num_features + 1):
            coords = np.array(np.where(labeled == label_id)).T
            latitudes = ds.latitude.values[coords[:, 0]]

            lat_min_AR, lat_max_AR = latitudes.min(), latitudes.max()
            required_length = dynamic_length_requirement(lat_min_AR, base_len, max_len, lat_min, lat_max)

            # Fit polynomial to determine AR length
            longitudes = ds.longitude.values[coords[:, 1]]
            if (longitudes.max() - longitudes.min()) > 180:
                longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)

            best_poly = None
            best_r_squared = -np.inf
        
            for degree in range(3, 0, -1):  # Test polynomial degrees from 3 down to 1
                try:
                    poly = Polynomial.fit(latitudes, longitudes, deg=degree)
                    predicted_lons = poly(latitudes)
        
                    # Compute R-squared
                    residual_sum_of_squares = np.sum((longitudes - predicted_lons) ** 2)
                    total_sum_of_squares = np.sum((longitudes - longitudes.mean()) ** 2)
                    r_squared = 1 - (residual_sum_of_squares / total_sum_of_squares)
        
                    # Check if this degree has the best goodness of fit
                    if r_squared > best_r_squared:
                        best_r_squared = r_squared
                        best_poly = poly
                except np.linalg.LinAlgError:
                    continue
        
            if best_poly is None:
                continue  # Skip if no valid polynomial was found
        
        ar_length = polynomial_length_in_km(best_poly, lat_min_AR, lat_max_AR)
    
        if ar_length >= required_length:
            out[t, coords[:, 0], coords[:, 1]] = 1

    unique_longitudes = np.unique(ds['longitude'].values)
    outy = out[:, :, :len(unique_longitudes)] + out[:, :, len(unique_longitudes):]
    outy[outy > 1] = 1
    
    ds_short = ds.groupby('longitude').max(dim='longitude')
    out_xr = xr.DataArray(
        data=outy,
        dims=ds_short.dims,
        coords=ds_short.coords,
        attrs=ds_short.attrs,
    )

    return out_xr


def filter_ARs_meridional_length(ds, min_lat_extent=20):
    """
    Filter ARs based on their meridional extent (latitude coverage).

    Parameters
    ----------
    ds : xarray.DataArray
        Binary dataset of AR indicators (1s).
    min_lat_extent : int
        Minimum meridional extent of ARs in degrees latitude.

    Returns
    -------
    xarray.DataArray
        Filtered AR binary dataset with the same shape as the input.
    """
    from skimage.measure import label, regionprops
    

    lat_resolution = abs(ds.latitude.diff(dim='latitude').min().values)
    required_pixels = int(min_lat_extent / lat_resolution)

    out = np.zeros_like(ds, dtype=int)

    for t in range(ds.sizes['time']):
        labeled, num_features = label(ds.isel(time=t).values, connectivity=2)

        for region in regionprops(labeled):
            if (region.bbox[2] - region.bbox[0]) >= required_pixels:
                coords = region.coords
                out[t, coords[:, 0], coords[:, 1]] = 1

    # Handle longitude wrapping
    unique_longitudes = np.unique(ds['longitude'].values)
    outy = out[:, :, :len(unique_longitudes)] + out[:, :, len(unique_longitudes):]
    outy[outy > 1] = 1
     
    ds_short = ds.groupby('longitude').max(dim='longitude')
    out_xr = xr.DataArray(
        data=outy,
        dims=ds_short.dims,
        coords=ds_short.coords,
        attrs=ds.attrs,
    )
    return out_xr


def extract_ARs(scheme, year, percentile, scan_extent, hemisphere='ant', dir_path='.', source_path='.', adaptative_threshold=False):
    """
    Extract and filter atmospheric rivers (ARs) based on a scheme and percentile threshold.

    Parameters
    ----------
    scheme : str
        AR scheme ('vIVT' or 'IWV curved').
    year : int
        Year of the data to process.
    percentile : float
        Percentile threshold for AR detection.
    scan_extent : list
        Latitude extent to scan [min, max].
    hemisphere : str, optional
        Hemisphere ('ant' for Antarctic, 'arc' for Arctic).
    dir_path : str, optional
        Output directory path.
    source_path : str, optional
        Directory containing the input data.

    adaptative_threshold : bool, optional
        Whether to use an adaptative humidity threshold to take into account Clausius Clapeyron with global warming
    Returns
    -------
    None
    """
    scan_extent_sorted = validate_latitudes(scan_extent)
    scheme_map = {
    'vIVT': ('vertical_integral_of_northward_water_vapour_flux', 'p72.162', filter_ARs_meridional_length),
    'IWV curved': ('total_column_water', 'tcw', filter_ARs_varying_length)
}
    if scheme not in scheme_map:
        raise ValueError(f"Unknown scheme: {scheme}")

    long_varname, varname, filter_function = scheme_map[scheme]
    combined_ds = None

    for month in range(1, 13):
        try:
            ds_open = xr.open_dataset(
                os.path.join(source_path, f"{long_varname}_{year}_reanaHS.nc")
            )[varname]
            ds_open = ds_open.sel(time=ds_open['time.month'] == month)

            # Load percentile thresholds
            ds_per_path = os.path.join(source_path, f"{scheme}_per{percentile}_{month}_{hemisphere}.nc")
            if not os.path.exists(ds_per_path):
                raise FileNotFoundError(f"Percentile file not found: {ds_per_path}")
            ds_per = xr.open_dataset(ds_per_path)['per']

            # Subset data by hemisphere
            if hemisphere == 'ant':
                ds_open = ds_open.sel(latitude=slice(-scan_extent_sorted[0], -scan_extent_sorted[1]))
                ds_per = ds_per.sel(latitude=slice(-scan_extent_sorted[0], -scan_extent_sorted[1]))
            elif hemisphere == 'arc':
                ds_open = ds_open.sel(latitude=slice(scan_extent_sorted[1], scan_extent_sorted[0]))
                ds_per = ds_per.sel(latitude=slice(scan_extent_sorted[1], scan_extent_sorted[0]))

            # Apply scheme-specific adjustments
            if scheme == 'vIVT' and hemisphere == 'ant':
                ds_open = xr.where(ds_open < 0, -ds_open, np.nan)
            elif scheme == 'vIVT' and hemisphere == 'arc':
                ds_open = xr.where(ds_open > 0, ds_open, np.nan)
            
            if adaptative_threshold:
                adaptative_threshold_ds = xr.open_dataset(f'{source_path}"adaptative_threshold_IWV.nc')
                factor = adaptative_threshold_ds['smoothed_yearly_mean'].sel(year=year)/adaptative_threshold_ds['smoothed_yearly_mean'].mean(dim='year')
                threshold = ds_per * factor
                del(adaptative_threshold_ds, factor)
            else: 
                threshold = ds_per
                
            # Calculate binary result
            excess = ds_open - threshold
            binary_result = xr.where(excess > 0, 1, 0)
            concat_bin_res = xr.concat([binary_result, binary_result], dim='longitude')

            # Apply the filter function from the scheme_map
            filtered_ars = filter_function(concat_bin_res)

            # Combine results for all months
            filtered_ds = xr.Dataset(
                {"ar_binary_tag": (["time", "latitude", "longitude"], filtered_ars.data)},
                coords={"time": filtered_ars.time, "latitude": filtered_ars.latitude, "longitude": filtered_ars.longitude}
            )
            combined_ds = xr.concat([combined_ds, filtered_ds], dim="time") if combined_ds else filtered_ds

        except Exception as e:
            print(f"Error processing month {month}: {e}")
    
    if adaptative_threshold:
        adapt='_adaptative'
    else:
        adapt=''
    # Save combined dataset
    output_filename = os.path.join(
        dir_path,
        f"ar_tag_{year}_era5_{scheme.replace(' ', '_')}_{percentile}p_{scan_extent[0]}-{scan_extent[-1]}_{hemisphere}{adapt}.nc4"
    )
    save_to_netcdf(combined_ds, output_filename)
