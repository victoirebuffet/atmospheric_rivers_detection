#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:53:13 2025

@author: vickybuffy
"""

import json
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ar_detection_funcs import *

def load_config(config_path):
    """Load configuration from a JSON file."""
    with open(config_path, 'r') as file:
        config = json.load(file)
    
    # Handle years: explicitly listed or range
    if 'years' not in config and 'start_year' in config and 'end_year' in config:
        start_year = config.pop('start_year')
        end_year = config.pop('end_year')
        config['years'] = list(range(start_year, end_year + 1))
    elif 'years' in config:
        # Ensure 'years' is a list
        if not isinstance(config['years'], list):
            raise ValueError("The 'years' field must be a list if provided explicitly.")
    else:
        raise ValueError("The configuration must specify either 'years' or both 'start_year' and 'end_year'.")
    
    return config

def main(config_path):
    # Load the configuration
    config = load_config(config_path)
    
    # Extract configuration parameters
    scheme = config.get('scheme')
    years = config.get('years')
    scan_extent = config.get('scan_extent')
    percentile = config.get('percentile')
    source_path = config.get('source_path', '.')
    dir_path = config.get('dir_path', '.')
    hemisphere = config.get('hemisphere', 'ant')
    adaptative_threshold = config.get('adaptative_threshold', False)
    
    if adaptative_threshold:
        print("Computing the adaptive threshold.")
        compute_adaptative_threshold(years, scan_extent, source_path=source_path)
    
    print("Computing the percentile threshold.")
    get_percentile(scheme=scheme, years=years, scan_extent=scan_extent, percentile=percentile, hemisphere=hemisphere, source_path=source_path)
    
    print("Extracting ARs.")
    # Process each year in the configuration
    for year in years:
        print(f"Processing year: {year}")
        extract_ARs(
            scheme=scheme,
            year=year,
            percentile=percentile,
            scan_extent=scan_extent,
            hemisphere=hemisphere,
            dir_path=dir_path,
            source_path=source_path,
            adaptative_threshold=adaptative_threshold
        )

if __name__ == "__main__":
    # Specify the path to the configuration file
    config_file_path = "config.json"
    main(config_file_path)
