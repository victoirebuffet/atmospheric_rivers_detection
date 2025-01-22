# Atmospheric River Detection Tool

## Overview
This tool is designed to detect Atmospheric Rivers (ARs) from the ERA5 dataset (Hersbach et al., 2020) using two distinct schemes:

1. **vIVT AR Detection**:
   - Follows the methodology outlined in Wille et al. (2019, 2021).
   - Utilizes vertical integral of northward water vapor transport (vIVT) as the input variable.
   - Retains regions where poleward vIVT exceeds the monthly climatological Xth percentile.
   - Labels as ARs only those continuous transports with a meridional extent greater than 20°.

2. **IWV Curved Detection**:
   - A novel approach based on total column water (IWV) or integrated water vapor, with an absolute length requirement that varies by latitude.
   - Retains regions where IWV exceeds the monthly climatological Xth percentile.
   - Labels as ARs only those continuous transports that exceed a required length. The required length is 2000 km for transports with their most poleward point below 70° and decreases linearly to 1000 km for objects reaching 85°.
   - Represents an upcoming AR Detection Tool designed to characterize intense and long water vapor transports not restricted to Antarctica.

The threshold characterized by the Xth percentile can be modified, allowing the tool to be applied to either hemisphere (Antarctic or Arctic), within user-specified latitude ranges, and over a set of years.

Additionally, the tool includes an **adaptive threshold** feature, which accounts for the increased humidity in the atmosphere due to global warming, following the methodology described in Barthélemy et al. (to be published).

---

## Features
- **Detection Schemes**: Choose between the vIVT or IWV Curved scheme for AR detection.
- **Hemisphere Support**: Apply the detection to either the Antarctic or Arctic regions.
- **Custom Latitude Ranges**: Specify latitude bounds to focus on specific areas.
- **Yearly Analysis**: Process data over multiple years, as specified by the user.
- **Adaptive Threshold**: Incorporates an adaptative humidity threshold to adjust for the Clausius-Clapeyron relationship under global warming conditions.

---

## Requirements
- **Programming Language**: Python 3.6+
- **Dependencies**:
  - `xarray`
  - `numpy`
  - `scipy`
  - `dask`
  - `glob`
  - `geopy`
  - `skimage`
  - `matplotlib` (optional for visualization)

---

## Usage

### Input Configuration
The tool uses a JSON configuration file to specify parameters. Example:

```json
{
    "scheme": "vIVT",
    "start_year": 1979,
    "end_year": 2023,
    "scan_extent": [-85, -15],
    "percentile": 98,
    "source_path": "/path/to/source/files",
    "dir_path": "/path/to/output/files",
    "hemisphere": "ant",
    "adaptative_threshold": true
}
```
### Running the Tool
After placing the `config.json` file in the same directory as the Python scripts, execute the following command:
` python main.py`

### Output
Percentile Threshold Files: NetCDF files containing computed thresholds for each scheme. ***can be deleted***
AR Tagging Files: NetCDF files with AR detection tags for each year and region processed.

---
## Contact

Victoire Buffet  
PhD student at Université Grenoble Alpes  
victoire.buffet@univ-grenoble-alpes.fr  

---
## References
Wille, J. D., Favier, V., Dufour, A., Gorodetskaya, I. V., Turner, J., Agosta, C., & Codron, F. (2019). West Antarctic surface melt triggered by atmospheric rivers. Nature Geoscience, 12(11), 911–916. https://doi.org/10.1038/s41561-019-0460-1  
Wille, J. D., Favier, V., Gorodetskaya, I. V., Agosta, C., Kittel, C., Beeman, J. C., et al. (2021). Antarctic atmospheric river climatology and precipitation impacts. Journal of Geophysical Research: Atmospheres, 126(8), e2020JD033788. https://doi.org/10.1029/2020JD033788
