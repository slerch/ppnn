Note: The following has only been tested under Linux, might require additional work under Windows

All data is available publicly from the TIGGE archives at ECMWF (requires registration, see http://apps.ecmwf.int/datasets/)

## Prerequisites
- installation of `Python` and `R`
- installation of `ecCodes` library (https://software.ecmwf.int/wiki/display/ECC/ecCodes+Home)
- installation of `ECMWF API` to access public data sets (https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets), including registration of a user account

## Download ECMWF forecast data
- run `data_retrieval/forecasts/retrieve_ecmwf_fc_data.py` (requires ECMWF API, may take minutes to hours for each of the 10 requested years)
- convert grib data to netCDF for better handling in Python and R:
    - navigate to folder with downloaded grib files (specified in `retrieve_ecmwf_fc_data.py`)
    - run `data_retrieval/forecasts/convert_grib_to_nc.sh` (requires ecCodes)
    
## Download DWD station observation data
- ...

## Interpolate forecasts to observation locations, and save data to netCDF file
- run `data_processing/interpolation.R`

## Download auxiliary model data
- run `data_retrieval/forecasts/retrieve_ecmwf_auxiliary_ZZZ_data.py`, with `ZZZ` replaced by
    - `geo` (to download land-sea mask and orography information)
    - `surface` (to download cloud cover, surface pressure and CAPE on the surface level)
    - `pl500` (to download u,v wind components and geopotential height at pressure level 500 hPa)
    - `pl850` (to download u,v wind components and specific humidity at pressure level 850 hPa)