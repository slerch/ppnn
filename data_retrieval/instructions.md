Note: The following has only been tested under Linux, might require additional work under Windows

## Prerequisites
- installation of Python and R
- installation of ecCodes library (https://software.ecmwf.int/wiki/display/ECC/ecCodes+Home)
- installation of ECMWF API to access public data sets (https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets), including registration of a user account

## Download ECMWF forecast data
- run data_retrieval/forecasts/retrieve_ecmwf_fc_data.py (requires ECMWF API)
- download grib data from http://apps.ecmwf.int/webmars/joblist/
- convert grib data to netCDF for better handling in Python and R:
    - navigate to folder with downloaded grib files
    - run data_retrieval/forecasts/convert_grib_to_nc.sh (requires ecCodes)