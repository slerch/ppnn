for file in *.grib
do
  grib_to_netcdf -o "${file%.grib}.nc" $file
done
