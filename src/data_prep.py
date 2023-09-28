"""
This module will:

- Read .nc files
- Deal with fill values
- Convert times (convert_time(sec->dt))
- Take user input?? (maybe should be in main)

Outputs:
- Array / dictonary of 'prepped' data (IE fill values corrected, in EPN)
 N x 3 array
- Time N x 1 array
- User input???
"""
import netCDF4 as nc


file = 'path/to/file.nc'
data = nc.Dataset(file)

