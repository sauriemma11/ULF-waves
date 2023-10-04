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
from datetime import datetime, timedelta
import datetime
import numpy as np


def time_convert(seconds_2000):
    date_original = datetime(2000, 1, 1, 12, 0)
    return date_original + timedelta(seconds=int(seconds_2000))

#---- TODO : 1) make this part of an earlier step when we read in the l2 file the first time
#             2) make concatenate able to take in unlimited dimension
#note: this is only here because I was only reading in the pickle file and needed to re-read the L2 to get time


file = 'path/to/file.nc'
data = nc.Dataset(file)


t1 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc')
time_ob_native1 = t1['time']
time_ob_1 = np.asarray(time_ob_native1).flatten()

t2 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event/dn_magn-l2-hires_g16_d20230228_v1-0-1.nc')
time_ob_native2 = t2['time']
time_ob_2 = np.asarray(time_ob_native2).flatten()


time_g16 = np.concatenate((time_ob_1,time_ob_2),axis=None)