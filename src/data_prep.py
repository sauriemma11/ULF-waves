"""
Module to read in the data file and create outputs in the necessary format
- Inputs:
    * goes l2 hires data file (format *.nc)
- Outputs:
    * dictionary with
        ['b_epn'] nx3 array of B-field in EPN (with nan values from fill
        value -9999)
        ['time'] nx1 array of time in datetime formatted sec
- Used in:
    * main.py
- Functions:
    * read_nc_file -- read an nc file and extract its variables
    * read_variables -- read and store all variables from
                    an nc dataset into a dictionary
    * deal_with_fills -- replaces variables that are -9999
                    with NaN values in the 'b_epn' variable
    * time_convert -- convert time values from seconds
                    since 2000-01-01 to datetime objects
- Example usage:
    prepped_data = output_data_prepped_dict(
    '../data/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc')
    print(prepped_data)
"""

import netCDF4 as nc
from datetime import datetime, timedelta
import datetime
import numpy as np


def read_nc_file(file_path):
    """
    Read an nc file and extract its variables

    Parameters
    ----------
    file_path: str
        Full path to the nc file

    Returns
    -------
    variables: dict
        NetCDF dataset
    """
    dataset = nc.Dataset(file_path)
    variables = read_variables(dataset)
    dataset.close()
    return variables


def read_variables(dataset):
    """
    Read and store all variables from an nc dataset into a dictionary

    Parameters
    ----------
    dataset: NetCDF data
        NetCDF dataset

    Returns
    -------
    data_dict: dict
        Dictionary containing variables from the dataset
    """
    data_dict = {}
    for variable_name in dataset.variables.keys():
        data_dict[variable_name] = dataset.variables[variable_name][:]
    return data_dict


def deal_with_fills(data_dict):
    """
    Replaces variables that are -9999 with NaN values in the 'b_epn' variable.

    Parameters
    ----------
    data_dict: dict
        Dictionary containing data that needs to be filtered
        It has variables 'time' (nx1) and 'b_epn' (nx3)

    Returns
    -------
    new_epn: Dict
        A new dictionary with -9999 values in 'b_epn' replaced by NaNs.
    """
    b_epn = data_dict['b_epn']
    b_epn = np.where(b_epn == -9999, np.nan, b_epn)
    new_epn = {'b_epn': b_epn}
    return new_epn


def time_convert(seconds_2000):
    """
    Convert time values from seconds since 2000-01-01 to datetime objects

    Parameters
    ----------
    seconds_2000: numpy.ndarray
        Array of time values in seconds since 2000-01-01

    Returns
    -------
    dt_array: numpy.ndarray
        Array of datetime objects
    """
    date_original = datetime.datetime(2000, 1, 1, 12, 0)
    dt_array = np.array([date_original +
                         timedelta(seconds=int(seconds))
                         for seconds in seconds_2000])
    return dt_array


def output_data_prepped_dict(file_path):
    """
    Process a data file, including reading, dealing with fill values,
    and converting time. Puts everything into a new dictionary

    Parameters
    ----------
    file_path: str
        Full path to the data file

    Returns
    -------
    dict_prepped: dict
        Processed data dictionary containing 'time' and 'b_epn' with
        fill values corrected
    """

    original_dict = read_nc_file(file_path)
    epn_with_fills_fixed = deal_with_fills(original_dict)
    converted_time = time_convert(original_dict['time'])

    dict_prepped = {'time': converted_time,
                    'b_epn': epn_with_fills_fixed['b_epn']}

    return dict_prepped
