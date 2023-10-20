'''
Module to calculate and combine tau values for specified time interval
-Inputs:
    * dictionary (output from mfa_transform.py) with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
    * user input for time interval (or a default configuation)
-Outputs:
    * individual returns from functions, particularly
        mx1 tau
        mx1 D_LL
        mx1 PSD
        mx1 time
        mx1 frequencies
-Used in:
    * main.py
-Functions:
    * define_window -- creates a window to split into smaller intervals
    * concat_tau -- gathers all B-field and time arrays into a dictionary
    * call_tau -- calculate tau for each time interval by looping over data set
'''

import calc_tau
import mfa_transform  #
import data_prep


def define_window(time_list, timespan_hrs=1):
    """
    Creates lists of indices each with a length matching input timespan

    Parameters
    ----------
    time_list: list
        List of datetime timesteps for when each data entry is received
        Typically 86400 secs/day with 10 entries/sec --> 864000 entries
    timespan_hrs: int (default = 1)
        Number of hours over which to find average Tau, PSD, and DLL
        Must be divisible by 24 (i.e. 1, 2, 3, 4, 6, 8, 12 hrs)

    Returns
    -------
    time_sublists: list
        List of lists; each list is indices of data entries for each timespan
    """
    # Find number of data entries over defined input (timespan in hours)
    num_timespan_data_entries = timespan_hrs * 3600 * 10
    time_sublists = []
    for i in range(0, len(time_list), num_timespan_data_entries):
        sublist = time_list[i:i + num_timespan_data_entries]
        time_sublists.append(sublist)
    return time_sublists


def concat_tau(bfields, time_list):
    """
    Concatenates individual data into one dictionary

    Parameters
    ----------


    Returns
    -------
    """
    for i in range(len(time_list)):
        bfield = bfields[i]
        time = time_list[i]
        tau_dict = 0
        return tau_dict


def call_tau(bfields, time_list):
    """
    Calls get_tau from calc_tau.py module to find average values over each
    window

    Parameters
    ----------


    Returns
    -------
    """
    windows = define_window(time_list, timespan_hrs=1)
    window_avg_tau = calc_tau.get_tau()

    concat_tau = []
    for time_interval in ?:
        one_tau = call_tau.get_tau()
        concat_tau.append(one_tau)
    return concat_tau
