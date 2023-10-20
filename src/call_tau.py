'''
Module to calculate and combine tau values for specified time interval
-Requires:
    * function get_tau from calc_tau.py
    * user input for time interval (or a default configuation)
-Outputs:
    * dictionary with a list of lists for each key;
      each list of lists has 24 hours / defined timespan in hours entries;
      keys described in concat_tau function below
-Used in:
    * main.py
-Functions:
    * define_window -- creates a window to split into smaller intervals
    * concat_tau -- gathers calculated values from each window into one dict
'''

import calc_tau


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
        sublist = time_list[i:i+num_timespan_data_entries]
        time_sublists.append(sublist)
    return time_sublists


def concat_tau(b_mfa, fband, ftype, comp, timespan_hrs=1):
    """
    Concatenates individual window data into one dictionary

    Parameters
    ----------
    b_mfa: tx1 list [nT]
        Magnetic field values in MFA with background field subtracted
    time: tx1 list [s]
        Datetime array
    fband: 1x2 list of ints [Hz]
        Low and high frequency for the band of interest
    comp:
        Component of magnetic field
    timespan_hrs: int (default = 1) [Hrs]
        Number of hours over which to find average Tau, PSD, and DLL
        Must be divisible by 24 (i.e. 1, 2, 3, 4, 6, 8, 12 hrs)

    Returns
    -------
    all_windows_tau_dict: dict
        Dictionary with information from all windows concatenated
        Keys: "tau", "D_LL", "psd", "Sxx", "time", "freqs", "b_filt"
    """
    num_windows = 24/timespan_hrs
    all_windows_tau_dict = {"tau": [], "D_LL": [], "psd": [], "Sxx": [],
                            "time": [], "freqs": [], "b_filt": []}
    for i in range(num_windows):
        tau_dict_for_window = calc_tau.get_tau(b_mfa, fband, ftype, comp)
        all_windows_tau_dict["tau"].append(tau_dict_for_window["tau"])
        all_windows_tau_dict["D_LL"].append(tau_dict_for_window["D_LL"])
        all_windows_tau_dict["psd"].append(tau_dict_for_window["psd"])
        all_windows_tau_dict["Sxx"].append(tau_dict_for_window["Sxx"])
        all_windows_tau_dict["time"].append(tau_dict_for_window["time"])
        all_windows_tau_dict["freqs"].append(tau_dict_for_window["freqs"])
        all_windows_tau_dict["b_filt"].append(tau_dict_for_window["freqs"])
    return all_windows_tau_dict
