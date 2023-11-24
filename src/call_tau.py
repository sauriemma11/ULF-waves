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
from tqdm import tqdm


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
        Dimension 24/timespan_hrs x 864000/(24/timespan_hrs)
    """
    # Check that 24 is divisible by entered timespan before proceeding
    if 24 % timespan_hrs != 0:
        raise ValueError("24 is not divisible by timespan used.")
        return None
    else:
        # Find number of data entries over defined input (timespan in hours)
        num_timespan_data_entries = timespan_hrs * 3600 * 10
        time_sublists = []
        for i in range(0, len(time_list), num_timespan_data_entries):
            sublist = time_list[i:i+num_timespan_data_entries]
            time_sublists.append(sublist)
        return time_sublists


def concat_tau(b_mfa, time_sublists,
               fband, comp=2, ftype='highpass', timespan_hrs=1):
    """
    Concatenates individual window data into one dictionary

    Parameters
    ----------
    b_mfa: tx1 list [nT]
        Magnetic field values in MFA with background field subtracted
    time_sublists: list
        List of lists; each list is indices of data entries for each timespan
        Dimension 24/timespan_hrs x 864000/(24/timespan_hrs)
    fband: 1x2 list of ints [Hz]
        Low and high frequency for the band of interest
    comp: int (default = 2)
        Componenet of the magnetic field to filter -- 0=radial,1=phi,2=paralell
    ftype: str
        Frequency type; options are 'high', 'low', or 'bandpass'

    timespan_hrs: int (default = 1) [Hrs]
        Number of hours over which to find average Tau, PSD, and DLL
        Must be divisible by 24 (i.e. 1, 2, 3, 4, 6, 8, 12 hrs)
        Must be the same timespan used in define_windows function

    Returns
    -------
    all_windows_tau_dict: dict
        Dictionary with information from all windows concatenated
        Keys: "tau", "D_LL", "psd", "Sxx", "time", "freqs", "b_filt"
    """

    num_windows = int(24/timespan_hrs)
    len_one_window = len(time_sublists[0])

    all_windows_tau_dict = {"tau": [], "D_LL": [], "psd": [], "Sxx": [],
                            "time": [], "freqs": [], "b_filt": []}

    dict_keys = ["tau", "D_LL", "psd", "Sxx", "time", "freqs", "b_filt"]

    for i in tqdm(range(num_windows)):
        b_mfa_window = b_mfa[i*len_one_window:(i+1)*len_one_window]
        fband_window = fband[i*len_one_window:(i+1)*len_one_window]

        tau_dict_for_window = calc_tau.get_tau(b_mfa_window, fband_window,
                                               ftype, comp)
        for key in dict_keys:
            all_windows_tau_dict[key].append(tau_dict_for_window[key])
    return all_windows_tau_dict
