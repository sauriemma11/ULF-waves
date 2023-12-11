'''
Module to calculate and combine tau values for specified time interval
- Requires:
    * function get_tau from calc_tau.py
    * user input for time interval (or a default configuation)
- Outputs:
    * dictionary with a list of lists for each key;
      each list of lists has 24 hours / defined timespan in hours entries;
      keys described in concat_tau function below
- Used in:
    * main.py
- Functions:
    * concat_tau -- gathers calculated values from each window into one dict
'''

import calc_tau
from tqdm import tqdm
import sys


def concat_tau(b_mfa, num_data_entries, fband, comp, ftype, timespan_hrs):
    """
    Concatenates individual window data into one dictionary

    Parameters
    ----------
    b_mfa: tx1 list [nT]
        Magnetic field values in MFA with background field subtracted
    num_data_entries: int
        Number of data entries in file; typically 864000, or
        10 samples/sec for 24 hours (86400 seconds)
    fband: 1x2 list of ints [Hz]
        Low and high frequency for the band of interest
    comp: int
        Componenet of the magnetic field to filter -- 0=radial,1=phi,2=paralell
    ftype: str
        Frequency type; options are 'high', 'low', or 'bandpass'

    timespan_hrs: int (default = 1) [Hrs]
        Number of hours over which to find average Tau, PSD, and DLL
        Must be a multiple of 24 (i.e. 1, 2, 3, 4, 6, 8, 12 hrs)
        Must be the same timespan used in define_windows function

    Returns
    -------
    all_windows_tau_dict: dict
        Dictionary with information from all windows concatenated
        Keys: "tau", "D_LL", "psd", "b_filt"
    """
    if 24 % timespan_hrs != 0:
        raise ValueError("Timespan must be a multiple of 24.")
        sys.exit(11)
    num_windows = int(24/timespan_hrs)
    len_one_window = int(num_data_entries/num_windows)

    all_windows_tau_dict = {"tau": [], "D_LL": [], "psd": [], "b_filt": []}

    dict_keys = ["tau", "D_LL", "psd"]
    for i in tqdm(range(num_windows)):
        b_mfa_window = b_mfa[i*len_one_window:(i+1)*len_one_window]
        tau_dict_for_window = calc_tau.get_tau(b_mfa_window, fband,
                                               ftype, comp)
        b_mfa_window_comp = [entry[comp] for entry in b_mfa_window]
        all_windows_tau_dict["b_filt"].append(b_mfa_window_comp)
        for key in dict_keys:
            all_windows_tau_dict[key].append(tau_dict_for_window[key])
    return all_windows_tau_dict
