'''
Main module
-Inputs:
    * data file (format *.nc)
-Outputs:
    * plots
    * data file
'''

import data_prep
import mfa_transform
import call_tau
import argparse
import os.path
import sys
import plots
import pickle
import numpy as np
from itertools import chain

parser = argparse.ArgumentParser(
    description="""Pass in parameters for calculating
      Tau for a given dataset.""")

parser.add_argument('--filename',
                    type=str,
                    help="""Path to file of interest.
                    Must be '.nc' file type.""",
                    required=True)

parser.add_argument('--timespan',
                    type=int,
                    default=1,  # default is 1 hour
                    help="""Number of hours over which to find average
                         Tau, PSD, and DLL values. Must be divisible by 24
                         (i.e. 1, 2, 3, 4, 6, 8, 12 hrs).""",
                    required=False)

parser.add_argument('--num_entries',
                    type=int,
                    # default is 86400 secs with 10 samples/sec = 864000
                    default=864000,
                    help="""Number of entries in data file.
                         Typically 864000 (10 ssamples/sec for 24 hours)""",
                    required=False)

parser.add_argument('--fband',
                    nargs='+',
                    type=float,
                    default=[0.001, 0.01],  # default is 0.001 - 0.01 Hz
                    help="""Values to define the frequency band
                         of interest. First element is lower frequency,
                         and second element is upper frequency [Hz].""",
                    required=False)

parser.add_argument('--comp',
                    choices=[0, 1, 2],
                    default=2,  # default is 2 -- parallel
                    help="""Component of the magnetic field to filter
                         (0=radial, 1=phi, 2=parallel)""",
                    required=False)

parser.add_argument('--ftype',
                    choices=['highpass', 'lowpass', 'bandpass'],
                    default='highpass',  # default is highpass filter
                    help="""Frequency type;
                    (options are 'high', 'low', or 'bandpass')""",
                    required=False)

args = parser.parse_args()


def main(filename, timespan, num_entries, fband, comp, ftype):
    """
    Parameters
    ----------
    filename: str
      data file to work with

    timespan: integer
    `hours to split up data by

    freq_band: list of 2 integers
      first element is lower frequency,
      and second element is upper frequency

    Returns
    -------
    tau_dict: dict
    """
    # checking that file exists -- 0
    check_file = os.path.exists(filename)
    if not check_file:
        print(f'Could not find {filename}. Try a different file.')
        sys.exit(1)

    # check that file is '.nc' -- 0
    is_nc = filename.endswith('.nc')
    if not is_nc:
        print("Wrong file type -- must be '.nc' file.")
        sys.exit(2)

    # checking `timespan` input is factor of 24 -- 1
    if not 24 % timespan == 0:
        print("'timespan' must be a factor of 24.")
        sys.exit(11)

    # checking `fband` length (must have 2 elements) -- 3
    if len(fband) != 2:
        print(f"""`fband` must have 2 elements
              but given `fband` has length {len(fband)}.""")
        sys.exit(31)

    # if `fband` elements are reversed, flip it
    if fband[0] > fband[1]:
        print('Flipping `fband` elements because given in wrong order')
        fband = [fband[1], fband[0]]

    print('Prepping data...')
    variable_dict = data_prep.read_nc_file(filename)

    print('Performing MFA transform...')
    b_mga = mfa_transform.main(variable_dict['b_epn'])

    print('Calculating tau...')
    tau_dict = call_tau.concat_tau(b_mga,
                                   num_data_entries=num_entries,
                                   fband=fband,
                                   comp=comp,
                                   ftype=ftype,
                                   timespan_hrs=timespan)

    # Save tau_dict as pickle:
    file_path = "../docs/tau_dict.pkl"
    # file_path = 'tau_dict.pkl'
    with open(file_path, 'wb') as file:
        pickle.dump(tau_dict, file)
    print(f'tau_dict saved to {file}')

    number_of_windows = int(24/timespan)
    window_start_times = np.linspace(0, num_entries, number_of_windows)

    # chain from itertools flattens the list
    magnetic_field_data = list(chain(*tau_dict['b_filt']))

    print('Plotting data...')
    plots.plot_data(variable_dict['time'],
                    magnetic_field_data,
                    window_start_times,
                    tau_dict['psd'],
                    tau_dict['tau'],
                    output_dir='../docs/output_plot.png')

    return tau_dict


if __name__ == '__main__':
    main(args.filename,
         args.timespan,
         args.num_entries,
         args.fband,
         args.comp,
         args.ftype)
