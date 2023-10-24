'''
Main module
-Inputs:
    * data file (format *.nc)
-Outputs:
    * plots
    * data file
'''

# Import statements here
# Or maybe we have a library file
import data_prep
import mfa_transform
import call_tau
import argparse

parser = argparse.ArgumentParser(
    description='Pass in parameters for calculating Tau for a given dataset.'
    )

parser.add_argument('--filename',
                    type=str,
                    help='Path to file of interest.',
                    required=True)

parser.add_argument('--timespan',
                    type=int,
                    default=1,  # default is 1 hour
                    help='Number of hours over which to find average Tau, PSD, and DLL values. Must be divisible by 24 (i.e. 1, 2, 3, 4, 6, 8, 12 hrs).',  #noqa
                    required=False)

parser.add_argument('--freq_band',
                    type=list,
                    default=[0.001, 0.01],  # default is 0.001 - 0.01 Hz
                    help='Values to define the frequency band of interest. First element is lower frequency, and second element is upper frequency [Hz].',  #noqa
                    required=False)

args = parser.parse_args()


def main(filename, timespan, freq_band):
    """
    Parameters
    ----------
    filename: str
        data file to work with

    timespan: integer
        hours to split up data by 

    freq_band: list of 2 integers
        first element is lower frequency, and second element is upper frequency

    Returns
    -------
    tau_dict: dict

    """
    # prep data
    variable_dict = data_prep.read_nc_file(filename) # dictionary of all variables

    # mfa_transform
    mfa_dict = mfa_transform.main(variable_dict)

    # calc Tau
    tau_dict = call_tau.call_tau(mfa_dict, timespan, freq_band) # calculates tau for each time window

    # CREATE PLOTS???

    return tau_dict


if __name__ == '__main__':
    main(args.filename, args.timespan, args.freq_band)
