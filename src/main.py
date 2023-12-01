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
import os.path
import sys

parser = argparse.ArgumentParser(
    description="""Pass in parameters for calculating
      Tau for a given dataset."""
)

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
                    default= 864000,
                    help="""Number of entries in data file.
                         Typically 864000 (10 ssamples/sec for 24 hours)""",
                    required=False)

parser.add_argument('--fband',
                    type=list,
                    default=[0.001, 0.01],  # default is 0.001 - 0.01 Hz
                    help="""Values to define the frequency band
                         of interest.First element is lower frequency,
                         and second element is upper frequency [Hz].""",
                    required=False)

parser.add_argument('--comp',
                    type=list,
                    default=2,  # default is 2 --- parallel
                    help="""Componenet of the magnetic field to filter
                         (0=radial, 1=phi, 2=parallel)""",
                    required=False)

parser.add_argument('--ftype',
                    type=str,
                    default='highpass',
                    help="""Frequency type;
                         options are 'high', 'low', or 'bandpass'""",
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
    # checking that file exists
    check_file = os.path.exists(filename)
    if not check_file:
        print(f'Could not find {filename}. Try a different file.')
        sys.exit(1)

    # check that file is '.nc'
    is_nc = filename.endswith('.nc')
    if not is_nc:
        print("Wrong file type -- must be '.nc' file.")
        sys.exit(2)

    # checking 'timespan' input
    if not 24 % timespan:
        print("'timespan' must be a factor of 24.")
        sys.exit(21)


    variable_dict = data_prep.read_nc_file(filename)

    b_mga = mfa_transform.main(variable_dict['b_epn'])

    tau_dict = call_tau.concat_tau(b_mga,
                                   num_data_entries=num_entries,
                                   fband=fband,
                                   comp=comp,
                                   ftype=ftype,
                                   timespan_hrs=timespan)

    # CREATE PLOTS???

    return tau_dict


if __name__ == '__main__':
    try:
        main(args.filename,
            args.timespan,
            args.num_entries,
            args.fband,
            args.comp,
            args.ftype)
    except FileNotFoundError:
      print("The file does not exist. Please check your path.")
