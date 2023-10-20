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


def main(file_path, timespan, freq_band):
    """
    Parameters
    ----------
    file_path: str
        data file to work with

    timespan: integer
        hours to split up data by 

    freq_band: list of 2 integers
        first element is lower frequency
        second element is upper frequency

    Returns
    -------
    tau_dict: dict

    """
    # prep data
    variable_dict = data_prep.read_nc_file(
        file_path)  # dictionary of all variables

    # mfa_transform
    mfa_dict = mfa_transform.main(variable_dict)

    # calc Tau
    tau_dict = call_tau.call_tau(mfa_dict, timespan,
                                 freq_band)  # calculates tau for each time
    # window

    # CREATE PLOTS???

    return tau_dict


if __name__ == '__main__':
# Do the main stuff
# main()
