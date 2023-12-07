'''
Module to find the Magnetic Field Aligned (MFA) coordinates
- Inputs:
    * dictionary (output from data_prep.py) with
        nx3 array of background-subtracted B-field in EPN (with nan values)
        nx1 array of time in dt sec
- Outputs:
    * dictionary with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
- Used in:
    * calc_tau.py
    * call_tau.py
    * main.py
- Functions:
    * get_bav -- calcualate average background field
    * compute_mfa -- perform EPN to MFA transformation
    * background_sub -- filter out (subtract) the background field from data
'''
import utils as u
import numpy as np


def get_bav(b_in, fs=10, fc=1/(30*60), N=2, btype='lowpass'):
    """
    Step 1: get average background field. This is done by taking a 30 min
            lowpass butterworth filter of the entire dataset. This calls
            butter_filter and apply butter.

    Parameters
    ----------
    b_epn: list of int nx3
        Magnetic field msmts in EPN coordinates where n = times
        and ea. col is a componenet.
        Data should be accounted for fill values already.
    b_time: list of int size n
        Timestamps corresponding to b_epn, converted to dt
    fs: [Hz] (10 Hz default)
        Sampling frequency
    fc: (30 minutes)
        Cut off frequency
    N : filter order
    btype: str
        Filter type, high/low/band pass (lowpass default)

    Outputs
    -------
    b_av : list of int size 1x3
        Average background field
    """
    b, a = u.butter_filter(fs, fc, N, btype)
    b_avx = u.apply_butter(b_in[:][:, 0], b, a)
    b_avy = u.apply_butter(b_in[:][:, 1], b, a)
    b_avz = u.apply_butter(b_in[:][:, 2], b, a)

    b_av = np.column_stack((b_avx, b_avy, b_avz))
    return (b_av)


def compute_mfa(b_av, b_epn):
    """
    Step 2: compute MFA coordinates. See wave_analysis_steps.pdf for details

    Parameters
    ----------
    b_av: float
        average background field from get_bav
    b_epn: list of int size nx3
        (Data prep already accoutned for) Magnetic field observations in EPN

    Outputs
    -------
    b_mfa: list of int size nx3
        Magnetic field measurements in MFA coordinates
    """
    b_mfa = np.zeros((len(b_epn), 3))
    for n in np.arange(len(b_epn)):
        e_par2 = b_av[n]/np.linalg.norm(b_av[n])
        e_phi2 = e_par2 * [-1, 0, 0]
        e_phi3 = e_phi2 / np.linalg.norm(e_phi2)
        e_r = e_phi2*e_par2
        e_r3 = e_r / np.linalg.norm(e_r)

        mfa_trans = np.row_stack((e_r, e_phi2, e_par2))
        b_mfa_temp = np.matmul(b_epn[n], mfa_trans.T)
        b_mfa[n] = b_mfa_temp

    return (b_mfa)


def background_sub(b_in):
    """
    Step 3: subtrack the background signal so that it is easier to see the
            waves in the frequency band of interest

    Parameters
    ----------
    b_in: list of ints size nx3
        Magnetic field measurements in MFA coordinates

    Outputs
    -------
    b_mfa_bsub: List of int size nx3
        Magnetic field measurements in MFA coordinates with
        the average background field subtracted

    Constants
    ---------
    fs: [Hz] (default = 10)
        Sampling frequency
    fc: [Hz] (30 min)
        Cutoff frequency
    btype: str
        Filter type (highpass)
    """

    fs = 10
    fc = 1/(30*60)
    N = 2
    btype = 'highpass'
    b, a = u.butter_filter(fs, fc, N, btype)
    b_x = u.apply_butter(b_in[:][:, 0], b, a)
    b_y = u.apply_butter(b_in[:][:, 1], b, a)
    b_z = u.apply_butter(b_in[:][:, 2], b, a)

    b_mfa_bsub = np.column_stack((b_x, b_y, b_z))
    return (b_mfa_bsub)


def main(b_epn):
    b_avg = get_bav(b_epn)

    b_mfa = compute_mfa(b_avg, b_epn)

    b_mga_bsub = background_sub(b_mfa)
    return (b_mga_bsub)
