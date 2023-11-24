'''
Module to calculate the tau value
-Inputs:
    * dictionary (output from mfa_transform.py) with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
-Outputs:
    * individual returns from functions, particularly
        tx1 tau
        tx1 D_LL
        tx1 PSD
        tx1 time
        tx1 frequencies
        tx1 filtered componenet of magnetic field

-Used in:
    * call_tau.py
    * main.py
-Functions:
    * highps_b -- highpass filter that returns tx1 B-field in MFA
    * spect -- spectrogram that returns tx1 frequencies, tx1 times, tx1
      spectrum
    * get_fband -- retrieves indexes for frequencies in band of interest and
                   returns start (low freq.) index and stop (high freq. index)
    * avg_psd -- returns average Power Spectral Density within freq. band
    * calc_dll -- calculate and return D_LL
    * calc_tau -- calculate and return tau
'''

import utils as u
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
# TODO: GET RID OF GLOBAL VARIABLES HERE
# THINGS TO ADD TO CONFIG: frequency/Pc band
#   -highpass or bandpass for tau filtering


def filter_b(b_mfa, ftype='highpass', comp=2, fs=10, N=1, fc=0.001):
    # TO DO: option for which components to get waves for?
    # TO DO: optoin for frequency band? dpdds on application..
    """
    Filter the mangetic field component of interest

    Parameters
    ----------
    b_mfa: size nx3 [nT]
        Magnetic field measurements in MFA coordinates, assumed to
        be one time window
    ftype: str
        Frequency type options are 'high' or 'low' or 'bandpass'
    comp: int
        Component of magnetic field to filter 0=radial, 1=phi, 2=paralell
    fs: int [Hz] (default = 10Hz)
        Sampling frequency
    N: int (default = 1)
        Filter order
    fc: int or 1x2 list
        int if cutoff frequency is highpass or lowpass filter
        1x2 list of int high and low frequencies if cutoff is bandpass filter

    Returns
    -------
    highps_z: n [nT]
        Highpass filtered componenent of b_mfa
    """

    if ftype == 'highpass' or 'lowpass':
        b, a = u.butter_filter(fs, fc, N, ftype)
    if ftype == 'bandpass':
        w1 = fs[0] / (fs/2)  # normalize the frequency
        w2 = fs[1] / (fs/2)
        b, a = u.butter_filter(fs, np.asarray([w1, w2]),
                               N, ftype)  # bandpass filter

    filt_comp = u.apply_butter(b_mfa[:][:, comp], b, a)  # noqa

    return (filt_comp)


def spect(b_in, fs=10):
    """
    Create spectrogram of magnetic field timeseries

    Parameters
    ----------
    b_in: [nT]
        magnetic field in MFA coordinates, filtered, size n (i.e.
        just one componenet of magnetic field)
    fs: int [Hz] (default = 10Hz)
        Sampling frequency

    Returns
    -------
    f_sps: [Hz]
        Frequencies present in the spectrogram, [Hz]
    t_s: [s]
        Times in the spectrogram
    Sxx_s: [nT^2/Hz]
        Power in the spectrogram
    """
    NFFT = 4500

    f_sps, t_s, Sxx_s = signal.spectrogram(b_in,
                                           fs, nperseg=int(NFFT),
                                           noverlap=int(NFFT / 4),
                                           scaling='spectrum',
                                           return_onesided=True, mode='psd')
    return (f_sps, t_s, Sxx_s)


# TO DO : ADD IN CONDITION FOR CHECKING FBAND[1] > FBAND[0]
def get_f_band(f_hp_spectrum, fband):
    # TO DO: band option in config?
    """
    Find the indicies for the low and high frequencies

    Parameters
    ----------
    inputs: f_hp_spectrum: frequencies in the spectrogram, mx1 [Hz]
            fband: low and high frequencies of band, 1x2 list i.e.
                   [low,high], [Hz]

    Returns
    -------
    band_st: int
        Index of the starting frequency
    band_sp: int
        Index of the stop frequency
    delta_f: int [Hz]
        Difference between the high and low frequencies
    """
    # TODO: repetitive
    # print('fband: ')
    # print(type(fband[0]), type(fband[1]))
    band = [fband[0], fband[1]]
    # print('band: ')
    # print(band, type(band[0]), type(band[1]))
    # band = frequency band limits,: 1-10mHz,
    # use 0.011 so includes 10mHz
    delta_f = float(band[1]-band[0])
    # get start and stop index for frequency band
    # TODO: THINK ABOUT WHAT THE TOLERANCE SHOULD BE HERE
    band_st = u.find_nearest(f_hp_spectrum, band[0], tolerance=100)
    band_sp = u.find_nearest(f_hp_spectrum, band[1], tolerance=100)
    return (band_st, band_sp, delta_f)


def avg_psd(Sxx_hp_spectrum,  band_st, band_sp, delta_f, t_xx):
    """
    Calculate the average power spectral density in the frequency band

    Parameters
    ----------
    Sxx_hp_spectrum: txm [nT^2/Hz]
        Power specrum from spectrogram
    band_st: int
        Starting index of frequencies in frequency band
    band_sp: int
        Stopping index of frequencies in frequency band
    delta_f: int [Hz]
        difference between the high and low frequencies
    t_xx: [s]
        Times from spectrum

    Returns
    -------
    psd_avg: int [nT^2/Hz]
        Average power spectral density
    """

    band_pwrs = []
    for t in np.arange(0, len(t_xx)):
        band_sum = np.sum(Sxx_hp_spectrum[band_st:band_sp, t])
        band_pwr_tmp = band_sum / delta_f
        band_pwrs.append(band_pwr_tmp)
    psd_avg = np.mean(band_pwrs)
    return (psd_avg)


def calc_Dll(psd_avg):
    # TO DO: MAKE F A VARIABLE / JUST PASS IN FBAND FROM INPUTS
    """
    Calculate D_LL

    Parameters
    ----------
    psd_avg: int [nT^2/Hz]
        Average power spectral density for the frequency band of interest

    Constants
    ---------
    L: Lshell (Earth Radius) [RE]
    B_e: equitorial magnetic field strength at surface of earth [nT]

    Outputs
    -------
    D_LL:

    """
    L = 6.6
    P_b = psd_avg
    f = (0.001 + 0.01) / 2  # noralizes
    B_e = 31200
    Dll = ((L**8 * 4 * np.pi**2) / (9 * (8*B_e**2))) * (P_b * f**2)
    return (Dll)


def calc_tau(D_ll):
    """
    Calculate tau (the time it takes for an electron to diffuse one Lshell)

    Parameters
    ----------
    D_ll:

    Returns
    -------
    tau: int [min]
    """
    tau_s = 1/D_ll
    tau_m = tau_s / 60
    return (tau_m)


def get_tau(b_mfa, fband=[0.001, 0.01], ftype='highpass', comp=2):
    # function to be called by main module
    """
    Calculate tau

    Parameters
    ----------
    b_mfa: tx1 [nT]
        Magnetic field values in MFA with background field subtracted,
    time: tx1 [s]
        Datetime array,
    fband: 1x2 list of ints [Hz]
        Low and high frequency for the band of interest
    comp: int
        Componenet of the magnetic field to filter, 0=radial, 1=phi, 2=paralell

    Returns
    -------
    outs: dict, containing
        tau: int [min]
            Timescale to diffuse one Lshell
        D_LL: int
            D_LL
        t_hp_spectrum: [s]
            Times corresponding to frequency domain
        f_hp_spectrum: [Hz]
            Frequencies corresponding to frequency domain
        Sxx_hp_spectrum: [nT^2/Hz]
            Power spectrum corresponding to frequency domain
        b_filt: [nT]
            highpass (?) filtered componenet of magnetic field
    """

    b_filt = filter_b(b_mfa, ftype, comp)

    f_spect, t_spect, Sxx_spect = spect(b_filt)

    low_f, high_f, df = get_f_band(f_spect, fband)

    psd_av = avg_psd(Sxx_spect, low_f, high_f, df, t_spect)

    D_LL = calc_Dll(psd_av)

    tau = calc_tau(D_LL)

    outs = {}

    outs['tau'] = tau
    outs['D_LL'] = D_LL
    outs['freqs'] = f_spect
    outs['time'] = t_spect
    outs['Sxx'] = Sxx_spect
    outs['psd'] = psd_av
    outs['b_filt'] = b_filt

    return (outs)
