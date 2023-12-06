'''
Module to calculate the tau value
- Inputs:
    * dictionary (output from mfa_transform.py) with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
- Outputs:
    * individual returns from functions, particularly
        tx1 tau
        tx1 D_LL
        tx1 PSD
        tx1 time
        tx1 frequencies
        tx1 filtered componenet of magnetic field
- Used in:
    * call_tau.py
    * main.py
- Functions:
    * filter_b -- highpass filter that returns tx1 B-field in MFA
    * spect -- spectrogram that returns tx1 frequencies, tx1 times, tx1
      spectrum
    * get_fband_ind -- retrieve indexes for frequencies in band of interest and
                    returns start (low freq.) index and stop (high freq. index)
    * avg_psd -- returns average Power Spectral Density within freq. band
    * calc_Dll -- calculate and return D_LL
    * calc_tau -- calculate and return tau
    * get_tau -- calculate tau; the function to be called in call_tau.py
'''

import utils as u
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np


def filter_b(b_mfa, ftype='highpass', comp=2, fs=10, N=1, fc=0.001):
    """
    Filter the mangetic field component of interest

    Parameters
    ----------
    b_mfa: size nx3 [nT]
        Magnetic field measurements in MFA coordinates, assumed to
        be one time window
    ftype: str
        Frequency type options are 'highpass' or 'lowpass' or 'bandpass'
    comp: int
        Component of magnetic field to filter 0=radial, 1=phi, 2=paralell
    fs: int [Hz] (default = 10Hz)
        Sampling frequency
    N: int (default = 1)
        Filter order
    fc: float or 1x2 list
        float if cutoff frequency is highpass or lowpass filter
        1x2 list of int high and low frequencies if cutoff is bandpass filter

    Returns
    -------
    filt_comp: n [nT]
        Highpass filtered componenent of b_mfa
    """
    if str(ftype) == str('highpass') or str('lowpass'):
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


def get_fband_ind(f_hp_spectrum, fband):
    """
    Find the indicies for the low and high frequencies

    Parameters
    ----------
    f_hp_spectrum: mx1 [Hz]
        frequencies in the spectrogram,
    fband: 1x2 list of floats [Hz]
        low and high frequencies of band, [low freq, high freq]

    Returns
    -------
    band_st: int
        Index of the starting frequency
    band_sp: int
        Index of the stop frequency
    delta_f: int [Hz]
        Difference between the high and low frequencies
    """

    if fband[1] < fband[0] or fband[1] == fband[0]:
        print('Incorrect order for frequency band. should be [low, high]'
              'switching order.')
        fband_temp = np.copy(fband)
        fband[0] = fband_temp[1]
        fband[1] = fband_temp[0]

    delta_f = float(fband[1]-fband[0])
    band_st = u.find_nearest(f_hp_spectrum, fband[0], tolerance=0.005)
    band_sp = u.find_nearest(f_hp_spectrum, fband[1], tolerance=0.05)
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
    D_LL: float
        Diffusion coefficient

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
    D_ll: float
        Diffusion coefficient

    Returns
    -------
    tau_m: float [min]
        The time it takes for an electron to diffuse one Lshell
    """
    tau_s = 1/D_ll
    tau_m = tau_s / 60
    return (tau_m)


def get_tau(b_mfa, fband=[0.001, 0.01], ftype='highpass', comp=2):
    """
    Calculate tau -- function to be called in call_tau.py

    Parameters
    ----------
    b_mfa: tx1 [nT]
        Magnetic field values in MFA with background field subtracted,
    fband: 1x2 list of ints [Hz]
        Low and high frequency for the band of interest
    ftype: str
        Frequency type options are 'highpass' or 'lowpass' or 'bandpass'
    comp: int
        Componenet of the magnetic field to filter, 0=radial, 1=phi, 2=paralell

    Returns
    -------
    outs: dict, containing
        tau: int [min]
            Timescale to diffuse one Lshell
        D_LL: float
            Diffusion coefficient
        t_hp_spectrum: [s]
            Times corresponding to frequency domain
        f_hp_spectrum: [Hz]
            Frequencies corresponding to frequency domain
        Sxx_hp_spectrum: [nT^2/Hz]
            Power spectrum corresponding to frequency domain
        b_filt: [nT]
            Filtered componenet of magnetic field
    """

    b_filt = filter_b(b_mfa, ftype, comp)

    f_spect, t_spect, Sxx_spect = spect(b_filt)

    low_f, high_f, df = get_fband_ind(f_spect, fband)

    psd_av = avg_psd(Sxx_spect, low_f, high_f, df, t_spect)

    D_LL = calc_Dll(psd_av)

    tau = calc_tau(D_LL)
    outs = {}

    outs['tau'] = tau
    outs['D_LL'] = D_LL
    outs['psd'] = psd_av
    outs['b_filt'] = b_filt

    return (outs)
