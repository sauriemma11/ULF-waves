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
    * spect -- spectrogram that returns tx1 frequencies, tx1 times, tx1 spectrum
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
# TO DO: GET RID OF GLOBAL VARIABLES HERE
# THINGS TO ADD TO CONFIG: frequency/Pc band
#   -highpass or bandpass for tau filtering

def filter_b(b_mfa,ftype='highpass',comp=2,fs=10,N=1, fc=0.001):
    # TO DO: option for which components to get waves for?
    # TO DO: optoin for frequency band? dpdds on application..
    """
    Filter the mangetic field component of interest.
    input: b_mfa: magnetic field measurements in MFA coordinates, assumed to
                  be one time window, size nx3 [nT]
           ftype: 'high' or 'low' or 'bandpass', str
           comp: component of magnetic field to filter 0=radial, 1=phi, 2=paralell
           fs: sampling frequency, default 10Hz, int [Hz]
           N: filter order, deault 1, int
           fc: cutoff frequency IF highpass or lowpass filter, int
               if bandpass, 1x2 list of int with low and high frequencies

    output: highps_z: highpass filtered componenent of b_mfa. size 1xn [nT]
    """

    if ftype == 'highpass' or 'lowpass':
        b,a = u.butter_filter(fs ,fc,N, ftype)
    if ftype == 'bandpass':
        w1 = fs[0] / (fs/2) #normalize the frequency
        w2 = fs[1] / (fs/2)
        b,a = u.butter_filter(fs ,np.asarray([w1,w2]),N, ftype) #bandpass filter

    filt_comp = u.apply_butter(b_mfa[:][:,comp],b,a)

    return(filt_comp)

def spect(b_in, fs=10):
    """
    create spectrogram of magnetic field timeseries.
    input: b_in: magnetic field in MFA coordinates, filtered [nT]
           fs: sampling frequency, default is 10 Hz [Hz]
    output: f_sps: frequencies present in the spectrogram, [Hz]
            t_s: times in the spectrogram [s]
            Sxx_s: power in the spectrogram [nT^2/Hz]
    """
    NFFT = 4500

    f_sps, t_s, Sxx_s = signal.spectrogram(b_in,
            fs, nperseg=int(NFFT),noverlap = int(NFFT / 4),scaling='spectrum',
            return_onesided=True,mode='psd')
    return(f_sps,t_s,Sxx_s)

def get_f_band(f_hp_spectrum,fband):
    # TO DO: band option in config?
    """
    Find the indicies for the low and high frequencies

    """

    band = [fband[0],fband[1]] #band = frequency band limits,: 1-10mHz, use 0.011 so includes 10mHz
    delta_f = band[1]-band[0]
    #get start and stop index for frequency band
    band_st = u.find_nearest(f_hp_spectrum,band[0])
    band_sp = u.find_nearest(f_hp_spectrum,band[1])
    return(band_st,band_sp,delta_f)

def avg_psd(Sxx_hp_spectrum,  band_st, band_sp, delta_f, t_xx):

    band_pwrs = []
    for t in np.arange(0,len(t_xx)):
        band_sum = np.sum(Sxx_hp_spectrum[band_st:band_sp,t])
        band_pwr_tmp = band_sum / delta_f
        band_pwrs.append(band_pwr_tmp)
    psd_avg = np.mean(band_pwrs)
    return(psd_avg)

def calc_Dll(psd_avg):
    L = 6.6
    P_b = psd_avg#psd #psd_time #average power of field pertubations over all frequencies being used
                        #using paralell componenet bc leads to E accelerating electrons in radial directions
    f = (0.001 + 0.01) / 2#(0.0012 + 0.01) / 2 #central frequency, using middle of frequency band 1-10mHz
    B_e = 31200  #B field at equator : equitaorial magnetic field strength at surface of earth
                # = 0.312 G = 31200 nT
    Dll =( (L**8 * 4 * np.pi**2 ) / (9 * (8* B_e**2)) ) * (P_b * f**2)
    return(Dll)

def calc_tau(D_ll):

    tau_s = 1/D_ll
    tau_m = tau_s / 60
    return(tau_m)

def get_tau(b_mfa, fband,ftype,comp): #function to be called by main module
    """
    calculate tau.
    input: b_mfa: magnetic field values in MFA with background field
                   subtracted, tx1 [nT]
           time: datetime array, tx1 [s]
           fband: low and high frequency for the band of interest,
                  1x2, list of ints [Hz]
    output: tau: timescale to diffuse one Lshell, int [min]
            Dll: D_LL, int
            t_hp_spectrum: times corresponding to frequency domain,  [s]
            f_hp_spectrum: frequencies corresponding to frequency domain [Hz]
            Sxx_hp_spectrum: power spectrum corresponding to frequency domain [nT^2/Hz]
            b_filt: highpass (?) filtered componenet of magnetic field [nT]
    """

    b_filt = filter_b(b_mfa, ftype, comp)

    f_spect,t_spect,Sxx_spect = spect(b_filt)

    low_f, high_f, df = get_f_band(f_spect,fband)

    psd_av = avg_psd(Sxx_spect, low_f, high_f,df, t_spect)

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

    return(outs)
