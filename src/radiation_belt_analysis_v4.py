#!/usr/bin/env python
# coding: utf-8

# ---- TODO : move imports and functions to a library file

##from utils import *
import os
import csv
from scipy import signal
import numpy as np
from datetime import datetime, timedelta
import datetime
from scipy.interpolate import interp1d
from typing import Dict, Optional, List
import deepdish as dd
from glob import glob
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from MAG_PLT_003 import *
##import spacepy.time as sptime
##import spacepy.coordinates as spcoords
import datetime as dt
import pickle

'''
legacy is easier_to_run_radiation_belt_analys.ipynb
updates to hopefully run even better and with updated Dll from Osmane 2023
'''


def time_convert(seconds_2000):
    date_original = datetime(2000, 1, 1, 12, 0)
    return date_original + timedelta(seconds=int(seconds_2000))


def find_nearest(array, value):
    """"
    Finds closest frequency and outputs the closest frequency. Workaround
    for how python stores decimals, there isn't going to have an exact match
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def getTau(chunk, b_filt_nT_mfai):
    b_hp_filt_short = b_filt_nT_mfai[chunk[0]:chunk[1]]
    fs = 10  # sampling frequncy, Hz, usually 10 Hz for mag
    fc = 0.001  # cut off frequency 10mHz = 0.01 Hz, 1mHz = 0.001 Hz
    w1 = 0.001 / (fs / 2)  # normalize the frequency
    w2 = 0.01 / (fs / 2)
    N = 1  # filter order
    # btype = 'bandpass'# 'high' #high or low pass filterbb\\
    btype = 'high'

    # b,a = butter_filter(fs ,np.asarray([w1,w2]),N, btype) #bandpass filter
    b, a = butter_filter(fs, fc, N, btype)
    highps_z = apply_butter(b_hp_filt_short[:][:, 2], b, a)

    NFFT = 4500

    f_hp_spectrum, t_hp_spectrum, Sxx_hp_spectrum = signal.spectrogram(
        highps_z, fs=10, nperseg=int(NFFT), noverlap=int(NFFT / 4),
        scaling='spectrum', return_onesided=True, mode='psd')

    f_hp_z, t_hp_z, Sxx_hp_z = signal.spectrogram(highps_z, fs=10,
                                                  nperseg=int(NFFT),
                                                  noverlap=int(NFFT / 4))
    sgdb_hp_z = 10 * np.log10(Sxx_hp_z)

    fig, ax = plt.subplots(figsize=(12, 6))
    cmap = plt.get_cmap('rainbow')
    im = plt.pcolormesh(t_hp_z, f_hp_z, sgdb_hp_z, vmin=-20, vmax=20,
                        cmap=cmap)
    formatter = matplotlib.ticker.FuncFormatter(timeTicks)
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.ylim([0, 0.05])
    fig.colorbar(im, ax=ax, label='Power Spectral Density (dB)')

    band = [0.001,
            0.01]  # band = frequency band limits,: 1-10mHz, use 0.011 so
    # includes 10mHz
    delta_f = band[1] - band[0]
    # get start and stop index for frequency band
    band_st = find_nearest(f_hp_spectrum, band[0])
    band_sp = find_nearest(f_hp_spectrum, band[1])

    band_pwrs = []
    for t in np.arange(0, len(t_hp_z)):
        band_sum = np.sum(Sxx_hp_spectrum[band_st:band_sp, t])
        band_pwr_tmp = band_sum / delta_f
        band_pwrs.append(band_pwr_tmp)
    psd_avg = np.mean(band_pwrs)
    print(psd_avg)
    # new Dll, updated from how we interperet Osmane 2023

    L = 6.6
    P_b = psd_avg  # psd #psd_time #average power of field pertubations over
    # all frequencies being used
    # using paralell componenet bc leads to E accelerating electrons in
    # radial directions
    f = (
                    0.001 + 0.01) / 2  # (0.0012 + 0.01) / 2 #central
    # frequency, using middle of frequency band 1-10mHz
    B_e = 31200  # B field at equator : equitaorial magnetic field strength
    # at surface of earth
    # = 0.312 G = 31200 nT
    Dll = ((L ** 8 * 4 * np.pi ** 2) / (9 * (8 * B_e ** 2))) * (P_b * f ** 2)
    print('Sandhu', Dll)

    tau_s = 1 / Dll
    tau_m = tau_s / 60
    print('Sandhu tau : ', tau_m)

    return (highps_z, f_hp_z, t_hp_z, sgdb_hp_z, psd_avg, Dll, tau_m)


def butter_filter(fs, fc, N, btype):
    # fs = sampling frequency [Hz], fc = cut off frequency [Hz], N = filter
    # order, btype = high / low
    w = fc / (fs / 2)  # normalize the frequency
    b, a = signal.butter(N, w, btype)  # design filter
    return b, a


def apply_butter(siggy, b, a):
    # b,a = output from butter_filter, signal = 1d array
    filtered = signal.filtfilt(b, a, siggy[
        ~np.isnan(siggy)])  # apply filter forwards and back, ignore nans
    return filtered


def parse_10Hz_L2(file, prod='b_epn'):
    data = nc.Dataset(file)
    ob_native = data[prod][:]

    ob_x = np.asarray(ob_native[:, 0]).flatten()
    ob_y = np.asarray(ob_native[:, 1]).flatten()
    ob_z = np.asarray(ob_native[:, 2]).flatten()

    time_ob_native = data['time']
    time_ob_flat = np.asarray(time_ob_native).flatten()

    ob_full = np.column_stack((ob_x, ob_y, ob_z))

    return (ob_full, time_ob_flat)


def get_bav(b_in):
    fs = 10  # sampling frequncy, Hz, usually 10 Hz for mag
    fc = 1 / (30 * 60)  # cut off frequency, 30 min
    N = 2  # filter order
    btype = 'lowpass'  # high or low pass filter
    b, a = butter_filter(fs, fc, N, btype)
    b_avx = apply_butter(b_in[:][:, 0], b, a)
    b_avy = apply_butter(b_in[:][:, 1], b, a)
    b_avz = apply_butter(b_in[:][:, 2], b, a)

    b_av = np.column_stack((b_avx, b_avy, b_avz))
    return (b_av)


def background_sub(b_in):
    fs = 10  # sampling frequncy, Hz, usually 10 Hz for mag
    fc = 1 / (30 * 60)  # cut off frequency, 30 min
    N = 2  # filter order
    btype = 'high'  # high or low pass filter
    b, a = butter_filter(fs, fc, N, btype)
    b_x = apply_butter(b_in[:][:, 0], b,
                       a)  # applying highpass filter subtracts the
    # background field(lower freq)
    b_y = apply_butter(b_in[:][:, 1], b,
                       a)  # same as subtracting the low pass filter
    b_z = apply_butter(b_in[:][:, 2], b, a)

    b_backsub = np.column_stack((b_x, b_y, b_z))
    return (b_backsub)


def get_bav_bsub(bfile, date):
    # may have to ammend this based on what the fill values are like...
    g_epn, g_time = parse_10Hz_L2(bfile, 'b_epn')
    # convert times to datetimes
    g_dt = []
    for i, time_new in enumerate(g_time.data):
        g_dt.append(time_convert(time_new))
        # set fill values to nan
    g_epn = np.where((g_epn != -9.999e+03), g_epn,
                     np.nan)  # should work for nx3
    # get average background field
    b_av = get_bav(g_epn)

    # calc mfa coordinates
    b_mfa = np.zeros((len(g_epn), 3))
    for n in np.arange(len(g_epn)):
        e_par2 = b_av[n] / np.linalg.norm(b_av[n])
        e_phi2 = e_par2 * [-1, 0,
                           0]  # e_phi = -1 points radially away from earth
        e_phi3 = e_phi2 / np.linalg.norm(e_phi2)
        e_r = e_phi2 * e_par2
        e_r3 = e_r / np.linalg.norm(e_r)

        mfa_trans = np.row_stack(
            (e_r, e_phi2, e_par2))  # create transformation matrix
        b_mfa_temp = np.matmul(g_epn[n], mfa_trans.T)  # rotate EPN to MFA

        b_mfa[n] = b_mfa_temp
        # subtract the background field
    b_hp_filt_nT_mfa_o = background_sub(b_mfa)
    # save for later in a pickle
    waves_out = {}
    waves_out['b_hp_filt_nT_mfa'] = b_hp_filt_nT_mfa_o
    waves_out['mfa'] = b_mfa
    waves_out['epn'] = g_epn
    waves_out['times'] = g_dt

    name = '/Users/aspen.davis/Documents/radiation_belt/mfa_bav_' + date + \
           '.pickle'
    f = open(name, 'wb')
    pickle.dump(waves_out, f)
    f.close()
    return (b_hp_filt_nT_mfa_o)


'''
get field in MFA coordinates and subtract background field
'''
# confirmed, works
# ---- TODO : 1) break get_bav_bsub into smaller functions
#             2) do not create a file in get_bav_bsub, just return array for 
#                b_hp_filt_nT_mfa
#             3) maybe accept a command line entry for the file? possible
#             accept mulitple files?
date = '20230227'
bfile = '/Users/aspen.davis/Desktop/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc'
b_hp_filt_nT_mfa = get_bav_bsub(bfile, date)


# ---- TODO : move this to library with all functions
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


'''
call getTau for each 3 hour "chunk" portion of b_hp_filt_nt_mfa
'''

# ---- TODO : 1) user input for time length (chunk) to calculate Tau for
#             2) remove the steps to open and load the pickle (or keep if we
#             decide thats best)
#             3) create easier return type rather than a bunch or arrays (
#             dict? named tuple?)
#             4) break getTau into smaller functions 
#             5) create a loop or a calling function to call getTau rather
#             than repeated line calls

date = '20230227'
# psd_store1 = '/Users/aspen.davis/Documents/radiation_belt/20220809
# /mfa_bav_g17_20220808.pickle'
# psd_store1 = '/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/mfa_bav_g18_20230227.pickle'
psd_store1 = '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event' \
             '/mfa_bav_g16_20230227.pickle'
file1 = open(psd_store1, 'rb')
psd_analysis1 = pickle.load(file1)
b_hp_filt_nT_mfa1 = psd_analysis1['b_hp_filt_nT_mfa']

highps_z1_01, f_hp_z1_01, t_hp_z1_01, sgdb_hp_z1_01, psd_avg1_01, Dll1_01, \
    tau_m1_01 = getTau(
    chunk_1, b_hp_filt_nT_mfa1)
highps_z2_01, f_hp_z2_01, t_hp_z2_01, sgdb_hp_z2_01, psd_avg2_01, Dll2_01, \
    tau_m2_01 = getTau(
    chunk_2, b_hp_filt_nT_mfa1)
highps_z3_01, f_hp_z3_01, t_hp_z3_01, sgdb_hp_z3_01, psd_avg3_01, Dll3_01, \
    tau_m3_01 = getTau(
    chunk_3, b_hp_filt_nT_mfa1)
highps_z4_01, f_hp_z4_01, t_hp_z4_01, sgdb_hp_z4_01, psd_avg4_01, Dll4_01, \
    tau_m4_01 = getTau(
    chunk_4, b_hp_filt_nT_mfa1)
highps_z5_01, f_hp_z5_01, t_hp_z5_01, sgdb_hp_z5_01, psd_avg5_01, Dll5_01, \
    tau_m5_01 = getTau(
    chunk_5, b_hp_filt_nT_mfa1)
highps_z6_01, f_hp_z6_01, t_hp_z6_01, sgdb_hp_z6_01, psd_avg6_01, Dll6_01, \
    tau_m6_01 = getTau(
    chunk_6, b_hp_filt_nT_mfa1)
highps_z7_01, f_hp_z7_01, t_hp_z7_01, sgdb_hp_z7_01, psd_avg7_01, Dll7_01, \
    tau_m7_01 = getTau(
    chunk_7, b_hp_filt_nT_mfa1)
highps_z8_01, f_hp_z8_01, t_hp_z8_01, sgdb_hp_z8_01, psd_avg8_01, Dll8_01, \
    tau_m8_01 = getTau(
    chunk_8, b_hp_filt_nT_mfa1)

date = '20230228'
# psd_store2 = '/Users/aspen.davis/Documents/radiation_belt/20220809
# /mfa_bav_g17_20220809.pickle'
# psd_store2 = '/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/mfa_bav_g18_20230228.pickle'
psd_store2 = '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event' \
             '/mfa_bav_g16_20230228.pickle'
file2 = open(psd_store2, 'rb')
psd_analysis2 = pickle.load(file2)
b_hp_filt_nT_mfa2 = psd_analysis2['b_hp_filt_nT_mfa']
highps_z1_02, f_hp_z1_02, t_hp_z1_02, sgdb_hp_z1_02, psd_avg1_02, Dll1_02, \
    tau_m1_02 = getTau(
    chunk_1, b_hp_filt_nT_mfa2)
highps_z2_02, f_hp_z2_02, t_hp_z2_02, sgdb_hp_z2_02, psd_avg2_02, Dll2_02, \
    tau_m2_02 = getTau(
    chunk_2, b_hp_filt_nT_mfa2)
highps_z3_02, f_hp_z3_02, t_hp_z3_02, sgdb_hp_z3_02, psd_avg3_02, Dll3_02, \
    tau_m3_02 = getTau(
    chunk_3, b_hp_filt_nT_mfa2)
highps_z4_02, f_hp_z4_02, t_hp_z4_02, sgdb_hp_z4_02, psd_avg4_02, Dll4_02, \
    tau_m4_02 = getTau(
    chunk_4, b_hp_filt_nT_mfa2)

highps_z5_02, f_hp_z5_02, t_hp_z5_02, sgdb_hp_z5_02, psd_avg5_02, Dll5_02, \
    tau_m5_02 = getTau(
    chunk_5, b_hp_filt_nT_mfa2)
highps_z6_02, f_hp_z6_02, t_hp_z6_02, sgdb_hp_z6_02, psd_avg6_02, Dll6_02, \
    tau_m6_02 = getTau(
    chunk_6, b_hp_filt_nT_mfa2)
highps_z7_02, f_hp_z7_02, t_hp_z7_02, sgdb_hp_z7_02, psd_avg7_02, Dll7_02, \
    tau_m7_02 = getTau(
    chunk_7, b_hp_filt_nT_mfa2)
highps_z8_02, f_hp_z8_02, t_hp_z8_02, sgdb_hp_z8_02, psd_avg8_02, Dll8_02, \
    tau_m8_02 = getTau(
    chunk_8, b_hp_filt_nT_mfa2)

date = '20230228'
# psd_store3 = '/Users/aspen.davis/Documents/radiation_belt/20220809
# /mfa_bav_g17_20220810.pickle'
# psd_store3 = '/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/mfa_bav_g18_20230301.pickle'
psd_store3 = '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event' \
             '/mfa_bav_g16_20230301.pickle'
file3 = open(psd_store3, 'rb')
psd_analysis3 = pickle.load(file3)
b_hp_filt_nT_mfa3 = psd_analysis3['b_hp_filt_nT_mfa']

highps_z1_03, f_hp_z1_03, t_hp_z1_03, sgdb_hp_z1_03, psd_avg1_03, Dll1_03, \
    tau_m1_03 = getTau(
    chunk_1, b_hp_filt_nT_mfa3)
highps_z2_03, f_hp_z2_03, t_hp_z2_03, sgdb_hp_z2_03, psd_avg2_03, Dll2_03, \
    tau_m2_03 = getTau(
    chunk_2, b_hp_filt_nT_mfa3)
highps_z3_03, f_hp_z3_03, t_hp_z3_03, sgdb_hp_z3_03, psd_avg3_03, Dll3_03, \
    tau_m3_03 = getTau(
    chunk_3, b_hp_filt_nT_mfa3)
highps_z4_03, f_hp_z4_03, t_hp_z4_03, sgdb_hp_z4_03, psd_avg4_03, Dll4_03, \
    tau_m4_03 = getTau(
    chunk_4, b_hp_filt_nT_mfa3)

highps_z5_03, f_hp_z5_03, t_hp_z5_03, sgdb_hp_z5_03, psd_avg5_03, Dll5_03, \
    tau_m5_03 = getTau(
    chunk_5, b_hp_filt_nT_mfa3)
highps_z6_03, f_hp_z6_03, t_hp_z6_03, sgdb_hp_z6_03, psd_avg6_03, Dll6_03, \
    tau_m6_03 = getTau(
    chunk_6, b_hp_filt_nT_mfa3)
highps_z7_03, f_hp_z7_03, t_hp_z7_03, sgdb_hp_z7_03, psd_avg7_03, Dll7_03, \
    tau_m7_03 = getTau(
    chunk_7, b_hp_filt_nT_mfa3)
highps_z8_03, f_hp_z8_03, t_hp_z8_03, sgdb_hp_z8_03, psd_avg8_03, Dll8_03, \
    tau_m8_03 = getTau(
    chunk_8, b_hp_filt_nT_mfa3)

'''
concatenate all of the individual tau / highpass filtered data, etc. 

'''
# ---- TODO : do this better lol. allow for unlimited dimensions based on
#  days / "chunl" lengths

highps_z_01 = np.concatenate((highps_z1_01, highps_z2_01, highps_z3_01,
                              highps_z4_01, highps_z5_01, highps_z6_01,
                              highps_z7_01, highps_z8_01), axis=None)
highps_z_02 = np.concatenate((highps_z1_02, highps_z2_02, highps_z3_02,
                              highps_z4_02, highps_z5_02, highps_z6_02,
                              highps_z7_02, highps_z8_02), axis=None)
highps_z_03 = np.concatenate((highps_z1_03, highps_z2_03, highps_z3_03,
                              highps_z4_03, highps_z5_03, highps_z6_03,
                              highps_z7_03, highps_z8_03), axis=None)
# highps_z_04 = np.concatenate((highps_z1_04, highps_z2_04,highps_z3_04,
# highps_z4_04,highps_z5_04,highps_z6_04,highps_z7_04,highps_z8_04), axis=None)
# highps_z_all = np.concatenate((highps_z_01,highps_z_02,highps_z_03,
# highps_z_04),axis=None)
highps_z_all = np.concatenate((highps_z_01, highps_z_02, highps_z_03),
                              axis=None)

tau_m_02 = np.concatenate((tau_m1_02, tau_m2_02, tau_m3_02, tau_m4_02,
                           tau_m5_02, tau_m6_02, tau_m7_02, tau_m8_02),
                          axis=None)
tau_m_03 = np.concatenate((tau_m1_03, tau_m2_03, tau_m3_03, tau_m4_03,
                           tau_m5_03, tau_m6_03, tau_m7_03, tau_m8_03),
                          axis=None)
# tau_m_04 = np.concatenate((tau_m1_04,tau_m2_04,tau_m3_04,tau_m4_04,
# tau_m5_04,tau_m6_04,tau_m7_04,tau_m8_04),axis=None)
tau_m_all = np.concatenate((tau_m_01, tau_m_02, tau_m_03), axis=None)

psd_avg_01 = np.concatenate((
                            psd_avg1_01, psd_avg2_01, psd_avg3_01, psd_avg4_01,
                            psd_avg5_01, psd_avg6_01, psd_avg7_01,
                            psd_avg8_01), axis=None)
psd_avg_02 = np.concatenate((
                            psd_avg1_02, psd_avg2_02, psd_avg3_02, psd_avg4_02,
                            psd_avg5_02, psd_avg6_02, psd_avg7_02,
                            psd_avg8_02), axis=None)
psd_avg_03 = np.concatenate((
                            psd_avg1_03, psd_avg2_03, psd_avg3_03, psd_avg4_03,
                            psd_avg5_03, psd_avg6_03, psd_avg7_03,
                            psd_avg8_03), axis=None)
# psd_avg_04 = np.concatenate((psd_avg1_04,psd_avg2_04,psd_avg3_04,
# psd_avg4_04,psd_avg5_04,psd_avg6_04,psd_avg7_04,psd_avg8_04),axis=None)
psd_avg_all = np.concatenate((psd_avg_01, psd_avg_02, psd_avg_03), axis=None)

dll_01 = np.concatenate(
    (Dll1_01, Dll2_01, Dll3_01, Dll4_01, Dll5_01, Dll6_01, Dll7_01, Dll8_01),
    axis=None)
dll_02 = np.concatenate(
    (Dll1_02, Dll2_02, Dll3_02, Dll4_02, Dll5_02, Dll6_02, Dll7_02, Dll8_02),
    axis=None)
dll_03 = np.concatenate(
    (Dll1_03, Dll2_03, Dll3_03, Dll4_03, Dll5_03, Dll6_03, Dll7_03, Dll8_03),
    axis=None)
# dll_04 = np.concatenate((Dll1_04,Dll2_04,Dll3_04,Dll4_04,Dll5_04,Dll6_04,
# Dll7_04,Dll8_04),axis=None)
Dll_all = np.concatenate((dll_01, dll_02, dll_03), axis=None)

'''
Grab time from L2 files
'''

# ---- TODO : 1) make this part of an earlier step when we read in the l2
#  file the first time
#             2) make concatenate able to take in unlimited dimension
# note: this is only here because I was only reading in the pickle file and
# needed to re-read the L2 to get time


# t1 = nc.Dataset('/Users/aspen.davis/Documents/radiation_belt/20220809
# /dn_magn-l2-hires_g17_d20220808_v1-0-1.nc')
# t1 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/dn_magn-l2-hires_g18_d20230227_v1-0-1.nc')
t1 = nc.Dataset(
    '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event/dn_magn-l2'
    '-hires_g16_d20230227_v1-0-1.nc')
time_ob_native1 = t1['time']
time_ob_1 = np.asarray(time_ob_native1).flatten()

# t2 = nc.Dataset('/Users/aspen.davis/Documents/radiation_belt/20220809
# /dn_magn-l2-hires_g17_d20220809_v1-0-1.nc')
# t2 = nc.Dataset('/Users/aspen.davis/Documents/radiation_belt/data/G16L2
# /dn_magn-l2-hires_g16_d20220402_v1-0-1.nc')
# t2 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/dn_magn-l2-hires_g18_d20230228_v1-0-1.nc')
t2 = nc.Dataset(
    '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event/dn_magn-l2'
    '-hires_g16_d20230228_v1-0-1.nc')
time_ob_native2 = t2['time']
time_ob_2 = np.asarray(time_ob_native2).flatten()

# t3 = nc.Dataset('/Users/aspen.davis/Documents/radiation_belt/20220809
# /dn_magn-l2-hires_g17_d20220810_v1-0-1.nc')
# t3 = nc.Dataset('/Users/aspen.davis/Documents/radiation_belt/data/G16L2
# /dn_magn-l2-hires_g16_d20220403_v1-0-1.nc')
# t3 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/dn_magn-l2-hires_g18_d20230301_v1-0-1.nc')
t3 = nc.Dataset(
    '/Users/aspen.davis/Documents/GOES/G16/G16l1b/20230226_event/dn_magn-l2'
    '-hires_g16_d20230301_v1-0-1.nc')
time_ob_native3 = t3['time']
time_ob_3 = np.asarray(time_ob_native3).flatten()

# t4 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G18/G18_data/L1b_daily
# /20230226_event/dn_magn-l2-hires_g18_d20230301_v1-0-1.nc')
# t4 = nc.Dataset('/Users/aspen.davis/Documents/GOES/G16/G16l1b
# /20230226_event/dn_magn-l2-hires_g16_d20230301_v1-0-1.nc')
# time_ob_native4 = t4['time']
# time_ob_4 = np.asarray(time_ob_native4).flatten()

time_g16 = np.concatenate((time_ob_1, time_ob_2, time_ob_3), axis=None)

'''
convert the time from seconds since J2000 epoch to a julian date time
'''

dt_g16 = []

for i, time_new in enumerate(time_g16.data):
    dt_g16.append(time_convert(time_new))

'''
create time "chunks" or time lengths to calculate tau for 
'''

# ---- TODO : make this better : accept length of time from user, not hard
#  coded

chunk_9 = [864000, 972000]
chunk_10 = [972000, 1080000]
chunk_11 = [1080000, 1188000]
chunk_12 = [1188000, 1296000]
chunk_13 = [1296000, 1404000]
chunk_14 = [1404000, 1512000]
chunk_15 = [1512000, 1620000]
chunk_16 = [1620000, 1728000]
chunk_17 = [1728000, 1836000]
chunk_18 = [1836000, 1944000]
chunk_19 = [1944000, 2052000]
chunk_20 = [2052000, 2160000]
chunk_21 = [2160000, 2268000]
chunk_22 = [2268000, 2376000]
chunk_23 = [2376000, 2484000]
chunk_24 = [2484000, 2592000]
# chunk_25 = [2592000,2700000]
# chunk_26 = [2700000,2808000]
# chunk_27 = [2808000,2916000]
# chunk_28 = [2916000,3024000]
# chunk_29 = [3024000,3132000]
# chunk_30 = [3132000,3240000]
# chunk_31 = [3240000,3348000]
# chunk_32 = [3348000,3456000]

# this re-fills the tau, again, need to make better and more flexible


long_tau = np.zeros(len(dt_g16))
long_tau[chunk_1[0]:chunk_1[1]] = tau_m_all[0]
long_tau[chunk_2[0]:chunk_2[1]] = tau_m_all[1]
long_tau[chunk_3[0]:chunk_3[1]] = tau_m_all[2]
long_tau[chunk_4[0]:chunk_4[1]] = tau_m_all[3]
long_tau[chunk_5[0]:chunk_5[1]] = tau_m_all[4]
long_tau[chunk_6[0]:chunk_6[1]] = tau_m_all[5]
long_tau[chunk_7[0]:chunk_7[1]] = tau_m_all[6]
long_tau[chunk_8[0]:chunk_8[1]] = tau_m_all[7]
long_tau[chunk_9[0]:chunk_9[1]] = tau_m_all[8]
long_tau[chunk_10[0]:chunk_10[1]] = tau_m_all[9]
long_tau[chunk_11[0]:chunk_11[1]] = tau_m_all[10]
long_tau[chunk_12[0]:chunk_12[1]] = tau_m_all[11]

long_tau[chunk_13[0]:chunk_13[1]] = tau_m_all[12]
long_tau[chunk_14[0]:chunk_14[1]] = tau_m_all[13]
long_tau[chunk_15[0]:chunk_15[1]] = tau_m_all[14]
long_tau[chunk_16[0]:chunk_16[1]] = tau_m_all[15]
long_tau[chunk_17[0]:chunk_17[1]] = tau_m_all[16]
long_tau[chunk_18[0]:chunk_18[1]] = tau_m_all[17]
long_tau[chunk_19[0]:chunk_19[1]] = tau_m_all[18]
long_tau[chunk_20[0]:chunk_20[1]] = tau_m_all[19]
long_tau[chunk_21[0]:chunk_21[1]] = tau_m_all[20]
long_tau[chunk_22[0]:chunk_22[1]] = tau_m_all[21]
long_tau[chunk_23[0]:chunk_23[1]] = tau_m_all[22]
long_tau[chunk_24[0]:chunk_24[1]] = tau_m_all[23]

# long_tau[chunk_25[0]:chunk_25[1]] = tau_m_all[24]
# long_tau[chunk_26[0]:chunk_26[1]] = tau_m_all[25]
# long_tau[chunk_27[0]:chunk_27[1]] = tau_m_all[26]
# long_tau[chunk_28[0]:chunk_28[1]] = tau_m_all[27]
# long_tau[chunk_29[0]:chunk_29[1]] = tau_m_all[28]
# long_tau[chunk_30[0]:chunk_30[1]] = tau_m_all[29]
# long_tau[chunk_31[0]:chunk_31[1]] = tau_m_all[30]
# long_tau[chunk_32[0]:chunk_32[1]] = tau_m_all[31]
# long_tau[chunk_21[0]:chunk_21[1]] = tau_m_all[32]


# same for power spectral density


long_psd = np.zeros(len(dt_g16))
long_psd[chunk_1[0]:chunk_1[1]] = psd_avg_all[0]
long_psd[chunk_2[0]:chunk_2[1]] = psd_avg_all[1]
long_psd[chunk_3[0]:chunk_3[1]] = psd_avg_all[2]
long_psd[chunk_4[0]:chunk_4[1]] = psd_avg_all[3]
long_psd[chunk_5[0]:chunk_5[1]] = psd_avg_all[4]
long_psd[chunk_6[0]:chunk_6[1]] = psd_avg_all[5]
long_psd[chunk_7[0]:chunk_7[1]] = psd_avg_all[6]
long_psd[chunk_8[0]:chunk_8[1]] = psd_avg_all[7]
long_psd[chunk_9[0]:chunk_9[1]] = psd_avg_all[8]
long_psd[chunk_10[0]:chunk_10[1]] = psd_avg_all[9]
long_psd[chunk_11[0]:chunk_11[1]] = psd_avg_all[10]
long_psd[chunk_12[0]:chunk_12[1]] = psd_avg_all[11]

long_psd[chunk_13[0]:chunk_13[1]] = psd_avg_all[12]
long_psd[chunk_14[0]:chunk_14[1]] = psd_avg_all[13]
long_psd[chunk_15[0]:chunk_15[1]] = psd_avg_all[14]
long_psd[chunk_16[0]:chunk_16[1]] = psd_avg_all[15]
long_psd[chunk_17[0]:chunk_17[1]] = psd_avg_all[16]
long_psd[chunk_18[0]:chunk_18[1]] = psd_avg_all[17]
long_psd[chunk_19[0]:chunk_19[1]] = psd_avg_all[18]
long_psd[chunk_20[0]:chunk_20[1]] = psd_avg_all[19]
long_psd[chunk_21[0]:chunk_21[1]] = psd_avg_all[20]
long_psd[chunk_22[0]:chunk_22[1]] = psd_avg_all[21]
long_psd[chunk_23[0]:chunk_23[1]] = psd_avg_all[22]
long_psd[chunk_24[0]:chunk_24[1]] = psd_avg_all[23]

# long_psd[chunk_25[0]:chunk_25[1]] = psd_avg_all[24]
# long_psd[chunk_26[0]:chunk_26[1]] = psd_avg_all[25]
# long_psd[chunk_27[0]:chunk_27[1]] = psd_avg_all[26]
# long_psd[chunk_28[0]:chunk_28[1]] = psd_avg_all[27]
# long_psd[chunk_29[0]:chunk_29[1]] = psd_avg_all[28]
# long_psd[chunk_30[0]:chunk_30[1]] = psd_avg_all[29]
# long_psd[chunk_31[0]:chunk_31[1]] = psd_avg_all[30]
# long_psd[chunk_32[0]:chunk_32[1]] = psd_avg_all[31]


print(len(long_psd))
print(len(psd_avg_all))

NFFT = 4500  # 6000
f_hp_zT, t_hp_zT, Sxx_hp_zT = signal.spectrogram(highps_z_all, fs=10,
                                                 nperseg=int(NFFT),
                                                 noverlap=int(NFFT / 4))
sgdb_hp_zT = 10 * np.log10(Sxx_hp_zT)

'''
plotting
'''

# ---- TODO : 1) make these into functions
#             2) accept user input for what kind of plots to generate? 
#             3) auto save figures 


from matplotlib.offsetbox import AnchoredText

cmap = plt.get_cmap('rainbow')
# 12,9
fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, figsize=(18, 14))
fig.tight_layout(pad=-0.5)
ax0.plot(dt_g16, highps_z_all)
ax0.set_ylabel('B_par MFA \n [nT]', fontsize=24)
ax0.xaxis.set_ticklabels([])
ax0.tick_params(axis='y', labelsize=20)
at = AnchoredText(
    "Highpass filtered, 1 mHz", prop=dict(size=16), frameon=True,
    loc='upper right')
# at = AnchoredText(
#     "Bandpass filtered, 1-10 mHz", prop=dict(size=16), frameon=True,
#     loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax0.add_artist(at)
im = ax1.pcolormesh(t_hp_zT, f_hp_zT, sgdb_hp_zT, vmin=-20, vmax=20, cmap=cmap)
ax1.axes.get_xaxis().set_visible(False)
ax1.set_ylabel('Frequency \n [Hz]', fontsize=24)
ax1.tick_params(axis='y', labelsize=20)
ax0.margins(x=0)
ax1.set_ylim([0, 0.05])
# fig.colorbar(im, ax=ax1, label='Power Spectral Density (dB)')
ax2.plot(dt_g16, long_psd, linewidth=2)
ax2.margins(x=0)
ax2.set_ylabel('Average $P^{B}$ \n [$nT^{2}$/Hz]', fontsize=24)
ax2.xaxis.set_ticklabels([])
ax2.tick_params(axis='y', labelsize=20)
at2 = AnchoredText(
    "3 hour windows", prop=dict(size=16), frameon=True, loc='upper right')
at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax2.add_artist(at2)
ax3.semilogy(dt_g16, (long_tau / 60), linewidth=2)
ax3.set_ylabel('Tau \n [hrs]', fontsize=24)
ax3.margins(x=0)
ax3.set_xlabel('Time [UT]', fontsize=24)
ax3.tick_params(axis='both', labelsize=22)

'''
*** THIS SECTION STARTS SOME OPTIONAL THINGS THAT WE SHOULD HAVE THE ABILITY 
TODO, BUT 
WILL ONLY NEED TODO ONCE 


Minimum resonant energies
'''

# ---- TODO : 1) make this into functions : calculate resoancnce output value
#             2) maybe make it an optional part of the program to run? 

import math

# note on 11/16 : pretty sure I have this correct now. Getting what seems to
# be same values
# as Paul gets on his paper for his values

w = 0.001  # *2*np.pi#0.01#*2*np.pi#0.01#*2*np.pi #f #azmuthal wave mode
# frequency, for now using the central wave frequency we settled on above
m = 2  # azmuthal wave mode
q = 1.6021766208e-19  # electron charge, C
alpha = 90  # degrees
P = 0.35 + 0.15  # *np.sin(alpha)) #sin(90) = 1 so simplify it out
B_s = 3.12E-5  # Tesla #31200#B_e #magnetic field at Earths surface
R_e = 6371000 ** 2  # earth radius, meters
L = 6.6  # L shell
m_0 = scipy.constants.m_e  # 510998#scipy.constants.m_e#0.51099895 # #MeV
c = scipy.constants.c  # speed of light, m/s

p_one = w / m * ((q * B_s * R_e) / (6 * P * L))
p_two = np.sqrt((p_one ** 2) + (m_0 ** 2 * c ** 4))
p_three = p_one + p_two
p_four = p_three - (m_0 * c ** 2)

w_res = p_four / q
w_res_MeV = w_res * 10E-6


# ---- TODO : move this function to function library


def j2000_sec_to_datetime(unix_ms):
    dt = np.array(
        [dtm.datetime(2000, 1, 1, 12) + dtm.timedelta(seconds=unix_ms[i]) for i
         in np.arange(len(unix_ms))])
    return (dt)


# In[ ]:


'''
seiss plot : electron integral flux levels 
this plots the electron flux levels from the GOES SEISS instrument 
currenlty, this loops through a folder containing seiss files, concatenates 
an array with all the 
days to run for and then plots 
'''

# --- TODO : 1) make this into small functions
#            2) think about if reading all files from a folder is the best
#            way to plot mulitple days?
#            3) make this an optional thing to run? 
#            4) include seiss integral fulx in the final output file? 

import glob

dir = '/Users/aspen.davis/Documents/radiation_belt/20230226_ext/G18_seiss/'
##dir = '/Users/aspen.davis/Documents/radiation_belt/G17/seiss/'
files = glob.glob(dir + '*.nc')
files = np.sort(files)
days = len(files)
var1 = 'AvgIntElectronFlux'
var3 = 'L2_SciData_TimeStamp'
var2 = 'time'
old = files[0]
ufiles = []
fluxes = []
times = []

# d = Dataset('/Users/aspen.davis/Documents/radiation_belt/G17/sci_mpsh-l2
# -avg5m_g17_d20220331_v1-0-3.nc')
# time1_temp=d.variables[var3][:]
# times.append(time1_temp)
# flux1_temp=d.variables[var1][:]
# fluxes.append(flux1_temp)
for i in files:
    # print(i)
    # if i.split('_')[-2] != old.split('_')[-2]:
    d = Dataset(i)
    flux1 = d.variables[var1][:]
    fluxes.append(flux1)
    time1 = d.variables[var2][:]
    times.append(time1)
    ufiles.append(i)
    old = i

nfluxes = np.asarray(fluxes)
print(np.shape(fluxes))
q = np.where(nfluxes < -1e30)
nfluxes[q] = np.nan
print(np.shape(nfluxes))
afluxes = np.reshape(nfluxes, (
288 * 16, 5))  # multiply by how many days you are doing :) (288*16,5)
times = np.asarray(times)
times = np.reshape(times, (288 * 16))  # same
atimes = j2000_sec_to_datetime(times)
# use channel 4, which is index 3


import glob

dir = '/Users/aspen.davis/Documents/radiation_belt/20230226_ext/G16_seiss/'
##dir = '/Users/aspen.davis/Documents/radiation_belt/G16/seiss/'
files16 = glob.glob(dir + '*.nc')
files16 = np.sort(files)
days16 = len(files16)
var1 = 'AvgIntElectronFlux'
var3 = 'L2_SciData_TimeStamp'
var2 = 'time'
old16 = files16[0]
ufiles16 = []
fluxes16 = []
times16 = []

# d16 = Dataset('/Users/aspen.davis/Documents/radiation_belt/G16/sci_mpsh-l2
# -avg5m_g16_d20220331_v1-0-3.nc')
# time16_temp=d16.variables[var2][:]
# times16.append(time16_temp)
# flux16_temp=d.variables[var1][:]
# fluxes16.append(flux16_temp)
for i in files16:
    # print(i)
    # if i.split('_')[-2] != old.split('_')[-2]:
    d16 = Dataset(i)
    flux16 = d16.variables[var1][:]
    fluxes16.append(flux16)
    time16 = d16.variables[var2][:]
    times16.append(time16)
    ufiles16.append(i)
    old16 = i

nfluxes16 = np.asarray(fluxes16)
q16 = np.where(nfluxes16 < -1e30)
nfluxes16[q16] = np.nan
afluxes16 = np.reshape(nfluxes16, (288 * 16, 5))
times16 = np.asarray(times16)
times16 = np.reshape(times16, (288 * 16))
atimes16 = j2000_sec_to_datetime(times16)

# ---- TODO : make plotting into functions


# %matplotlib notebook
from matplotlib.offsetbox import AnchoredText

fig, (ax0, ax1) = plt.subplots(2, figsize=(12, 6))
fig.tight_layout(pad=-0.5)
ax0.plot(atimes, afluxes[:, 3])
ax0.set_yscale('log')
ax0.set_ylim([10, 10000])
date_form = DateFormatter('%m-%d')
ax0.xaxis.set_major_formatter(date_form)
at = AnchoredText(
    "GOES-18", prop=dict(size=20), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax0.add_artist(at)
ax0.tick_params(axis='y', labelsize=20)
ax0.xaxis.set_ticklabels([])
ax0.margins(x=0)
ax1.margins(x=0)
ax0.axhline(y=1000, color='r', linestyle='--')
ax1.plot(atimes16, afluxes16[:, 3])
ax1.set_yscale('log')
ax1.set_ylim([10, 10000])
ax1.xaxis.set_major_formatter(date_form)
ax1.set_xlabel('Date [UT]', fontsize=20)
at2 = AnchoredText(
    "GOES-16", prop=dict(size=20), frameon=True, loc='upper right')
at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax1.add_artist(at2)
ax1.tick_params(axis='both', labelsize=20
                )
ax1.axhline(y=1000, color='r', linestyle='--')
ax1.set_ylabel(
    'Particle Flux \n [particles $\cdot cm^{-2} \cdot s^{-1} \cdot sr^{-1}$]',
    fontsize=14)
ax0.set_ylabel(
    'Particle Flux \n [particles $\cdot cm^{-2} \cdot s^{-1} \cdot sr^{-1}$]',
    fontsize=14)

# ---- TODO : the plot below is probably the type of comprehensive figure we
#  would want to output


# %matplotlib notebook
'''
Figures for SWWS lightning talk
'''
from matplotlib.offsetbox import AnchoredText

fig, (ax0, ax1) = plt.subplots(2, figsize=(18, 10))
fig.tight_layout(pad=0.5)
fig.autofmt_xdate()
ax0.plot(atimes, afluxes[:, 3])
ax0.set_yscale('log')
ax0.set_ylim([10, 10000])
date_form = DateFormatter('%m-%d')
ax0.xaxis.set_major_formatter(date_form)
at = AnchoredText(
    "GOES-18", prop=dict(size=24), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax0.add_artist(at)
ax0.tick_params(axis='y', labelsize=26)
ax0.xaxis.set_ticklabels([])
ax0.margins(x=0)
ax1.margins(x=0)
ax0.axhline(y=1000, color='r', linestyle='--', linewidth=2)
ax1.plot(atimes16, afluxes16[:, 3])
ax1.set_yscale('log')
ax1.set_ylim([10, 10000])
ax1.xaxis.set_major_formatter(date_form)
ax1.set_xlabel('Date [UT]', fontsize=24)
at2 = AnchoredText(
    "GOES-16", prop=dict(size=24), frameon=True, loc='upper right')
at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax1.add_artist(at2)
ax1.tick_params(axis='both', labelsize=26
                )
ax1.axhline(y=1000, color='r', linestyle='--', linewidth=2)
ax1.set_ylabel(
    'Particle Flux \n [particles $\cdot cm^{-2} \cdot s^{-1} \cdot sr^{-1}$]',
    fontsize=20)
ax0.set_ylabel(
    'Particle Flux \n [particles $\cdot cm^{-2} \cdot s^{-1} \cdot sr^{-1}$]',
    fontsize=20)

'''
Save outputs to a file
'''

# ---- TODO : 1) decide if a pickle file is if we want to output, or an h5 /
#  netCDF / etc.
#             2) save following parameters to a file : highpass_z_all,
#             b_mfa_nT_*, psd_avg_all, time,tau, psd,
#                dt_*, t_hp, f_hp, sgdb, Dll, seis values, resontant energies
#             3) make this into a function 

name = '/Users/aspen.davis/Documents/radiation_belt/20230226_ext' \
       '/g16_0227_0301_3hr_highpass.pickle'

analysis_out = {}
analysis_out['highps'] = highps_z_all
analysis_out['psd_avg'] = psd_avg_all
analysis_out['time'] = time_g16
analysis_out['tau'] = long_tau
analysis_out['tau_m'] = tau_m_all
analysis_out['psd'] = long_psd
# analysis_out['dst'] = dst_all
# analysis_out['days'] = days
analysis_out['g18_dt'] = dt_g16
analysis_out['t_hp_zT'] = t_hp_zT
analysis_out['f_hp_zT'] = f_hp_zT
analysis_out['sgdb_hp_zT'] = sgdb_hp_zT
analysis_out['Dll'] = Dll_all

f = open(name, 'wb')
pickle.dump(analysis_out, f)
f.close()
