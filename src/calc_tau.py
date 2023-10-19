'''
Module to calculate the tau value
-Inputs:
    * dictionary (output from mfa_transform.py) with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
-Outputs:
    * individual returns from functions, particularly
        mx1 tau
        mx1 D_LL
        mx1 PSD
        mx1 time
        mx1 frequencies
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

# calculate the tau value. this assumes that the call_tau
# module is set up to repeatedly call this module for every "chunk"
# of time that we need to calculate tau for. i.e. this will do one tau cal


# step 1: filter chunk
# step 2: create spectrogram
# step 3: define frequency band
# step 4: avg psd
# step 5: D_LL calculation
# step 6: tau (invert D_LL)

# TO DO: GET RID OF GLOBAL VARIABLES HERE
# THINGS TO ADD TO CONFIG: frequency/Pc band
#   -highpass or bandpass for tau filtering

def filter_b(b_mfa,ftype,comp):
    # TO DO: option for which componen to get waves for?
    # TO DO: optoin for frequency band? dpdds on application..
    # could make separate config files for difference pc bands?
    # TO DO: option to highipass or banpass filter
    """
    **NOTE : ASSUMES B_MFA IS ALREADY "TIME CHUNKED"
    setp 1: highpass filter to get frequencies of interest
    input: b_mfa: magnetic field measurements in MFA coordinates, assumed to
           be one time window, size nx3
    output: highps_z: highpass filtered componenent of b_mfa. size 1xn
    parameters: fs: sampling frequency, Hz
                fc: cut off frequency

    """

    fs = 10 #sampling frequncy, Hz, usually 10 Hz for mag
    N = 1
    if ftype == 'high':
        fc = 0.001 #cut off frequency 10mHz = 0.01 Hz, 1mHz = 0.001 Hz
        b,a = butter_filter(fs ,fc,N, ftype)
        filt_z = apply_butter(b_mfa[:][:,comp],b,a)
    if ftype == 'bandpass'
        w1 = 0.001 / (fs/2) #normalize the frequency
        w2 = 0.01 / (fs/2)
        b,a = butter_filter(fs ,np.asarray([w1,w2]),N, ftype) #bandpass filter
        filt_z = apply_butter(b_mfa[:][:,comp],b,a)

    return(filt_z)

def spect(filt_z):

    NFFT = 4500

    f_hp_spectrum, t_hp_spectrum, Sxx_hp_spectrum = signal.spectrogram(filt_z,
            fs=10, nperseg=int(NFFT),noverlap = int(NFFT / 4),scaling='spectrum',
            return_onesided=True,mode='psd')
    return(f_hp_spectrum,t_hp_spectrum,Sxx_hp_spectrum)

def get_f_band(f_hp_spectrum):
    # TO DO: band option in config?

    band = [0.001,0.01] #band = frequency band limits,: 1-10mHz, use 0.011 so includes 10mHz
    delta_f = band[1]-band[0]
    #get start and stop index for frequency band
    band_st = find_nearest(f_hp_spectrum,band[0])
    band_sp = find_nearest(f_hp_spectrum,band[1])
    return(band_st,band_sp)

def avg_psd(Sxx_hp_spectrum):

    band_pwrs = []
    for t in np.arange(0,len(t_hp_z)):
        band_sum = np.sum(Sxx_hp_spectrum[band_st:band_sp,t])
        band_pwr_tmp = band_sum / delta_f
        band_pwrs.append(band_pwr_tmp)
    psd_avg = np.mean(band_pwrs)
    print(psd_avg)
    return(psd_avg)

def calc_Dll(psd_avg):
    L = 6.6
    P_b = psd_avg#psd #psd_time #average power of field pertubations over all frequencies being used
                        #using paralell componenet bc leads to E accelerating electrons in radial directions
    f = (0.001 + 0.01) / 2#(0.0012 + 0.01) / 2 #central frequency, using middle of frequency band 1-10mHz
    B_e = 31200  #B field at equator : equitaorial magnetic field strength at surface of earth
                # = 0.312 G = 31200 nT
    Dll =( (L**8 * 4 * np.pi**2 ) / (9 * (8* B_e**2)) ) * (P_b * f**2)
    print('Sandhu', Dll)
    return(psd_avg)

def calc_tau(D_ll):

    tau_s = 1/Dll
    tau_m = tau_s / 60
    print('Sandhu tau : ',tau_m )
    return(tau_m)

# ------------------------------------------ original getTau

def getTau(chunk,b_filt_nT_mfai):
    b_hp_filt_short = b_filt_nT_mfai[chunk[0]:chunk[1]]
    fs = 10 #sampling frequncy, Hz, usually 10 Hz for mag
    fc = 0.001 #cut off frequency 10mHz = 0.01 Hz, 1mHz = 0.001 Hz
    w1 = 0.001 / (fs/2) #normalize the frequency
    w2 = 0.01 / (fs/2)
    N = 1 #filter order
    #btype = 'bandpass'# 'high' #high or low pass filterbb\\
    btype = 'high'

   # b,a = butter_filter(fs ,np.asarray([w1,w2]),N, btype) #bandpass filter
    b,a = butter_filter(fs ,fc,N, btype)
    highps_z = apply_butter(b_hp_filt_short[:][:,2],b,a)

    NFFT = 4500

    f_hp_spectrum, t_hp_spectrum, Sxx_hp_spectrum = signal.spectrogram(highps_z,fs=10, nperseg=int(NFFT),noverlap = int(NFFT / 4),scaling='spectrum',return_onesided=True,mode='psd')

    f_hp_z, t_hp_z, Sxx_hp_z = signal.spectrogram(highps_z,fs=10, nperseg=int(NFFT),noverlap = int(NFFT / 4))
    sgdb_hp_z = 10 * np.log10(Sxx_hp_z)

    fig, ax = plt.subplots(figsize=(12,6))
    cmap = plt.get_cmap('rainbow')
    im = plt.pcolormesh(t_hp_z, f_hp_z, sgdb_hp_z,vmin=-20,vmax=20,cmap=cmap)
    formatter = matplotlib.ticker.FuncFormatter(timeTicks)
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.ylim([0,0.05])
    fig.colorbar(im, ax=ax, label='Power Spectral Density (dB)')

    band = [0.001,0.01] #band = frequency band limits,: 1-10mHz, use 0.011 so includes 10mHz
    delta_f = band[1]-band[0]
    #get start and stop index for frequency band
    band_st = find_nearest(f_hp_spectrum,band[0])
    band_sp = find_nearest(f_hp_spectrum,band[1])

    band_pwrs = []
    for t in np.arange(0,len(t_hp_z)):
        band_sum = np.sum(Sxx_hp_spectrum[band_st:band_sp,t])
        band_pwr_tmp = band_sum / delta_f
        band_pwrs.append(band_pwr_tmp)
    psd_avg = np.mean(band_pwrs)
    print(psd_avg)
    #new Dll, updated from how we interperet Osmane 2023
    L = 6.6
    P_b = psd_avg#psd #psd_time #average power of field pertubations over all frequencies being used
                        #using paralell componenet bc leads to E accelerating electrons in radial directions
    f = (0.001 + 0.01) / 2#(0.0012 + 0.01) / 2 #central frequency, using middle of frequency band 1-10mHz
    B_e = 31200  #B field at equator : equitaorial magnetic field strength at surface of earth
                # = 0.312 G = 31200 nT
    Dll =( (L**8 * 4 * np.pi**2 ) / (9 * (8* B_e**2)) ) * (P_b * f**2)
    print('Sandhu', Dll)

    tau_s = 1/Dll
    tau_m = tau_s / 60
    print('Sandhu tau : ',tau_m )

    return(highps_z,f_hp_z, t_hp_z,sgdb_hp_z,psd_avg,Dll,tau_m)
