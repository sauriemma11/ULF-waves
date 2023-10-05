# calculate the tau value. this assumes that the call_tau 
# module is set up to repeatedly call this module for every "chunk"
# of time that we need to calculate tau for. i.e. this will do one tau cal


#TO DO : break this into smaller chunks, write tests

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
