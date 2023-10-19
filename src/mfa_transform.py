# MFA transformation steps
# to be called after data_prep
# output will go into calc_tau

# step 1 : find average background field 

# TO DO: USE CONFIG FILE TO SET UP THE FILTER: i.e. set up fs/fc/N/btype
def get_bav(b_in):
    """
    step 1: get average background field. This is done by taking a 30 min
            lowpass butterworth filter of the entire dataset. This calls 
            butter_filter and apply butter. 
    input : b_epn: magnetic field msmts in EPN coordinates, list of int size
            nx3 where n = times and ea. col is a componenet. data should be 
            accounted for fill values already.
            b_time: timestamps corresponding to b_epn, converted to dt,
            list of int size n
    output: b_av : average background field. list of int size 1x3

    other parameters : fs: sampling frequency, Hz (10 Hz default)
                    fc: cut off frequency (30 minutes)
                    N : filter order
                    btype: filter type, high/low/band pass (lowpass default)
    
    """


    fs = 10
    fc = 1/(30*60)
    N = 2
    btype = 'lowpass'
    b,a = butter_filter(fs ,fc,N, btype)
    b_avx = apply_butter(b_in[:][:,0],b,a)
    b_avy = apply_butter(b_in[:][:,1],b,a)
    b_avz = apply_butter(b_in[:][:,2],b,a)

    b_av = np.column_stack((b_avx,b_avy,b_avz))
    return(b_av)

# setp 2 : compute mfa 
def compute_mfa(b_av,b_epn):
    """
    step 2: compute MFA coordinates. See wave_analysis_steps.pdf for all dtails
    input: b_av: average background field from get_bav
           b_epn: magnetic field observations in EPN. List of int size nx3 w/
                  data prep already accoutned for 
    output : magnetic field measurements in MFA coordinates, list of int size
            nx3.
    """


    b_mfa = np.zeros((len(b_epn),3))
    for n in np.arange(len(b_epn)):
        e_par2 = b_av[n]/np.linalg.norm(b_av[n])
        e_phi2 = e_par2 * [-1,0,0]  
        e_phi3 = e_phi2 / np.linalg.norm(e_phi2)
        e_r = e_phi2* e_par2
        e_r3 = e_r / np.linalg.norm(e_r)

        mfa_trans = np.row_stack((e_r,e_phi2,e_par2)) 
        b_mfa_temp = np.matmul(b_epn[n],mfa_trans.T)
        b_mfa[n] = b_mfa_temp

    return(b_mfa)

# step 3 : background subtraction
def background_sub(b_in):
    """
    step 3: subtrack the background signal so that it is easier to see the 
            waves in the frequency band of interest
    input: b_in: magnetic field measurements in MFA coordinates, list of
           ints size nx3
    output: b_mfa_bsub: magnetif field measurements in MFA coordinates with 
            the average background field subtracted, list of int size nx3
    parameters : fs: sampling frequency, Hz (10 Hz default)
                 fc: cutoff frequency, Hz (30 min)
                 btype: filter type (highpass)
    """

    fs = 10
    fc = 1/(30*60)
    N = 2
    btype = 'high'
    b,a = butter_filter(fs ,fc,N, btype)
    b_x = apply_butter(b_in[:][:,0],b,a)
    b_y = apply_butter(b_in[:][:,1],b,a)
    b_z = apply_butter(b_in[:][:,2],b,a)

    b_mfa_bsub = np.column_stack((b_x,b_y,b_z))
    return(b_mfa_bsub)

# TO DO: snake rule to link these variables to the other modules? 
def main():
    b_avg = get_bav(b_epn)

    b_mfa = compute_mfa(b_avg,b_epn)

    b_mga_bsub = background_sub(b_mfa)


if __name__ == '__main__':
    main()
