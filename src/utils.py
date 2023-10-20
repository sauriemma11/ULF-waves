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
