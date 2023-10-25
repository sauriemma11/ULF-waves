from scipy import signal
import numpy as np


def butter_filter(fs, fc, N, btype):
    # fs = sampling frequency [Hz], fc = cut off frequency [Hz],
    # N = filter order, btype = high / low
    w = fc / (fs/2)  # normalize the frequency
    b, a = signal.butter(N, w, btype)  # design filter
    return b, a


def apply_butter(siggy, b, a):
    # b,a = output from butter_filter, signal = 1d array
    filtered = signal.filtfilt(b, a, siggy[
        ~np.isnan(siggy)])  # apply filter forwards and back, ignore nans
    return filtered


# TO DO : find_nearest is quite fragile atm, add in conditioning
# TO DO : ADD IN CONDITION FOR CHECKING FBAND[1] > FBAND[0]
# TO DO : add in condition for if something is really not in the set
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
  


def read_txt(file_pth):
    data = np.zeros((1, 3))
    with open(file_pth) as lines:
        for line in lines:
            fields = line.rstrip().split(',')
            fields = np.float32(fields)
            data = np.vstack((data, fields))
    data = np.delete(data, (0), axis=0)
    return (data)
