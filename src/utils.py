from scipy import signal
import numpy as np


def butter_filter(fs, fc, N, btype):
    """
    Create butterworth filter.
    input: fs: sampling frequency, int [Hz]
           fc: cut off frequency of filter, int [Hz]
           N: filter order, int
           btype: filter type, 'highpass', 'lowpass', 'bandpass', str
    output: b, a: filter coefficients
    """
    # fs = sampling frequency [Hz], fc = cut off frequency [Hz],
    # N = filter order, btype = high / low
    w = fc / (fs/2)  # normalize the frequency
    b, a = signal.butter(N, w, btype)  # design filter
    return b, a


def apply_butter(siggy, b, a):
    """
    Apply a filter forwards and backwards
    input: siggy: signal, 1D array of ints / floats
           b: filter coefficient (see butter_filter)
           a: filter coefficient (see butter_filter)
    output: filtered signal, 1D array
    """
    # b,a = output from butter_filter, signal = 1d array
    filtered = signal.filtfilt(b, a, siggy[
        ~np.isnan(siggy)])  # apply filter forwards and back, ignore nans
    return filtered


# TO DO : find_nearest is quite fragile atm, add in conditioning
def find_nearest(array, value, tolerance=None):
    """
    Find the index of  closest element to value in the array. i.e. if you have
    [5, 5.5] and looking for 6, return 1
    input : array: 1D array of ints or floats
            value: search value to look for, int
            tolerance: maximum deviation of element in list from search value,
                        int
    output: idx: index of the closest value in array, int
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if np.abs(array - value).min() >= tolerance:
        return None
    else:
        return idx


def read_txt(file_pth):
    """
    Read a text file. NOTE: currently set up to take in 3 componenet data
    intput: file_path: path to the file, str
    output: data: data from file, nx3 array of floats
    """
    data = np.zeros((1, 3))
    with open(file_pth) as lines:
        for line in lines:
            fields = line.rstrip().split(',')
            fields = np.float32(fields)
            data = np.vstack((data, fields))
    data = np.delete(data, (0), axis=0)
    return (data)
