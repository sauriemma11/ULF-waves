'''
Library of functions used throughout code base
- Functions:
    * butter_filter -- create a Butterworth filter
    * apply_butter -- apply a filter forwards and backwards
    * find_nearest -- find the index of  closest element to value in the array
    * read_txt -- read a text file (of 3-component data)
'''
from scipy import signal
import numpy as np


def butter_filter(fs, fc, N, btype):
    """
    Create Butterworth filter

    Parameters
    ----------
    fs: int [Hz]
        Sampling frequency
    fc: int [Hz]
        Cut off frequency of filter
    N: int
        Filter order
    btype: str
        Filter type, 'highpass', 'lowpass', 'bandpass'

    Outputs
    -------
    b, a: filter coefficients
    """
    # normalize the frequency
    w = fc / (fs/2)

    # design filter
    b, a = signal.butter(N, w, btype)
    return b, a


def apply_butter(siggy, b, a):
    """
    Apply a filter forwards and backwards

    Parameters
    ----------
    siggy: 1D array of ints / floats
        Signal
    b: (see butter filter)
        Filter coefficient
    a: (see butter_filter)
        Filter coefficient

    Outputs
    -------
    Filtered: 1D array
        Signal
    """
    filtered = signal.filtfilt(b, a, siggy[
        ~np.isnan(siggy)])  # apply filter forwards and back, ignore nans
    return filtered


def find_nearest(array, value, tolerance=None):
    """
    Find the index of  closest element to value in the array.
    (i.e. if you have [5, 5.5] and looking for 6, return 1)

    Parameters
    ----------
    array: 1D array of ints or floats
    value: int
        Search value to look for
    tolerance: int
        Maximum deviation of element in list from search value

    Outputs
    -------
    idx: int
        Index of the closest value in array, int
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if tolerance is None:
        return idx
    if np.abs(array - value).min() >= tolerance:
        print('No value within tolerance range')
        return None
    else:
        return idx


def read_txt(file_pth):
    """
    Read a text file. NOTE: currently set up to take in 3-componenet data

    Parameters
    ----------
    file_path: str
        Full path to the file

    Outputs
    -------
    data: nx3 array of floats
        Data from file,
    """
    data = np.zeros((1, 3))
    with open(file_pth) as lines:
        for line in lines:
            fields = line.rstrip().split(',')
            fields = np.float32(fields)
            data = np.vstack((data, fields))
    data = np.delete(data, (0), axis=0)
    return (data)
