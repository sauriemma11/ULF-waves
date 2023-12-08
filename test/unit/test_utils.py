import sys
import numpy as np
import random
import unittest
from numpy.random import randint
sys.path.insert(0, '../../src')  # noqa
import utils as u
import scipy.signal as signal


class TestMathLib(unittest.TestCase):

    def test_butter_filter(self):
        fs = random.randint(1, 10000)
        fc = random.randint(1, 1000)
        N = random.randint(1, 5)
        self.assertIsNotNone(u.butter_filter(fs, fc, N, 'highpass'))
        self.assertIsNotNone(u.butter_filter(fs, fc, N, 'lowpass'))

    def test_apply_butter(self):
        siggy = np.random.randint(low=1, high=10000, size=(100, 3))
        w = 0.005 / (10/2)
        b, a = signal.butter(1, w, 'highpass')
        self.assertIsNotNone(u.apply_butter(siggy, b, a))
        self.assertRaises(TypeError, u.apply_butter, 'hello', b, a)

    def test_find_nearest(self):
        array = np.random.randint(low=1, high=1000, size=(100))
        value_in = np.random.randint(low=1, high=1000)
        value_not = np.random.randint(low=1001, high=10000)
        self.assertIsNotNone(u.find_nearest(array, value_in, tolerance=None))
        self.assertIsNone(u.find_nearest(array, value_not, tolerance=1))

    def test_read_txt(self):
        test_file = '../../test/unit/test_data/epn_test_set.txt'
        self.assertEqual(np.shape(u.read_txt(test_file)),  (10000, 3))
        self.assertRaises(FileNotFoundError, u.read_txt, 'nofile.txt')


if __name__ == '__main__':
    unittest.main()
