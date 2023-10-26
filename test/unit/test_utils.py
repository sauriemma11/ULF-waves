import sys
import numpy as np
import random
import os
import unittest
from numpy.random import randint
sys.path.insert(0, '../../src')  # noqa # must run from unit test directory
import utils as u  # noqa
import src.mfa_transform as mt
sys.path.insert(0, './test_data')
import mfa_transform as mt
import utils as u  # noqa
# import src.mfa_transform as mt

# TO DO: WHEN YOU HAVE RANDOMNESS, RUN THINGS FOR MANY ITTERATIONS
class TestMathLib(unittest.TestCase):

    def test_butter_filter(self):
        fs = random.randint(1, 10000)
        fc = random.randint(1, 1000)
        N = random.randint(1, 5)
        self.assertIsNotNone(u.butter_filter(fs, fc, N, 'highpass'))
        self.assertIsNotNone(u.butter_filter(fs, fc, N, 'lowpass'))
