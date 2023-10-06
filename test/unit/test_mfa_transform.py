import sys
import numpy as np
import random
import os
import unittest
from numpy.random import randint
sys.path.insert(0, '../../src')  # noqa # must run from unit test directory
import mfa_transform as t  # noqa
#  NOTE TO SELF : i think that this will go from the unit test directory back
#  one to the home directory where the utility functions are


class TestMathLib(unittest.TestCase):

    def test_mfa_basic(self):
        #TO DO: make this part of a set up maybe?
        b_epn_test_set = np.zeros((1,3))
        with open('../data/epn_test_set.txt', 'r') as f:
            for line in f: #f:
                line_data = line.split(',')
                line_data[0],line_data[1],line_data[2] = float(line_data[0]),float(line_data[1]),float(line_data[2])
                b_epn_test_set=np.vstack((b_epn_test_set,line_data))
        b_epn_test_set = np.delete(b_epn_test_set, (0), axis=0)
        
        out = t.get_bav(b_epn_test_set)
        self.assertIsNotNone(out)

if __name__ == '__main__':
    unittest.main()
