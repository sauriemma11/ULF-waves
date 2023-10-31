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
#
# TODO: WHEN YOU HAVE RANDOMNESS, RUN THINGS FOR MANY ITTERATIONS
class TestMathLib(unittest.TestCase):
    # def setUP(self):
    #     # create an nx3 random dataset for testing
    #     self.test_file_name = 'setup_test_file.txt'
    #     f = open(self.test_file_name, 'w')
    #
    #     for i in range(100):
    #         rand_int = random.randint(1, 100)
    #         f.write(str(rand_int) + str(rand_int) + str(rand_int) + '\n')
    #     f.close()

    def test_get_bav(self):
        test_dat = np.random.randint(low=1, high=10000, size=(100, 3))
        empty_dat = []
        av = np.mean(test_dat, axis=0)
        self.assertIsNotNone(mt.get_bav(test_dat))
        self.assertEqual(np.shape(mt.get_bav(test_dat)), np.shape(test_dat))
        self.assertEqual(np.all(mt.get_bav(test_dat)), np.all(av))
        self.assertRaises(TypeError, mt.get_bav, empty_dat)

    def test_get_bav_sci(self):
        # test consistency by comparison to a validated dataset
        # file to feed into get_bav function
        b_epn = u.read_txt('../../test/unit/test_data/epn_test_set.txt')
        # scientifically validated file
        b_av_sci = u.read_txt('../../test/unit/test_data/b_av_test_set.txt')
        self.assertAlmostEqual(np.all(mt.get_bav(b_epn)), np.all(b_av_sci))

    def test_compute_mfa(self):
        test_bav = np.random.randint(low=1, high=10000, size=(100, 3))
        test_dat = np.random.randint(low=1, high=10000, size=(100, 3))
        zero_bav = np.zeros((100, 3))
        wrong_sz = 2
        nan_dat = np.empty((100, 3))
        nan_dat[:] = np.nan

        self.assertIsNotNone(mt.compute_mfa(test_bav, test_dat))
        self.assertEqual(np.shape(mt.compute_mfa(test_bav, test_dat)),
                         np.shape(test_dat))
        self.assertRaises(TypeError, mt.compute_mfa, wrong_sz, zero_bav)
        self.assertIsNotNone(mt.compute_mfa(nan_dat, test_dat))

    def test_cmopute_mfa_sci(self):
        # test consistency by comparison to a validated dataset
        b_epn = u.read_txt('../../test/unit/test_data/epn_test_set.txt')
        b_av_sci = u.read_txt('../../test/unit/test_data/b_av_test_set.txt')
        b_mfa_sci = u.read_txt('../../test/unit/test_data/b_mfa_test_set.txt')
        self.assertAlmostEqual(np.all(mt.compute_mfa(b_av_sci, b_epn)),
                               np.all(b_mfa_sci))
    def test_background_sub(self):
        test_in = np.random.randint(low=1, high=10000, size=(100, 3))
        wrong_sz = 2
        nan_dat = np.empty((100, 3))
        nan_dat[:] = np.nan
        self.assertIsNotNone(mt.background_sub(test_in))
        self.assertEqual(np.shape(mt.background_sub(test_in)),
                         np.shape(test_in))
        self.assertRaises(TypeError, mt.background_sub, wrong_sz)
        self.assertIsNotNone(mt.compute_mfa(nan_dat, test_in))

    def test_background_sub_sci(self):
        # test consistency by comparison to a validated dataset
        b_mfa_sci = u.read_txt('../../test/unit/test_data/b_mfa_test_set.txt')
        b_mfa_bsub_sci = u.read_txt('../../test/unit/test_data/'+
                                    'b_mfa_bsub_test_set.txt')
        self.assertAlmostEqual(np.all(mt.background_sub(b_mfa_sci)),
                               np.all(b_mfa_bsub_sci))
    def test_main(self):
        test_in = np.random.randint(low=1, high=10000, size=(1000, 3))
        wrong_sz = 2
        nan_dat = np.empty((100, 3))
        nan_dat[:] = np.nan
        self.assertIsNotNone(mt.main(test_in))
        self.assertEqual(np.shape(mt.main(test_in)),
                         np.shape(test_in))
        self.assertRaises(TypeError, mt.main, wrong_sz)
        self.assertRaises(ValueError, mt.main, nan_dat)


