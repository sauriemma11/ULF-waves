import sys
import numpy as np
import random
import os
import unittest
from numpy.random import randint
sys.path.insert(0, '../../src')  # noqa # must run from unit test directory
import utils as u  # noqa
# import src.calc_tau as tau
import calc_tau as tau
sys.path.insert(0, './test_data')


class TestMathLib(unittest.TestCase):

    def test_filter_b(self):
        test_dat = np.random.randint(low=1, high=10000, size=(100, 3))
        out_dat = np.random.randint(low=1, high=10000, size=(100))
        empty_dat = []
        nan_dat = np.copy(test_dat)
        nan_dat[20:30][:, 2] = np.asarray(np.nan).astype(int)

        self.assertIsNotNone(tau.filter_b(test_dat))
        self.assertIsNotNone(tau.filter_b(nan_dat))
        self.assertEqual(np.shape(tau.filter_b(nan_dat)), np.shape(out_dat))
        self.assertRaises(TypeError, tau.filter_b, empty_dat)
        self.assertRaises(ZeroDivisionError, tau.filter_b, test_dat, fs=0)

    def test_spect(self):
        test_dat = np.random.randint(low=1, high=10000, size=(10000))
        empty_dat = []
        nan_dat = np.copy(test_dat)
        nan_dat[20:30] = np.asarray(np.nan).astype(int)

        # to test the shapes of the outputs, hand calculate FFT size
        NFFT = 4500  # matches constant defined in tau.spect
        noverlap = NFFT / 4  # matches default in tau.spect
        f_sz = (NFFT/2)+1
        t_sz = np.round((len(test_dat)/(NFFT-noverlap)) - 1)
        # Sxx_sz = np.empty((f_sz, t_sz))

        f_out, t_out, Sxx_out = tau.spect(nan_dat)

        self.assertIsNotNone(tau.spect(test_dat))
        self.assertIsNotNone(tau.spect(nan_dat))
        self.assertEqual(len(f_out), f_sz)
        self.assertEqual(len(t_out), t_sz)
        # note: if the above two test pass, Sxx_out has correct shape
        self.assertRaises(AttributeError, tau.spect, empty_dat)
        self.assertRaises(ZeroDivisionError, tau.spect, test_dat, fs=0)

    def test_get_fband_ind(self):
        NFFT = 4500  # matches constant defined in tau.spect
        noverlap = NFFT / 4  # matches default in tau.spect
        f_sz = (NFFT / 2) + 1

        f = np.random.randint(low=0, high=100, size=int(f_sz))
        nan_dat = np.copy(f)
        nan_dat[:] = np.asarray(np.nan).astype(int)
        fband = np.random.randint(low=0, high=1, size=2)
        fband_bad = np.random.randint(low=0, high=1, size=1)
        df_check = fband[1] - fband[0]
        fband_outrange = np.random.randint(low=1000, high=10000, size=2)
        fband_outrange.sort()
        df_outrange = float(fband_outrange[1] - fband_outrange[0])
        bd1, bd2, df = tau.get_fband_ind(f, fband)

        self.assertIsNotNone(tau.get_fband_ind(f, fband))
        self.assertEqual(df_check, df)
        self.assertRaises(IndexError, tau.get_fband_ind, f, fband_bad)
        self.assertEqual(tau.get_fband_ind(f, fband_outrange),
                         (None, None, df_outrange))

    def test_avg_psd(self):
        # lots of set up: prior outputs are inputs here
        test_dat = np.random.randint(low=1, high=10000, size=(10000))
        NFFT = 4500  # matches constant defined in tau.spect
        noverlap = NFFT / 4  # matches default in tau.spect
        f_sz = (NFFT / 2) + 1
        t_sz = np.round((len(test_dat) / (NFFT - noverlap)) - 1)
        Sxx = np.random.randint(low=1, high=100, size=((int(f_sz), int(t_sz))))
        Sxx_nans = np.copy(Sxx)
        Sxx_nans[20:30] = np.asarray(np.nan).astype(int)

        band_st = int(np.random.randint(low=0, high=10, size=1))
        band_sp = int(band_st+10)
        bad_band_st = band_st*(-1)

        delta_f = band_st - band_sp

        t_xx = np.random.randint(low=1, high=10000, size=int(t_sz))

        self.assertIsNotNone(tau.avg_psd(Sxx, band_st, band_sp, delta_f, t_xx))
        self.assertIsNotNone(tau.avg_psd(Sxx_nans, band_st, band_sp, delta_f,
                                         t_xx))
        self.assertRaises(ZeroDivisionError, tau.filter_b, Sxx_nans, band_st,
                          band_sp, 0, t_xx)
        self.assertRaises(AttributeError, tau.filter_b, Sxx_nans, bad_band_st,
                          band_sp, delta_f, t_xx)
        self.assertRaises(AttributeError, tau.filter_b, Sxx_nans, band_st,
                          band_st, delta_f, t_xx)

    def test_calc_DlL(self):
        test_avg = np.random.randint(low=1, high=10000, size=1)
        test_pb = 10000
        self.assertIsNotNone(tau.calc_Dll(test_avg))
        # hand calculated expected value with avg_psd = 10000 and 0.01
        self.assertAlmostEqual(tau.calc_Dll(10000), 0.0006134722)
        self.assertAlmostEqual(tau.calc_Dll(0.01), 6.13472240224591e-10)
        self.assertRaises(TypeError, tau.calc_Dll, [])
        # note: assertAlmostEqual looks for equality w/in 7 places, which is
        # fine for this application.

    def test_calc_tau(self):
        test_dll = np.random.randint(low=1, high=10000, size=1)

        self.assertIsNotNone(tau.calc_tau(test_dll))
        # hand calculated
        self.assertAlmostEqual(tau.calc_tau(10000), 1.6667e-06)
        self.assertAlmostEqual(tau.calc_tau(0.01), 1.66666667)
        self.assertRaises(TypeError, tau.calc_tau, [])
        self.assertRaises(ZeroDivisionError, tau.calc_tau, 0)

    def test_get_tau(self):
        test_dat = np.random.randint(low=1, high=10000, size=(10000, 3))

        self.assertIsNotNone(tau.get_tau(test_dat))
        self.assertIn('Sxx', tau.get_tau(test_dat))
        self.assertIn('tau', tau.get_tau(test_dat))
        self.assertIn('D_LL', tau.get_tau(test_dat))
        self.assertIn('freqs', tau.get_tau(test_dat))
        self.assertIn('time', tau.get_tau(test_dat))
        self.assertIn('psd', tau.get_tau(test_dat))
        self.assertIn('b_filt', tau.get_tau(test_dat))
        self.assertRaises(TypeError, tau.get_tau, [])
        self.assertRaises(IndexError, tau.get_tau, test_dat, comp=4)
        self.assertRaises(IndexError, tau.get_tau, test_dat, fband=[])


if __name__ == '__main__':
    unittest.main()
