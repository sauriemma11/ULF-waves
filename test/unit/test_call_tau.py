import random
import unittest
import numpy as np
import random
import sys
sys.path.insert(0, "../../src/")  # noqa
import call_tau
import utils as u  # noqa


class TestCallTau(unittest.TestCase):
    # Note that negative tests for error handling is done in functional tests;
    # there is not error handling in individual functions
    def test_concat_basic_func(self):
        rand_bfa = []
        i = 0
        while i < 864000:
            one_bfa_entry = []
            for j in range(3):
                one_bfa_entry.append(random.random())
            rand_bfa.append(one_bfa_entry)
            i += 1
        results = call_tau.concat_tau(np.array(rand_bfa), 864000,
                                      [0.001, 0.01], 2, "high", 1)
        self.assertEqual(len(results["tau"]), 24)
        self.assertEqual(len(results["D_LL"]), 24)
        self.assertEqual(len(results["psd"]), 24)
        self.assertEqual(len(results["b_filt"]), 24)

    def test_concat_bad_timespan(self):
        rand_bfa = []
        i = 0
        while i < 864000:
            one_bfa_entry = []
            for j in range(3):
                one_bfa_entry.append(random.random())
            rand_bfa.append(one_bfa_entry)
            i += 1
        with self.assertRaises(ValueError) as context:
            results = call_tau.concat_tau(np.array(rand_bfa), 864000,
                                          [0.001, 0.01], 2, "high", 9)
        self.assertEqual(str(context.exception),
                         "Timespan must be a multiple of 24.")
