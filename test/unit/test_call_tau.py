import random
import unittest
import numpy as np
import random
import datetime
import sys
sys.path.insert(0, "../../src/")  # noqa
import call_tau


class TestCallTau(unittest.TestCase):
    # Note that negative tests for error handling is done in functional tests;
    # there is not error handling in individual functions
    def test_concat_basic_func(self):
        rand_bfa = []
        i = 0
        while i < 864000:
            rand_bfa.append([random.random(), random.random(), random.random()])
            i += 1
        print(type(random.random()))
        results = call_tau.concat_tau(rand_bfa, 864000,
                                      [0.001, 0.01], 2, "high", 1)
        print(results["tau"])
        self.assertEqual(len(results["tau"]), 24)
