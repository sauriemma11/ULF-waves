import random
import unittest
import numpy as np
import random
import datetime
import sys
sys.path.insert(0, "../../src/")  # noqa
import call_tau


class TestCallTau(unittest.TestCase):
    def test_nondivisible_timespan(self):
        # Will a timespan that does not go evenly into 24 return a ValueError?
        t = 0
        rand_datetime_list = []
        while t < 25:
            rand_datetime_list.append(random.random() *
                                      datetime.timedelta(days=1))
            t += 1
        ts = 15
        self.assertRaises(ValueError, call_tau.define_window,
                          rand_datetime_list, ts)
