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

    def test_windows_basic_func(self):
        # Does define_window produce the expected number of windows?
        # Create a list of one day of time stamps that would match data input
        current_time = datetime.datetime(2018, 10, 20, 0, 0, 0)
        sample_hourly_timestamps = []
        while current_time < datetime.datetime(2018, 10, 21, 0, 0, 0):
            sample_hourly_timestamps.append(current_time)
            current_time += datetime.timedelta(microseconds=100000)
        ts = 3
        data_windows = call_tau.define_window(sample_hourly_timestamps, ts)
        self.assertEqual(len(data_windows), 24/ts)

    def test_concat_basic_func(self):
        rand_bfa = []
        rand_fband = []
        i = 0
        while i < 864000:
            rand_fband.append(random.sample(range(0, 10), 2))
            rand_bfa.append(random.random())
            i += 1
        results = call_tau.concat_tau(rand_bfa, rand_fband)
        print(results["tau"])
        self.assertEqual(len(results["tau"]), 24)
