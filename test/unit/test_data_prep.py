import unittest
import sys

import numpy as np

# sys.path.insert(0, '../../src')  #noqa
from src.data_prep import *
import os
import random

class TestDataPrep(unittest.TestCase):

    def test_read_nc_file_not_empty(self):
        file_name = '../../data/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc'

        variables = read_nc_file(file_name)

        # pass if returned variables are not empty
        self.assertTrue(variables)

    def test_read_nc_file_is_dict(self):
        file_name = '../../data/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc'
        variables = read_nc_file(file_name)

        # pass if the returned variables are a type dict
        self.assertIsInstance(variables, dict)

    def test_read_nc_file_missing(self):
        fake_file_path = 'doesnt_exist.nc'
        with self.assertRaises(FileNotFoundError):
            read_nc_file(fake_file_path)

    def test_time_convert(self):
        seconds_2000 = [0, 3600, 7200]  # 0 seconds, 1 hour, 2 hours
        converted_time = time_convert(seconds_2000)
        expected_time = ['2000-01-01 12:00:00', '2000-01-01 13:00:00',
                         '2000-01-01 14:00:00']

        # Convert datetime objects to strings for comparison
        converted_time_str = [time.strftime('%Y-%m-%d %H:%M:%S') for time in
                              converted_time]
        self.assertEqual(converted_time_str, expected_time)


if __name__ == '__main__':
    unittest.main()
