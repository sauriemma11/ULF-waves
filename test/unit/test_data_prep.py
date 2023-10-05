import unittest
import sys
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



if __name__ == '__main__':
    unittest.main()
