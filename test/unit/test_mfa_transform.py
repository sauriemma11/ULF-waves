import sys
import numpy as np
import random
import os
import unittest
from numpy.random import randint

sys.path.insert(0, '../../src')  # noqa # must run from unit test directory
import utils as u  # noqa
import src.mfa_transform as mt


class TestMathLib(unittest.TestCase):
    def setUP(self):
        # create an nx3 random dataset for testing
        self.test_file_name = 'setup_test_file.txt'
        f = open(self.test_file_name, 'w')

        for i in range(100):
            rand_int = random.randint(1, 100)
            f.write(str(rand_int) + str(rand_int) + str(rand_int) + '\n')
        f.close()

    def test_get_bav(self):
        test_dat = np.random.randint(low=1, high=10000, size=(100, 3))
        bav_out = mt.get_bav(test_dat)
        self.assertIsNotNone()

        #  def test_getcol_pos(self):
        ##      query_col = randint(0, 2)  # test file has 3 columns
        #     query_val = randint(0, 100)
        #     result_col = randint(0, 2)

        #      f = open('testing_data_file.txt', 'w')
        #     for i in range(100):
        #         rand_int = random.randint(1, 100)
        #         f.write(str(rand_int) + ',' + str(rand_int) + ',' + str(
        #         rand_int)
        #                 + '\n')
        #     f.close()
        #     val = u.get_column('testing_data_file.txt', query_col, query_val,
        result_col)
        #     self.assertIsNotNone(val)

        import sys
        import numpy as np
        import random
        import os
        import unittest
        from numpy.random import randint
        sys.path.insert(0,
                        '../../src')  # noqa # must run from unit test
        # directory
        import mfa_transform as t  # noqa
        #  NOTE TO SELF : i think that this will go from the unit test
        #  directory back
        #  one to the home directory where the utility functions are


class TestMathLib(unittest.TestCase):

    def test_mfa_basic(self):
        # TO DO: make this part of a set up maybe?
        b_epn_test_set = np.zeros((1, 3))
        with open('../data/epn_test_set.txt', 'r') as f:
            for line in f:  # f:
                line_data = line.split(',')
                line_data[0], line_data[1], line_data[2] = float(
                    line_data[0]), float(line_data[1]), float(line_data[2])
                b_epn_test_set = np.vstack((b_epn_test_set, line_data))
        b_epn_test_set = np.delete(b_epn_test_set, (0), axis=0)

        out = t.get_bav(b_epn_test_set)
        self.assertIsNotNone(out)


if __name__ == '__main__':
    unittest.main()
== == == =
# use epn test dataset to test the mfa transformation
