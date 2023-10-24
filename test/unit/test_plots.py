import unittest
import numpy as np
import sys

sys.path.insert(0, '../../src')  # noqa
from plots import plot_data


class TestPlots(unittest.TestCase):

    def test_plot_data(self):
        # Test random input data
        t_hp_z = np.linspace(0, 10, 100)
        f_hp_z = np.linspace(0, 0.05, 100)
        sgdb_hp_z = np.random.rand(100, 100) * 40 - 20
        dt_g16 = np.linspace(0, 10, 100)
        highps_z_all = np.random.uniform(-10, 10, len(dt_g16))

        try:
            plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all)
        except Exception as e:
            self.fail(f"plot_data raised an exception: {e}")

    def test_plot_data_inputs(self):
        # Test if inputs are non-array
        with self.assertRaises(Exception):
            plot_data("invalid", 1, None, [1, 2, 3], "invalid")

    def test_plot_data_empty(self):
        # Test with empty arrays
        t_hp_z = np.array([])
        f_hp_z = np.array([])
        sgdb_hp_z = np.array([])
        dt_g16 = np.array([])
        highps_z_all = np.array([])

        with self.assertRaises(Exception):
            plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all)


if __name__ == '__main__':
    unittest.main()
