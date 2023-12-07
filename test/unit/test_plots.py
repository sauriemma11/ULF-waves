import unittest
import numpy as np
import sys

sys.path.insert(0, '../../src')  # noqa
from plots import plot_data


class TestPlots(unittest.TestCase):

    def test_plot_data(self):
        # Test random input data
        dt_g16 = np.linspace(0, 24, 864000)
        avg_psd = np.linspace(0, 1000, 864000)
        t_hp_z = np.linspace(0, 24, 864000)
        avg_tau = np.linspace(0, 1000, 24)
        window_start_time = np.linspace(0, 24, 864000)
        highps_z_all = np.random.uniform(-10, 10, 864000)
        plot_data(t_hp_z, highps_z_all, window_start_time, avg_psd,
                  avg_tau, output_dir='test.png')

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
