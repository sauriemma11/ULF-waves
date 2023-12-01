import unittest
import numpy as np
import sys

sys.path.insert(0, '../../src')  # noqa
from plots import plot_data


class TestPlots(unittest.TestCase):

    def test_plot_data(self):
        # Test random input data
        windows_start_time = np.linspace(0, 864000, 24)
        time_by_data_entry = np.linspace(0, 10, 100)
        mag_field_data = np.random.uniform(-10, 10, 100)
        avg_psd = np.linspace(0, 300, len(windows_start_time))
        avg_tau = np.linspace(0, 10, len(windows_start_time))
        frequencies_2D = np.random.uniform(-20, 20, (len(mag_field_data), len(time_by_data_entry)))

        try:
            plot_data(time_by_data_entry, mag_field_data, frequencies_2D, windows_start_time, avg_psd, avg_tau)
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
