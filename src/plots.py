"""
Module to create final result file and visuals
-Inputs:
    * dictionary (output from call_tau.py) with
        mx1 tau
        mx1 D_LL
        mx1 PSD
        mx1 time : dt_g16
        mx1 frequencies : t_hp_z, f_hp_z, sg_hp_z
    * user input
-Outputs:
    * saved plot files
    * csv file with key result parameters
-Used in:
    * main.py

"""

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np


def plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all):
    """
        Create a 4-panel subplot

        Args:
            t_hp_z (numpy.ndarray): Array of time values.
            f_hp_z (numpy.ndarray): Array of frequency values.
            sgdb_hp_z (numpy.ndarray): Spectrogram data.
            dt_g16 (numpy.ndarray):
            highps_z_all (numpy.ndarray):

        Example:
            plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all)
    """
    fig, ax = plt.subplots(nrows=4)

    cmap = plt.get_cmap('rainbow')

    ax[0].plot(dt_g16, highps_z_all, linewidth=1)
    ax[0].set_ylabel('B_par MFA\n[nT]')
    ax[0].xaxis.set_ticklabels([])
    ax[0].tick_params(axis='y')

    at = AnchoredText("Highpass filtered, 1 mHz", frameon=True,
                      loc='upper right')
    at.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    ax[0].add_artist(at)

    ax[1].pcolormesh(t_hp_z, f_hp_z, sgdb_hp_z, vmin=-20, vmax=20, cmap=cmap)
    ax[1].axes.get_xaxis().set_visible(False)
    ax[1].set_ylabel('Frequency\n[Hz]')
    ax[1].tick_params(axis='y')

    ax[2].plot(dt_g16, np.random.rand(len(dt_g16)) * 100,
               linewidth=1)  # Replace with your data for long_psd
    ax[2].set_ylabel('Average $P^{B}$\n[$nT^2$/Hz]')
    ax[2].xaxis.set_ticklabels([])

    at2 = AnchoredText("3 hour windows", frameon=True, loc='upper right')
    at2.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    ax[2].add_artist(at2)

    ax[3].semilogy(dt_g16, np.random.rand(len(dt_g16)) * 100,
                   linewidth=1)  # Replace with your data for long_tau
    ax[3].set_ylabel('Tau\n[hrs]')
    ax[3].set_xlabel('Time [UT]')
    ax[3].tick_params(axis='both')

    plt.show()


# Make sample data for plotting:
t_hp_z = np.linspace(0, 10, 100)
f_hp_z = np.linspace(0, 0.05, 100)
sgdb_hp_z = np.random.rand(100, 100) * 40 - 20  # 2d array
dt_g16 = np.linspace(0, 10, 100)
highps_z_all = np.random.uniform(-10, 10, len(dt_g16))

plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all)

# TODO: have option to save plot?
# TODO: save csv file
