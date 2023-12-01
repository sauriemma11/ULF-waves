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
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import numpy as np


def plot_data(time_by_data_entry, mag_field_data, frequencies, windows_start_time, avg_psd, avg_tau):
    """
    time_by_data_entry, Magnetic_field_data, frequencies, windows_start_times, average_PSD, average_tau

    Create a 4-panel subplot
    """
    fig, ax = plt.subplots(nrows=4)

    cmap = plt.get_cmap('rainbow')
    print("ax1")
    ax[0].plot(time_by_data_entry, mag_field_data, linewidth=1)
    ax[0].set_ylabel('Comp. MFA\n[nT]')
    ax[0].xaxis.set_ticklabels([])
    ax[0].tick_params(axis='y')

    at = AnchoredText("Highpass filtered, 1 mHz", frameon=True,
                      loc='upper right')
    at.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    ax[0].add_artist(at)

    ax[1].pcolormesh(time_by_data_entry, mag_field_data, frequencies, vmin=-20, vmax=20, cmap=cmap)
    ax[1].axes.get_xaxis().set_visible(False)
    ax[1].set_ylabel('Frequency\n[Hz]')
    ax[1].tick_params(axis='y')

    ax[2].plot(windows_start_time, avg_psd, linewidth=1)
    ax[2].set_ylabel('Average $P^{B}$\n[$nT^2$/Hz]')
    ax[2].xaxis.set_ticklabels([])
    window_len_inhrs = 24 / len(avg_psd)
    at2 = AnchoredText(f'{window_len_inhrs} hour windows', frameon=True, loc='upper right')
    at2.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    ax[2].add_artist(at2)

    ax[3].semilogy(windows_start_time, avg_tau,
                   linewidth=1)
    ax[3].set_ylabel('Tau\n[hrs]')
    ax[3].set_xlabel('Time of day [hrs]')
    ax[3].tick_params(axis='both')
    fig.savefig('../docs/output_plot.png')

# **** Moved this all to unit test
# # Make sample data for plotting:
# t_hp_z = np.linspace(0, 10, 100)
# f_hp_z = np.linspace(0, 0.05, 100)
# sgdb_hp_z = np.random.rand(100, 100) * 40 - 20 #2d array
# dt_g16 = np.linspace(0, 10, 100)
# highps_z_all = np.random.uniform(-10, 10, len(dt_g16))
#
# plot_data(t_hp_z, f_hp_z, sgdb_hp_z, dt_g16, highps_z_all)


# TODO: have option to save plot -> in /docs
# TODO: save csv file
