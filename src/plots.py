# This is where I will be writing documentation for the plotting functionality
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.offsetbox import AnchoredText

# For this code chunk to work, need t_hp_z, f_hp_z, and sgdb_hp_z
# These come from signal.spectrogram -- which needs results from butter_filter and apply_butter (from get_tau)
fig, ax = plt.subplots(figsize=(12,6))
cmap = plt.get_cmap('rainbow')
im = plt.pcolormesh(t_hp_z, f_hp_z, sgdb_hp_z,vmin=-20,vmax=20,cmap=cmap)
formatter = matplotlib.ticker.FuncFormatter(timeTicks)
plt.gca().xaxis.set_major_formatter(formatter)
plt.ylim([0,0.05])
fig.colorbar(im, ax=ax, label='Power Spectral Density (dB)')

# For this code chunk to work, also need butter results AND dt_g16 (time data from data prep?)
cmap = plt.get_cmap('rainbow')
#12,9
fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,figsize=(18,14))
fig.tight_layout(pad=-0.5)
ax0.plot(dt_g16,highps_z_all)
ax0.set_ylabel('B_par MFA \n [nT]',fontsize=24)
ax0.xaxis.set_ticklabels([])
ax0.tick_params(axis='y',labelsize=20)
at = AnchoredText(
    "Highpass filtered, 1 mHz", prop=dict(size=16), frameon=True, loc='upper right')
# at = AnchoredText(
#     "Bandpass filtered, 1-10 mHz", prop=dict(size=16), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax0.add_artist(at)
im = ax1.pcolormesh(t_hp_zT, f_hp_zT, sgdb_hp_zT, vmin=-20, vmax=20, cmap=cmap)
ax1.axes.get_xaxis().set_visible(False)
ax1.set_ylabel('Frequency \n [Hz]',fontsize=24)
ax1.tick_params(axis='y',labelsize=20)
ax0.margins(x=0)
ax1.set_ylim([0,0.05])
#fig.colorbar(im, ax=ax1, label='Power Spectral Density (dB)')
ax2.plot(dt_g16,long_psd,linewidth=2)
ax2.margins(x=0)
ax2.set_ylabel('Average $P^{B}$ \n [$nT^{2}$/Hz]',fontsize=24)
ax2.xaxis.set_ticklabels([])
ax2.tick_params(axis='y',labelsize=20)
at2 = AnchoredText(
    "3 hour windows", prop=dict(size=16), frameon=True, loc='upper right')
at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax2.add_artist(at2)
ax3.semilogy(dt_g16,(long_tau/60),linewidth=2)
ax3.set_ylabel('Tau \n [hrs]',fontsize=24)
ax3.margins(x=0)
ax3.set_xlabel('Time [UT]',fontsize=24)
ax3.tick_params(axis='both', labelsize=22)
