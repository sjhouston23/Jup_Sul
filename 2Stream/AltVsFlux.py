#!/usr/bin/env python

# This code will read in the ion productions data file and plot the last run
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec



ElectronFlux="output/jelect/SE/ALTvsFLUX1000keV_test.dat"
ElectronFlux1="output/jelect/SE/ALTvsFLUX1000keV_test1.dat"
ElectronFlux2="output/jelect/SE/ALTvsFLUX1000keV_test2.dat"
ElectronFlux3="output/jelect/SE/ALTvsFLUX1000keV_test3.dat"
print ElectronFlux

altitude, Energy1, Energy2, Energy3, Energy4, Energy5 = np.loadtxt(ElectronFlux, unpack=True, skiprows=0, usecols=[0,101,161,201,231,241])
altitude1, Energy11, Energy21, Energy31, Energy41, Energy51 = np.loadtxt(ElectronFlux1, unpack=True, skiprows=0, usecols=[0,101,161,201,231,241])
altitude2, Energy12, Energy22, Energy32, Energy42, Energy52 = np.loadtxt(ElectronFlux2, unpack=True, skiprows=0, usecols=[0,101,161,201,231,241])
altitude3, Energy13, Energy23, Energy33, Energy43, Energy53 = np.loadtxt(ElectronFlux3, unpack=True, skiprows=0, usecols=[0,101,161,201,231,241])

fig1 = plt.figure()
gs = gridspec.GridSpec(2,2)

ax0 = plt.subplot(gs[0])
#plt.title('Upward Electron Fluxes', fontsize=30)
ax0.semilogy(altitude3, Energy53, label='9.75 eV', color='c', linewidth=2)
ax0.semilogy(altitude3, Energy43, label='20 eV', color='y', linewidth=2)
ax0.semilogy(altitude3, Energy33, label='50 eV', color='g', linewidth=2)
ax0.semilogy(altitude3, Energy23, label='100 eV', color='r', linewidth=2)
ax0.semilogy(altitude3, Energy13, label='975 eV', color='b', linewidth=2)
plt.grid()
font = {'size'   : 25}
mpl.rc('font', **font)
plt.text(altitude3[120], Energy53[120]+0.2, '9.75 eV')
#plt.text(altitude3[120], 1e-2, '20 eV')
#plt.text(altitude3[120], 2e-3, '50 eV')
#plt.text(altitude3[120], 2e-4, '100 eV')
plt.text(altitude3[120], 1e-8, '975 eV')
plt.annotate('20 eV', xy=(altitude3[120], Energy43[120]), xytext=(10,1e-1), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.annotate('50 eV', xy=(altitude3[120], Energy33[120]), xytext=(10,1e-3), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.annotate('100 eV', xy=(altitude3[120], Energy23[120]), xytext=(500,1e-6), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.text(2000, 1e-1, '50 keV/u')
plt.text(2750, 1e-11, 'a')
plt.ylim(10**-12,10**1)
plt.xlim(0,3000)
plt.tick_params(axis='x', which='major', length=6, width=2, labelsize=25, pad=7)
plt.tick_params(axis='x', which='minor', length=3, width=2)
plt.tick_params(axis='y', which='major', length=6, width=2, labelsize=25, pad=0)
plt.tick_params(axis='y', which='minor', length=3, width=2)
for label in ax0.yaxis.get_ticklabels()[::2]:
    label.set_visible(False)
#plt.ylabel('Upward Electron Flux [cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$]')
#plt.legend(fontsize=12)

###################################################################################################

ax1 = plt.subplot(gs[2], sharex=ax0)
ax1.semilogy(altitude2, Energy52, label='9.75 eV', color='c', linewidth=2)
ax1.semilogy(altitude2, Energy42, label='20 eV', color='y', linewidth=2)
ax1.semilogy(altitude2, Energy32, label='50 eV', color='g', linewidth=2)
ax1.semilogy(altitude2, Energy22, label='100 eV', color='r', linewidth=2)
ax1.semilogy(altitude2, Energy12, label='975 eV', color='b', linewidth=2)
plt.setp(ax0.get_xticklabels(), visible=False)
yticks = ax1.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
plt.grid()
font = {'size'  : 25}
mpl.rc('font', **font)
plt.text(altitude2[120], Energy52[120]+0.1, '9.75 eV')
#plt.text(altitude2[120], Energy42[120], '20 eV')
#plt.text(altitude2[120], Energy32[120], '50 eV')
#plt.text(altitude2[120], Energy22[120], '100 eV')
plt.text(altitude2[120], Energy12[220], '975 eV')
plt.annotate('20 eV', xy=(altitude2[120], Energy42[120]), xytext=(10,1e-1), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.annotate('50 eV', xy=(altitude2[120], Energy32[120]), xytext=(10,1e-3), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.annotate('100 eV', xy=(altitude2[120], Energy22[120]), xytext=(500,1e-5), arrowprops=dict(arrowstyle='->'))#,fontsize=30)
plt.text(2000, 1e-1, '300 keV/u')
plt.text(2750, 1e-8, 'b')
plt.xlabel('                                                                                       Altitude [km]', fontsize=30)
plt.ylabel('                            Upward Electron Flux [cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$]', fontsize=25)
plt.ylim(10**-9,10**1)
plt.xlim(0,3000)
plt.tick_params(axis='x', which='major', length=6, width=2, labelsize=25, pad=7)
plt.tick_params(axis='x', which='minor', length=3, width=2)
plt.tick_params(axis='y', which='major', length=6, width=2, labelsize=25, pad=0)
plt.tick_params(axis='y', which='minor', length=3, width=2)
for label in ax1.yaxis.get_ticklabels()[::2]:
    label.set_visible(False)
plt.subplots_adjust(hspace=.0)

###################################################################################################

ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.semilogy(altitude1, Energy51, label='9.75 eV', color='c', linewidth=2)
ax2.semilogy(altitude1, Energy41, label='20 eV', color='y', linewidth=2)
ax2.semilogy(altitude1, Energy31, label='50 eV', color='g', linewidth=2)
ax2.semilogy(altitude1, Energy21, label='100 eV', color='r', linewidth=2)
ax2.semilogy(altitude1, Energy11, label='975 eV', color='b', linewidth=2)
#plt.setp(ax1.get_xticklabels(), visible=False)
yticks = ax2.yaxis.get_major_ticks()
#yticks[7].label1.set_visible(False)
plt.grid()
font = {'size'  : 25}
mpl.rc('font', **font)
plt.text(altitude1[120], Energy51[120], '9.75 eV')
plt.text(altitude1[120], Energy41[120], '20 eV')
plt.text(altitude1[120], Energy31[120], '50 eV')
plt.text(altitude1[120], Energy21[120], '100 eV')
plt.text(altitude1[120], 3e-6, '975 eV')
plt.text(2000, 1e-1, '2 MeV/u')
plt.text(2750, 1e-5, 'c')
#plt.ylabel('                Upward Electron Flux [cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$]', fontsize=30)
plt.ylim(10**-6,10**0)
plt.xlim(0,3000)
plt.tick_params(axis='x', which='major', length=6, width=2, labelsize=25, pad=7)
plt.tick_params(axis='x', which='minor', length=3, width=2)
plt.tick_params(axis='y', which='major', length=6, width=2, labelsize=25, pad=0)
plt.tick_params(axis='y', which='minor', length=3, width=2)
plt.subplots_adjust(hspace=.0)

###################################################################################################

ax3 = plt.subplot(gs[3], sharex=ax2)
ax3.semilogy(altitude, Energy5, label='9.75 eV', color='c', linewidth=2)
ax3.semilogy(altitude, Energy4, label='20 eV', color='y', linewidth=2)
ax3.semilogy(altitude, Energy3, label='50 eV', color='g', linewidth=2)
ax3.semilogy(altitude, Energy2, label='100 eV', color='r', linewidth=2)
ax3.semilogy(altitude, Energy1, label='975 eV', color='b', linewidth=2)
plt.setp(ax2.get_xticklabels(), visible=False)
yticks = ax3.yaxis.get_major_ticks()
yticks[8].label1.set_visible(False)
plt.grid()
font = {'size'   : 25}
mpl.rc('font', **font)
plt.text(altitude[30], Energy5[30], '9.75 eV')
plt.text(altitude[30], Energy4[30], '20 eV')
plt.text(altitude[30], Energy3[30], '50 eV')
plt.text(altitude[30], Energy2[30], '100 eV')
plt.text(altitude[30], 1e-6, '975 eV')
plt.text(2000, 1e-1, '5 MeV/u')
plt.text(2750, 1e-6, 'd')
plt.ylim(10**-7,10**0)
plt.xlim(0,3000)
plt.tick_params(axis='x', which='major', length=6, width=2, labelsize=25, pad=7)
plt.tick_params(axis='x', which='minor', length=3, width=2)
plt.tick_params(axis='y', which='major', length=6, width=2, labelsize=25, pad=0)
plt.tick_params(axis='y', which='minor', length=3, width=2)
#plt.xlabel('Altitude [km]', fontsize=30)
#plt.ylabel('Upward Electron Flux [cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$]')
#plt.legend(fontsize=12)
plt.subplots_adjust(hspace=.0)
plt.show()
