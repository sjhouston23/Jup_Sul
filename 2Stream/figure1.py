#!/usr/bin/env python

# This code will read in the H production data files and plot the last run
import numpy as np
import matplotlib.pyplot as plt

#energy=1 #1,1500,2 ##1 for 1MeV, 1500 for 1500keV, 2 for 2MeV
#trial=1

#energy, trial = np.loadtxt("energytrial.dat", unpack=True, skiprows=0, usecols=[0,1])
#energy=int(energy)


    #energy=int(energies[e])
    #/SecElectrons/Output/1500keV/2strforward1500keV2.dat
felectrons="figure_1.dat"
print felectrons
    # Import the production data to be plotted
ebin, km300, km350, km400, km500, km750, km1000, km1500, km2000 = np.loadtxt(felectrons, unpack=True, skiprows=0, usecols=[0,1,2,3,4,5,6,7,8])
    #Columns, c, are based off of altitude where: Alt = 3000 - 2(c-1) -> c = 1500 - (Alt/2) + 1
    #example: 2000km -> c = 1500km - (2000km/2) + 1 = 1500km - 1000km + 1 = Coulmn 501 (which yields electron energies for forward electrons produced at 2000km)

fig1 = plt.figure()
plt.title('Forward Secondary Electron Production - [2 MeV/u Oxygen Ion]')

#plt.loglog(ebin, km300, label='300km',  color='k')
#plt.loglog(ebin, km350, label='350km',  color='darkmagenta')
plt.loglog(ebin, km400, label='400km',  color='mediumblue')
plt.loglog(ebin, km500, label='500km',  color='deepskyblue')
plt.loglog(ebin, km750, label='750km',  color='seagreen')
plt.loglog(ebin, km1000, label='1000km', color='limegreen')
plt.loglog(ebin, km1500, label='1500km', color='yellow')
plt.loglog(ebin, km2000, label='2000km', color='orange')
#print km2000

plt.grid()
#plt.ylim(10e-13,10e-4)
plt.xlim(1,10000)
plt.xlabel('Electron Energy [eV]')
plt.ylabel('Electron Production Rate [cm$^{-3}$ s$^{-1}$ eV$^{-1}$]')
plt.legend(fontsize=12)
#plt.savefig("./SecElectrons/Output/"+repr(energy)+"keV/Figure1_"+repr(energy)+"keV"+repr(trial)+".eps")
plt.show()
