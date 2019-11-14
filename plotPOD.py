################################################################################
# AccelerateCFD_CE
# Copyright 2018-2019 Illinois Rocstar LLC
# 
# This utility helps to visualize the cumulative energy content in POD modes.
# This utility is also provided under General Public License Version 3.
#
################################################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt

my_data = genfromtxt('podEnergy.csv', delimiter=',',skip_header=1)

POD = my_data[:,0]
Energy = my_data[:,2]
IndEnergy = my_data[:,1]

plt.figure(1)
plt.plot(POD,Energy, label='Basis Energy Line')
plt.title("Cummulative Energy contained in POD modes for velocity")
plt.xlabel("Number POD modes")
plt.ylabel("Energy (%)")
plt.plot(POD,np.linspace(100,100,POD.size), label='100% Energy Line')
plt.plot(POD,np.linspace(90,90,POD.size), label='90% Energy Line')
plt.legend()

plt.figure(2)
plt.plot(POD,IndEnergy, label='Basis Energy Line')
plt.title("Flow Energy contained in each POD mode")
plt.xlabel("Number POD modes")
plt.ylabel("Energy (%)")
plt.show()