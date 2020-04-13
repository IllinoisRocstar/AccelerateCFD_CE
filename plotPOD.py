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
EigVals = my_data[:,3]
pltNum = POD.size
color = 'tab:red'

fig, ax1 = plt.subplots()
ax1.set_xlabel('basis number')
ax1.set_ylabel('Cumulative energy (%)',color=color)
ax1.plot(POD[0:pltNum],Energy[0:pltNum],color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xlim([1,pltNum+1])
ax1.grid()
  
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Energy per basis (%)',color=color)
ax2.plot(POD[0:pltNum],IndEnergy[0:pltNum],color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_xlim([1,pltNum+1])
ax2.grid()
    
plt.title('Flow Energy content in POD basis')

plt.show()