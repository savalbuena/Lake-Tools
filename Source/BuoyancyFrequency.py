# Calculates buoyancy frequency.
#
# @description Calculate the buoyancy frequency (Brunt-Vaisala frequency) for a temperature
# profile.
#
#
# Param rho a numeric vector of water density previously obtained with wdensity() function in kg/m3
# Param depths a numeric vector corresponding to the depths (in m) of the rho
# measurements
# Return a vector of buoyancy frequency in units sec^-2.
# Return value has attribute "depths" which define buoyancy frequency depths (which differ from supplied depths).
# Examples
#
# A vector of water temperatures
# 	rho = wtr = c(22.51, 22.42, 22.4, 22.4, 22.4, 22.36, 22.3, 22.21, 22.11, 21.23, 16.42,
# 		15.15, 14.24, 13.35, 10.94, 10.43, 10.36, 9.94, 9.45, 9.1, 8.91, 8.58, 8.43)
#
# 	#A vector defining the depths
# 	depths = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
# 		17, 18, 19, 20)
#
# 	b.f = buoyancy.freq(wtr, depths)
#
# 	plot(b.f, attr(b.f, 'depths'), type='b',
# 		ylab='Depth', xlab='Buoyancy Frequency', ylim=c(max(depths), min(depths)))
#----Author: Jordan S Read 2011 ----
# Editied by:
# 1/30/2010 - Luke Winslow <lawinslow@gmail.com>
# Adapted for python and edited by:
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from WaterDensity import *

# Declaration of mathematical constants
g = 9.806           #[m/s2] Gravitational accelaration

def BuoyancyFrequency(rho,z):
    numdepths = len(z)
    N2 = np.zeros(numdepths-1)
    N2z = np.zeros(numdepths -1)
    N2[:] = np.nan
    N2z[:] = np.nan

    # Calculate the first derivative
    for i in range(0,numdepths-1):
        N2[i] = g/rho[i]*(rho[i+1] - rho[i])/(z[i+1] - z[i])
        N2z[i] = (z[i+1] + z[i])/2

    return N2, N2z
