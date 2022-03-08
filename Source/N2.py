# Calculates buoyancy frequency.
#
# Description Calculate the buoyancy frequency (Brunt-Vaisala frequency) for a density
# profile.
#
# Param rho a numeric vector of water density previously obtained with wdensity() function in kg/m3
# Param depths a numeric vector corresponding to the depths (in m) of the rho
# measurements
# Return a vector of buoyancy frequency in units sec^-2.
# Return value has attribute "depths" which define buoyancy frequency depths (which differ from supplied depths).

# Author:
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

import numpy as np

# Declaration of mathematical constants
g = 9.806           #[m/s2] Gravitational accelaration

def N2(rho,z):
    numdepths = len(z)
    N2z = np.zeros(numdepths -1)
    N2z[:] = np.nan

    # Calculate the gradient
    drhodz = np.gradient(rho,z,edge_order=1)
    # N2
    N2 = g/rho*drhodz

    # Calculate the new depth array
    for i in range(0,numdepths-1):
         N2z[i] = (z[i+1] + z[i])/2

    # return N2, N2z
    return N2, N2z
