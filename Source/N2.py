# Calculates buoyancy frequency.

# def N2(rho,z,axis)
# rho is the water density and can be a vector or a 2D matrix
# z is the depths at which info of water density is available. Must be same size as rho
# axis = 1,0 to indicate if the values of z change in the row direction (axis = 0) or in the column direction (axis = 1)

# Description Calculate the buoyancy frequency (Brunt-Vaisala frequency) for a density profile.
# Param rho a numeric vector of water density previously obtained with wdensity() function in kg/m3
# Param depths a numeric vector corresponding to the depths (in m) of the rho measurements
# Return a vector of buoyancy frequency in units sec^-2.

# Author:
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

import numpy as np

# Declaration of mathematical constants
g = 9.806           # [m/s2] Gravitational accelaration


def N2(rho, z, axis):

    drhodz = np.gradient(rho, z, axis=axis, edge_order=2)
    N2 = (g / rho) * (drhodz)

    return N2, drhodz
