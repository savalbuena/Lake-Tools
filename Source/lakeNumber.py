# lakeNumber.py LakeNumber estimation of a system following:
#     Ln = (g*St*(1-(ht/hm)))/(p0*(uStar^2)*(A0^1.5)*(1-hv/hm))
# ----Author: Jordan S Read 2009 ----
#  updated on 2 April 2010
#  updated on 25 January 2012
# Adapted to python and edited by:
# Sergio A. Valbuena 2022
#
# References:
#   -Imberger, Jorg, and John C. Patterson. "Physical Limnology." Advances in Applied Mechanics 27 (1990): 314-317.

# To import functions
import numpy as np

def lakeNumber(bthA,bthD,uSt,St,metaT,metaB,rhoHyp):
    g = 9.81
    dz = 0.1

    # if bathymetry has negative values, interpolate to 0
    if min(bthD) < 0:
        useI = bthD >= 0
        if not bthD >=0:
            depT = np.append(0,bthD)
        else:
            depT = bthD(useI)
        bthA = interp1(bthD,bthA,depT)
        bthD = depT
    else:
        bthA = bthA.to_numpy()
        bthD = bthD.to_numpy()

    Zo = np.min(bthD)
    Io, = np.where(bthD == Zo)[0]
    Ao = bthA[Io]
    if Ao == 0:
        error('Surface area cannot be zero, check *.bth file')

    # interpolates the bathymetry data
    layerD = np.arange(Zo,max(bthD)+dz,dz)
    layerA  = np.interp(layerD,bthD,bthA)

    # find depth to the center of volume
    Zv = layerD*layerA*dz
    Zcv = np.sum(Zv)/np.sum(layerA)/dz

    # Calculates the Lake Number according to the formula provided
    St_uC = St*Ao/g
    Ln = g*St_uC*(metaT+metaB)/(2*rhoHyp*uSt**2*Ao**(3/2)*Zcv)

    return Ln
