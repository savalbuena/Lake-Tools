# schmidtStability.py
# ----Author: Jordan S Read 2009 ----
#  updated 7 March 2011, adding salinity
#  updated 20 Nov 2011, adding timeseries bathymetric effects
#  Adapted for python and edited by:
# Sergio A Valbuena 2022
# equation provided by:
# °° Idso, S.B., 1973. On the concept of lake stability. °°
# °° Limnol. Oceanogr. 18: 681–683.                      °°
#
#    St = (g/A0)* [Integral from 0 to hm: (hv-h)*A(h)*p(h)*dh]
#
#  The Lake Analyzer program reverses the direction of this calculation so that z = 0 is the surface, and positive z is downwards

# schmidtStability(wtr,depths,bthA,bthD,sal):

import numpy as np
from wDensity import *

def schmidtStability(**kw):
    nargin = len(kw)

    if len(kw['wtr']) != len(kw['depths']):
        print(' ERROR! water temperature array must be the same length as the depth array')
        exit()
    elif nargin < 4:
        print('ERROR! not enough input arguments')
        exit()
    elif np.any(np.isnan(kw['wtr'])) or np.any(np.isnan(kw['depths'])):
        print('Error! input arguments must be numbers')
        exit()
    if nargin == 4:
        kw['sal'] = kw['wtr']*0

    g = 9.806   # [m s^-2]
    dz = 0.1   # [m]
    # if bathymetry has negative values, interpolate to 0
    if min(kw['bthD']) < 0:
        useI = kw['bthD'] >= 0
        if not kw['bthD'] >=0:
            depT = np.append(0,kw['bthD'])
        else:
            depT = kw['bthD'][useI]
        kw['bthA'] = np.interp(depT,kw['bthD'],kw['bthA'])
        kw['bthD'] = depT
    else:
        kw['bthA'] = kw['bthA']
        kw['bthD'] = kw['bthD']

    numD = len(kw['wtr'])
    if np.max(kw['bthD']) > kw['depths'][numD-1]:
        kw['wtr'] = np.append(kw['wtr'],kw['wtr'][numD-1])
        kw['sal'] = np.append(kw['sal'],kw['sal'][numD-1])
        kw['depths'] = np.append(kw['depths'],np.max(kw['bthD']))
    elif np.max(kw['bthD']) < kw['depths'][numD-1]:
        kw['bthD'] = np.append(kw['bthD'],kw['depths'][numD-1])
        kw['bthA'] = np.append(kw['bthA'],0)

    if np.min(kw['bthD']) < kw['depths'][0]:
        kw['wtr'] = np.append(kw['wtr'][0],kw['wtr'])
        kw['sal'] = np.append(kw['sal'][0],kw['sal'])
        kw['depths'] = np.append(np.min(kw['bthD']),kw['depths'])

    Zo = np.min(kw['bthD'])
    Zm = np.max(kw['depths'])
    Io, = np.where(kw['bthD'] == Zo)[0]
    Ao = kw['bthA'][Io]
    if Ao == 0:
        raise ValueError('Surface area cannot be zero, check *.bth file')

    rhoL = wDensity(T=kw['wtr'],P=kw['depths'],S=kw['sal'], pressureType = 'depth')

    # interpolates the bathymetry data
    layerD = np.arange(Zo,Zm+dz,dz)
    layerP  = np.interp(layerD,kw['depths'],rhoL)
    layerA = np.interp(layerD,kw['bthD'],kw['bthA'])

    # find depth to the center of volume
    Zv = layerD*layerA*dz
    Zcv = np.sum(Zv)/np.sum(layerA)/dz

    numInt = len(layerA)
    st = np.empty(numInt)*np.nan
    for i in range(0,numInt):
        z = layerD[i]
        A = layerA[i]
        st[i] = -(Zcv-z)*layerP[i]*A*dz
    St = g/Ao*np.sum(st)

    return St
