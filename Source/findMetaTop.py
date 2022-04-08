#findMetaTop.py
# ----Author: Jordan S Read 2009 ----
# updated 12/16/2009 with thermoD pass
# Adapted for python and edited by:
# Sergio A Valbuena 2022

import numpy as np

def findMetaTop(drho_dz,thermoD,depths,slope):
    numDepths = len(depths)
    metaTop_depth = 1/2*(depths[1] + depths[0])
    Tdepth = np.empty(numDepths-1)*np.nan
    for i in range(0,numDepths-1):
        Tdepth[i] = 1/2*(depths[i+1] + depths[i])
    sortDepth0 = np.append(Tdepth,thermoD - 1e-5)
    sortDepth = np.sort(sortDepth0)
    sortInd = np.argsort(sortDepth0)
    drho_dz = np.interp(sortDepth,Tdepth,drho_dz[:-1])
    thermo_index = 0
    thermoId = numDepths - 1
    for i in range(0,numDepths):
        if thermoId == sortInd[i]:
            thermo_index = i
            break

    for i in range(thermo_index,0,-1): #moving up from thermocline index
        if drho_dz[i] < slope: # top of metalimnion
            metaTop_depth = sortDepth[i]
            break

    #             % if interp can happen
    if thermo_index-i > 0 and drho_dz[thermo_index] > slope:
        metaTop_depth = np.interp(slope, drho_dz[i:thermo_index], sortDepth[i:thermo_index])

    if np.isnan(metaTop_depth):
        metaTop_depth = np.min(depths)
    if metaTop_depth > thermoD:
        metaTop_depth = thermoD

    return metaTop_depth
