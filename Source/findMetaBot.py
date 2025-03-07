# findMetaBot.py
# ----Author: Jordan S Read 2009 ----
# updated 12/15/2009 with thermoD pass
# Adapted for python and edited by:
# Sergio A Valbuena 2022


import numpy as np

def findMetaBot(drho_dz, thermoD, depths, slope):
    numDepths = len(depths)
    metaBot_depth = depths[numDepths - 1]
    Tdepth = np.empty(numDepths - 1) * np.nan
    for i in range(0, numDepths - 1):
        Tdepth[i] = 1/2*(depths[i+1] + depths[i])
    sortDepth0 = np.append(Tdepth,thermoD+1e-5)
    sortDepth = np.sort(sortDepth0)
    sortInd = np.argsort(sortDepth0)
    drho_dz = np.interp(sortDepth,Tdepth,drho_dz[:-1])
    thermo_index = 0
    thermoId = numDepths - 1
    for i in range(0,numDepths):
        if thermoId == sortInd[i]:
            thermo_index = i
            break

    for i in range(thermo_index,numDepths): #moving up from thermocline index
        if drho_dz[i] < slope: # bottom of metalimnion
            metaBot_depth = sortDepth[i]
            break

    #             % if interp can happen
    if i-thermo_index > 0 and drho_dz[thermo_index] > slope:
        metaBot_depth = np.interp(slope, drho_dz[thermo_index:i+1:1], sortDepth[thermo_index:i+1:1])

    if np.isnan(metaBot_depth):
        metaBot_depth = np.max(depths)

    return metaBot_depth
