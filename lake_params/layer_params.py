# layerTemperature.py
#----Author: Jordan S Read, 2013----
# updated 20 Nov 2013, adding timeseries bathymetric effects

# Finds the average temperature of a thermal layer, bounded by "top" and "bottom" where distances are measured from the surface of the water column.

#Input:
#   -temps: temperature values (celsius)
#   -depths: corresponding depth values (m)
#   -metaDef: critical slope used to define the metalimnion
#   -top: surface to top of metalimnion
#   -bottom: surface to bottom of metalimnion

#Output:
#   -aveTemperature: average temperature the layer bounded by top & bottom

import numpy as np

def layer_avg(top,bottom,wtr,depths,bthA,bthD):

    if top > bottom:
        raise ValueError('ERROR! bottom depth must be greater than top')
    if len(wtr) != len(depths):
        raise ValueError(' ERROR! water temperature array must be the same length as the depth array')
    elif np.any(np.isnan(wtr)) or np.any(np.isnan(depths)):
        raise ValueError('Error! input arguments must be numbers')

    # if bathymetry has negative values, interpolate to 0
    if min(bthD) < 0:
        useI = bthD >= 0
        if not bthD >=0:
            depT = np.append(0,bthD)
        else:
            depT = bthD[useI]
        bthA = np.interp(depT,bthD,bthA)
        bthD = depT
    else:
        bthA = bthA
        bthD = bthD
    dz = 0.0001
    numD = len(wtr)
    if np.max(bthD) > depths[numD-1]:
        wtr = np.append(wtr,wtr[numD-1])
        depths = np.append(depths,np.max(bthD))
    elif np.max(bthD) < depths[numD-1]:
        bthD = np.append(bthD,depths[numD-1])
        bthA = np.append(bthA,0)

    if np.min(bthD) < depths[0]:
        wtr = np.append(wtr[0],wtr)
        depths = np.append(np.min(bthD),depths)

    Zo = np.min(bthD)
    Zm = np.max(depths)
    Io, = np.where(bthD == Zo)[0]
    Ao = bthA[Io]
    if Ao == 0:
        raise ValueError('Surface area cannot be zero, check *.bth file')

    # interpolates the bathymetry data
    layerD = np.arange(top,bottom+dz,dz)
    layerT = np.interp(layerD,depths,wtr)
    layerA = np.interp(layerD,bthD,bthA)

    weightedT = layerA*layerT*dz;
    aveTemperature = sum(weightedT)/(sum(layerA))/dz

    return aveTemperature


def layer_avg_v2(top,bottom,wtr,depths,vol, bthD):

    if top > bottom:
        raise ValueError('ERROR! bottom depth must be greater than top')
    if len(wtr) != len(depths):
        raise ValueError(' ERROR! water temperature array must be the same length as the depth array')
    elif np.any(np.isnan(wtr)) or np.any(np.isnan(depths)):
        raise ValueError('Error! input arguments must be numbers')

    # if bathymetry has negative values, interpolate to 0
    if min(bthD) < 0:
        useI = bthD >= 0
        if not bthD >=0:
            depT = np.append(0,bthD)
        else:
            depT = bthD[useI]
        vol = np.interp(depT,bthD,vol)
        bthD = depT
    else:
        vol = vol.to_numpy()
        bthD = bthD
    dz = 0.0001
    numD = len(wtr)
    if np.max(bthD) > depths[numD-1]:
        wtr = np.append(wtr,wtr[numD-1])
        depths = np.append(depths,np.max(bthD))
    elif np.max(bthD) < depths[numD-1]:
        bthD = np.append(bthD,depths[numD-1])
        vol= np.append(vol,0)

    if np.min(bthD) < depths[0]:
        wtr = np.append(wtr[0],wtr)
        depths = np.append(np.min(bthD),depths)

    Zo = np.min(bthD)
    Zm = np.max(depths)
    Io, = np.where(bthD == Zo)[0]
    Ao = vol[Io] / dz
    if Ao == 0:
        raise ValueError('Surface area cannot be zero, check *.bth file')

    # interpolates the bathymetry data
    layerD = np.arange(top,bottom+dz,dz)
    layerT = np.interp(layerD,depths,wtr)
    # layerA = np.interp(layerD,bthD,bthA)
    layer_vol = np.interp(layerD, bthD, vol)

    weightedT = layer_vol*layerT;
    aveTemperature = sum(weightedT)/(sum(layer_vol))

    return aveTemperature
