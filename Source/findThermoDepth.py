# findThermDepth.py
# ----Author: Jordan S Read 2009 ----
#   updated 3 march 2011
# Adapted to python and edited by:
# Sergio A Valbuena 2022
# removed signal processing toolbox function 'findpeaks.m' replaced with 'locPeaks.m' which was written to provide the same functionality

import numpy as np
import sys
from locPeaks import *


def findThermoDepth(**kw):
    nargin = len(kw)
    rhoVar = kw['rhoVar']
    depths = kw['depths']
    ax = kw['axis']
    print(ax)
    if nargin <= 3:
        Smin = 0.1
        seasonal = False
    elif nargin >= 5:
        seasonal = True

    dRhoPerc = 0.15   # min percentage max for unique thermocline step
    numDepths = len(depths)
    # drho_dz = np.empty(numDepths-1)*np.nan
    drho_dz = np.gradient(rhoVar, depths, axis=ax, edge_order=2)

    if seasonal:
        # look for two distinct maximum slopes, lower one assumed to be seasonal
        mDrhoZ = np.max(drho_dz)
        thermoInd, = np.where(drho_dz == mDrhoZ)
        thermoD = 1 / 2 * (depths[thermoInd] + depths[thermoInd + 1])
        if thermoInd > 0 and thermoInd < numDepths - 1:
            # if within range,
            Sdn = -(depths[thermoInd + 1] - depths[thermoInd]) /     \
                (drho_dz[thermoInd + 1] - drho_dz[thermoInd])
            Sup = (depths[thermoInd] - depths[thermoInd - 1]) /      \
                (drho_dz[thermoInd] - drho_dz[thermoInd - 1])
            upD = depths[thermoInd]
            dnD = depths[thermoInd + 1]
            if not np.any(np.isinf(Sup)) and not np.any(np.isinf(Sdn)):
                thermoD = dnD * (Sdn / (Sdn + Sup)) + upD * (Sup / (Sdn + Sup))

        dRhoCut = np.max([dRhoPerc * mDrhoZ, Smin])
        # [pks,locs] = locPeaks(drho_dz,dRhoCut);
        # if isempty(pks)
        #     SthermoD = thermoD;
        #     SthermoInd = thermoInd;
        # else
        #     mDrhoZ = pks(length(pks));
        #     SthermoInd = locs(length(pks));
        #     if SthermoInd > thermoInd+1
        #         SthermoD = mean([depths(SthermoInd)...
        #             depths(SthermoInd+1)]);
        #         if SthermoInd > 1 && SthermoInd < numDepths-1
        #             Sdn = -(depths(SthermoInd+1)-depths(SthermoInd))/...
        #                 (drho_dz(SthermoInd+1)-drho_dz(SthermoInd));
        #             Sup = (depths(SthermoInd)-depths(SthermoInd-1))/...
        #                 (drho_dz(SthermoInd)-drho_dz(SthermoInd-1));
        #             upD  = depths(SthermoInd);
        #             dnD  = depths(SthermoInd+1);
        #             if ~any([isinf(Sup) isinf(Sdn)])
        #                 SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
        #             end
        #         end
        #     else
        #         SthermoD = thermoD;
        #         SthermoInd = thermoInd;
        #     end
        # end
        # if SthermoD < thermoD;
        #     SthermoD = thermoD;
        #     SthermoInd = thermoInd;
        # end
        # return thermoD, thermoInd, drho_dz, SthermoD, SthermoInd

    else:
        mDrhoZ = np.max(drho_dz)
        thermoInd = np.where(drho_dz == mDrhoZ)[0]
        thermoD = 1 / 2 * (depths[thermoInd] + depths[thermoInd + 1])  # depth of max slope
        if thermoInd > 0 and thermoInd < numDepths - 1:
            # if within range,
            Sdn = -(depths[thermoInd + 1] - depths[thermoInd]) /     \
                (drho_dz[thermoInd + 1] - drho_dz[thermoInd])
            Sup = (depths[thermoInd] - depths[thermoInd - 1]) /      \
                (drho_dz[thermoInd] - drho_dz[thermoInd - 1])
            upD = depths[thermoInd]
            dnD = depths[thermoInd + 1]
            if not np.any(np.isinf(Sup)) and not np.any(np.isinf(Sdn)):
                thermoD = dnD * (Sdn / (Sdn + Sup)) + upD * (Sup / (Sdn + Sup))
        return thermoD, thermoInd, drho_dz
