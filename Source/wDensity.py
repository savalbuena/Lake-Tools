# WaterDensity.py Density of water with supplied Temp and Salinity
# rho = wDensity(Temperature, Salinity,Pressure,pressureType,method) is the density of water with the given temperature, salinity, and/or pressure
# The Temperature and method input are mandatory, while Salinity, Pressure, and pressureType are optional. If Pressure is used, it must be accompanied by the pressureType.
# T is in °C
# S is in Pratical Salinity Scale units (dimensionless)
# P is in bars for pressureType = 'data'
# P is in m for pressureType = 'depth'
# pressureType = 'data' is for pressure data in bars
# pressureType = 'depth' uses the input as depth and uses hydrostatic approximation to estimate the pressure. A root finding tool is used then to estimate the density.
# method = 'UN' uses the UNESCO equation
# method = 'MM' uses the Martin & MacCutcheon 1999 equation rho(T)
# T, S, and P should be the same size, unless S is a scalar
# rho is in kg/m3
#
# <<<--- Effective range: 0-40°C, 0.5-43 Salinity --->>>
#
#----Author: Sergio A Valbuena 2021 ----
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

# Example of usage
# Example 1
# Temp = np.arange(24,4,-1)
# Pres = np.arange(0,100,5)
# rho = wDensity(T=Temp,P=Pres,pressureType='depth',method='UN')
# Example 2
# Temp = np.arange(24,4,-1)
# rho = wDensity(T=Temp,method='UN')
# Example 3
# Temp = np.arange(24,4,-1)
# Pres = np.arange(0,100,5)
# Sal = np.arange(0,21,1)
# rho = wDensity(T=Temp,P=Pres,S=Sal,pessureType='depth',method='UN')
# Example 4
# Temp = np.arange(24,4,-1)
# Sal = np.arange(0,21,1)
# rho = wDensity(T=Temp,S=Sal,method='UN')

import sys
import os
import numpy as np

def wDensity(**kw):
    keys = list(kw.keys())
    if 'T'in keys and not 'S' in keys:
        T = kw['T']
        S = 0
    elif 'T' in keys and 'S' in keys:
        T = kw['T']
        S = kw['S']

    if kw['method'] == 'MM':
        # <<equation provided by: Martin, J.L., McCutcheon, S.C., 1999. Hydrodynamics and Transport for Water Quality Modeling. Lewis Publications, Boca Raton, FL, 794pp.>>
        rho = (1000*(1-(T+288.9414)*(T-3.9863)**2/(508929.2*(T+68.12963))))
    elif kw['method'] == 'UN':
        # <<equations provided by: Millero, F.J., Poisson, A., 1981. International one-atmosphere equation of state of seawater. UNESCO Technical Papers in Marine Science. No. 36>>
        #--eqn (1):
        rho_0 = 999.842594+6.793952e-2*T-9.095290e-3*T**2+1.001685e-4*T**3-1.120083e-6*T**4+6.536335e-9*T**5
        # --eqn (2):
        A = 8.24493e-1-4.0899e-3*T+7.6438e-5*T**2-8.2467e-7*T**3+5.3875e-9*T**4
        # --eqn (3):
        B = -5.72466e-3+1.0227e-4*T-1.6546e-6*T**2
        # --eqn (4):
        C = 4.8314e-4
        # --eqn (5):
        rho_0 += A*S+B*S**(3/2)+C*S**2

        if 'P' in keys:
            P = kw['P']
            if kw['pressureType'] == 'depth':
                maxiter = 10000
                rhoguess = rho_0
                delta = 1
                iter = 0
                while np.mean(delta) > 1e-8 and iter < maxiter:
                    pressureh = rhoguess*9.806*P;
                    Pre = 1e-5*pressureh;
                    kw = 19652.21 + 148.4206*T - 2.327105*T**2           \
                        + 1.360477e-2*T**3 - 5.155288e-5*T**4

                    k = kw + S*(54.6746 - 0.603459*T + 1.09987e-2*T**2   \
                        - 6.1670e-5*T**3) \
                        + S**(3/2)*(7.944e-2 + 1.6483e-2*T - 5.3009e-4*T**2)

                    k = k + Pre*(3.239908 + 1.43713e-3*T + 1.16092e-4*T**2 \
                        - 5.77905e-7*T**3) \
                        + Pre*S*(2.2838e-3 - 1.0981e-5*T - 1.6078e-6*T**2) \
                        + 1.91075e-4*Pre*S**(3/2) \
                        + Pre**2*(8.50935e-5 - 6.12293e-6*T + 5.2787e-8*T**2) \
                        + Pre**2*S*(-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T**2)

                    rho = rho_0/(1-Pre/k)
                    delta = abs(rho - rhoguess)
                    rhoguess = rho
                    iter = iter + 1
            elif kw['pressureType'] == 'data':
                kw = 19652.21 + 148.4206*T - 2.327105*T**2           \
                    + 1.360477e-2*T**3 - 5.155288e-5*T**4

                k = kw + S*(54.6746 - 0.603459*T + 1.09987e-2*T**2   \
                    - 6.1670e-5*T**3) \
                    + S**(3/2)*(7.944e-2 + 1.6483e-2*T - 5.3009e-4*T**2)

                k = k + P*(3.239908 + 1.43713e-3*T + 1.16092e-4*T**2 \
                    - 5.77905e-7*T**3) \
                    + P*S*(2.2838e-3 - 1.0981e-5*T - 1.6078e-6*T**2) \
                    + 1.91075e-4*P*S**(3/2) \
                    + P**2*(8.50935e-5 - 6.12293e-6*T + 5.2787e-8*T**2) \
                    + P**2*S*(-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T**2)

                rho = rho_0/(1-P/k)
        else:
            rho = rho_0
    return rho
