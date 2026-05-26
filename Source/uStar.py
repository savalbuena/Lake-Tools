# uStar.py estimates the wind shear velocity
#----Author: Jordan S Read 2009 ----
# **** updated 19 Feb 2010 ****
# Adapted for python and edited by:
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

# Calculates the Wedderburn number for a particular system using the following equation:
#   ustar = (tau/)**(1/2)
#   tau = cd*rhoAir*wspd**2;
#
# where
#   cd = Wind drag coefficient
#   wspd = wind speed (m/s)
#   windHeight = height of wind measurement (m)
#   tau = surface wind shear stress due to wind Pa
#   uSt = water friction velocity due to wind stress
#   EpiRho = mean epilimniom density (kg/m3)

#equation (1) provided by:
#°° Hicks, B.B., 1972. A procedure for the formulation of bulk transfer coefficients over water bodies of different sizes. Boundary-Layer Meterology 3: 201-213

#equation (2) provided by:
#°° Amorocho, J., DeVries, J.J., 1980. A new evaluation of the wind stress coefficient over water surfaces. Journal of Geophysical Research 85: 433-442.

#equation (3) provided by:
#°° Fischer, H.B., List, E.J., Koh, R.C.Y., Imberger, J., Brooks, N.H., 1979. Mixing in inland and coastal waters. Academic Press.

#equation (4) provided by:
#°° Imberger, J., 1985. The diurnal mixed layer. Limnology and Oceanography 30: 737-770.

# Library Import
import numpy as np


rhoAir = 1.2;
vonK = 0.4; #von Karman k

def uStar(wspd,windHeight,EpiRho):
    # ~ eqn (1) ~
    wspd = np.atleast_1d(wspd)
    cd = 0.0015*np.ones(len(wspd))
    idata = wspd < 5
    cd[idata] = 0.001

    # ~ eqn (2) ~
    #If the windspeed is not measured at a height of 10 meters, corrected by
    #(eqn 21) from Amorocho and DeVries, 1980
    if windHeight != 10:
        wspd = wspd/(1-np.sqrt(cd)/vonK*np.log(10/windHeight));
    # ~ eqn (3) ~
    tau = cd*rhoAir*wspd**2;
    # ~ eqn (4) ~
    ustar = np.sqrt(tau/EpiRho);

    return ustar
