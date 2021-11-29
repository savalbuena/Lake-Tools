# Wedderburn.py Estimates the Wedderburn parameter
# ----Author: Jordan S Read 2009 ----
# updated 2 April 2010
# Adapted for python and edited by:
# 05-11-2021 - Sergio A Valbuena <savalbuena@ucdavis.edu>

# Calculates the Wedderburn number for a particular system using the following equation:
#
#   W = (g*delta_rho*(h^2))/(pHyp*(uSt^2)*Lo)
#
# where
#   g = force of gravity
#   delta_rho = density difference between the epilimnion and the hypolimnion
#   metaT = thickness of the surface layer
#   uSt = water friction velocity due to wind stress
#   Lo = fetch length in the direction of the wind.
#
# Reference:
# °° Imberger, Jorg, and John C. Patterson. "Physical Limnology." Advances in Applied Mechanics 27 (1990): 314-317.

import numpy as np
import math

g = 9.806
pi = math.pi

def WedderburnNumber(delta_rho,metaT,uSt,A0,AvHyp_rho):
    # Length at thermocline depth
    Lo = 2*(A0/pi)**(1/2)
    go = g*delta_rho/AvHyp_rho
    # Calculates W according to formula provided
    W = go*metaT**2/(uSt**2*Lo)
    return W
