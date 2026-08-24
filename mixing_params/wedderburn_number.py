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

def wedderburn(rho_epi, rho_hypo, h_1, uSt, A0, Lo):
    # Length at thermocline depth
    if Lo is None:
        L = 2 * (A0 / pi)**(1 / 2)
    else:
        L = Lo
    
    gp = g * (rho_hypo - rho_epi) / rho_epi
    # Calculates W according to formula provided
    W = gp * h_1**2 / (uSt**2 * L)

    return W
