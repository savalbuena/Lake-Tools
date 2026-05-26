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
        bthA = np.interp(bthD,bthA,depT)
        bthD = depT
    else:
        bthA = bthA
        bthD = bthD

    Zo = np.min(bthD)
    Io, = np.where(bthD == Zo)[0]
    Ao = bthA[Io]
    if Ao == 0:
        raise ValueError('Surface area cannot be zero, check *.bth file')

    # interpolates the bathymetry data
    layerD = np.arange(Zo,max(bthD)+dz,dz)
    layerA  = np.interp(layerD,bthD,bthA)

    # find depth to the center of volume
    Zv = layerD*layerA*dz
    Zcv = np.sum(Zv)/np.sum(layerA*dz)

    print(Zcv)

    # Calculates the Lake Number according to the formula provided
    St_uC = St * Ao / g
    Ln = g * St_uC * (metaT + metaB) / (2 * rhoHyp * uSt**2 * Ao ** (3 / 2) * Zcv)
    return Ln

def lakeNumber_v2(st, uSt, Ao, rho, zm, zt, zg):
    """
    Compute the Lake Number (Ln), a dimensionless index describing the
    relative importance of stratification versus wind forcing in a lake.

    The Lake Number compares the stabilizing effect of density
    stratification (quantified by Schmidt stability) to the destabilizing
    effect of wind stress acting on the lake surface. It is commonly used
    to assess whether wind forcing is sufficient to tilt the thermocline,
    induce internal seiches, or cause upwelling and vertical mixing.

    Physically, the Lake Number represents the ratio of restoring buoyancy
    forces to wind-driven shear forces.

        Ln >> 1  → stratification dominates (thermocline stable)
        Ln ~ 1   → wind and stratification are comparable
        Ln << 1  → wind dominates (mixing / upwelling likely)

    Parameters
    ----------
    st : float
        Schmidt stability [J m⁻²]. Represents the potential energy per unit
        surface area required to fully mix the lake to uniform density.

    uSt : float
        Water-side friction velocity [m s⁻¹], derived from wind stress at
        the lake surface. This term represents the intensity of wind
        forcing acting on the water column.

    Ao : float
        Lake surface area [m²].

    rho : float
        Representative water density [kg m⁻³], typically taken as the
        volumetric mean density of the epilimnion of the lake.

    zm : float
        Maximum depth of the lake [m]

    zt : float
        Depth of the thermocline [m], measured positive downward from the
        lake surface. This represents the depth of maximum density (or
        temperature) gradient separating the epilimnion and hypolimnion.

    zg : float
        Depth to the center of volume of the lake [m], measured positive
        downward from the lake surface.

    Returns
    -------
    Ln : float
        Lake Number (dimensionless). Higher values indicate a more stable,
        stratification-dominated system, while lower values indicate
        stronger wind influence and increased likelihood of mixing or
        upwelling.

    Notes
    -----
    - The Lake Number formulation assumes horizontally uniform density
      structure and neglects Coriolis effects.
    - Accurate estimation of Schmidt stability and thermocline depth is
      critical for reliable Lake Number estimates.
    - This formulation follows classical limnological definitions and is
      consistent with implementations used in LakeAnalyzer and related
      lake stability studies.

    References
    ----------
    Imberger, J., & Patterson, J. C. (1990). Physical limnology.
    *Advances in Applied Mechanics*, 27, 303–475.

    Robertson, D. M., & Imberger, J., Boland, K. (1990). An indication of reservoir mixing: a water quality management tool.
    """

    g = 9.81
    st = st * Ao / g    # to transform Schmidt Stability Index from Idso to Stability Index
    zt = zm - zt        # Estimate height from bottom of thermocline
    zg = zm - zg        # Estimate height from bottom of center of volume
    Ln = (g * st * (1 - zt/zm)) / (rho * uSt**2 * Ao**(3/2) * (1 - zg/zm))

    return Ln
