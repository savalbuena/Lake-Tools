# To import functions
def lakenumber(st, uSt, Ao, rho, zm, zt, zg):
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
    *Advances in Applied Mechanics*, 27, 303-475.

    Robertson, D. M., & Imberger, J., Boland, K. (1990). An indication of reservoir mixing: a water quality management tool.
    """

    g = 9.81
    st = st * Ao / g    # to transform Schmidt Stability Index from Idso to Stability Index
    zt = zm - zt        # Estimate height from bottom of thermocline
    zg = zm - zg        # Estimate height from bottom of center of volume

    Ln = (g * st * (1 - zt/zm)) / (rho * uSt**2 * Ao**(3/2) * (1 - zg/zm))

    return Ln
