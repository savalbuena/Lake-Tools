import numpy as np

def schmidtstability(z, rho, A, A0, z_bar_in, z_bar_input, g=9.806):
    """
    Compute the Stratification Stability Index (SSI).
    
    Parameters:
        z   : array-like, depths at layer centers (m)
        rho : array-like, densities at each depth (kg/m^3)
        Az  : array-like, horizontal area at each depth (m^2)
        dz  : float, layer thickness (m)
        A0  : float, lake surface area (m^2)
        g   : float, gravity (m/s^2)

    Returns:
        SS : float, stratification stability (kgm)
    """
    z = np.asarray(z)
    rho = np.asarray(rho)
    A = np.asarray(A)
    dz = z[1] - z[0]

    # Mean density
    V = np.trapezoid(A, z)
    
    if z_bar_in:
        zcv = z_bar_input
    else:
        zcv = np.trapezoid(z * A, z) / V

    # Compute SSI
    integrand = (z - zcv) * rho * A
    ss = g / A0 * np.trapezoid(integrand, z)
    return ss

def schmidt_stability_index(z, rho, A, A0, z_bar_in, z_bar_input, g=9.806):
    """
    Compute the Stratification Stability Index (SSI).
    
    Parameters:
        z   : array-like, depths at layer centers (m)
        rho : array-like, densities at each depth (kg/m^3)
        Az  : array-like, horizontal area at each depth (m^2)
        dz  : float, layer thickness (m)
        A0  : float, lake surface area (m^2)
        g   : float, gravity (m/s^2)

    Returns:
        SSI : float, stratification stability index (J/m^2)
    """
    z = np.asarray(z)
    rho = np.asarray(rho)
    A = np.asarray(A)
    dz = z[1] - z[0]

    # Mean density
    V = np.trapezoid(A, z)
    if z_bar_in:
        zcv = z_bar_input
    else:
        zcv = np.trapezoid(z * A, z) / V
    
    rho_mean = np.trapezoid(rho * A[:, None], z, axis=0) / V
    # Depth of center of gravity or where mean density is found
    drho = rho_mean - rho
    idx = np.nanargmin(np.abs(drho), axis=0)
    z_bar = z[idx]
    z = np.repeat(z[:, None], len(z_bar), axis=1)

    # Compute SSI
    integrand = (z - zcv) * (rho - rho_mean) * A[:,None]
    ssi = (g / A0) * np.trapezoid(integrand, z, axis=0)
    return ssi, (z - z_bar), (rho - rho_mean)