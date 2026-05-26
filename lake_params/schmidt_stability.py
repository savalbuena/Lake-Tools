import numpy as np

def schmidt_stability(z, rho, A, A0, z_bar_in, z_bar_input, g=9.81):
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
    V = np.trapz(A, z)
    
    if z_bar_in:
        zcv = z_bar_input
    else:
        zcv = np.trapz(z * A, z) / V

    # Compute SSI
    integrand = (z - zcv) * rho * A
    ss = g / A0 * np.trapz(integrand, z)
    return ss

def schmidt_stability_index(z, rho, A, A0, z_bar_in, z_bar_input, g=9.81):
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
    V = np.trapz(A, z)
    if z_bar_in:
        zcv = z_bar_input
    else:
        zcv = np.trapz(z * A, z) / V
    
    rho_mean = np.trapz(rho * A, z) / V

    # Depth of center of gravity or where mean density is found
    drho = rho_mean - rho
    idx = np.nanargmin(np.abs(drho))
    z_bar = z[idx]

    # Compute SSI
    integrand = (z - zcv) * (rho - rho_mean) * A
    ssi = (g / A0) * np.trapz(integrand, z)
    # ssi = (g / A0) * np.sum((z - z_bar_input) * (rho - rho_mean) * A * dz)
    return ssi, (z - z_bar), (rho - rho_mean)


# def schmidt_stability(z, rho, vol, A0, z_bar_in, z_bar_input, g=9.81):
#     """
#     Compute the Stratification Schmidt Stability (SS).
    
#     Parameters:
#         z   : array-like, depths at layer centers (m)
#         rho : array-like, densities at each depth (kg/m^3)
#         Az  : array-like, horizontal area at each depth (m^2)
#         dz  : float, layer thickness (m)
#         A0  : float, lake surface area (m^2)
#         g   : float, gravity (m/s^2)

#     Returns:
#         SS : float, stratification stability index (J/m^2)
#     """
#     z = np.asarray(z)
#     rho = np.asarray(rho)
#     vol = np.asarray(vol)

#     # Mean density
#     rho_mean = np.sum(rho * vol) / np.sum(vol)

#     drho = rho_mean - rho
#     idx = np.nanargmin(np.abs(drho))
#     z_bar = z[idx]
    
#     # Depth of center of volume (z̄)
#     # if z_bar_in:
#     #     z_bar = z_bar_input
#     # else:
#     #     z_bar = np.sum(rho * z * vol) / np.sum(rho * vol)
#     # Compute SSI
#     ss = (g / A0) * np.sum((z - z_bar) * (rho - rho_mean) * vol)
#     return ss, (z-z_bar), (rho - rho_mean)