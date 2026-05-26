"""
-----------------------------------------------------------------
hydrothermal.py
This script contains the functions used for STOL estimates related
to volumentric water temperature, deep temeprature, and stability
index for estimates of the stratification period.
It allows to obtain:
    - Mean annual temperature (Tavg)
    - Maximum annual temperature (Tmax)
    - Deep temperature at reference depth (T_ref_m)
----------------------------------------------------------------
Version 1.0
April 2025
Sergio A. Valbuena
"""
# Definition of functions to be used.
import numpy as np
import matplotlib.pyplot as plt

def compute_ssi(z, rho, vol, A0, z_bar_in, z_bar_input, g=9.81):
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

    # Mean density
    rho_mean = np.sum(rho * vol) / np.sum(vol)
    
    # Depth of center of volume (z̄)
    if z_bar_in:
        z_bar = z_bar_input
    else:
        z_bar = np.sum(rho * z * vol) / np.sum(rho * vol)
    # Compute SSI
    ssi = (g / A0) * np.sum((z - z_bar) * (rho - rho_mean) * vol)
    return ssi

def get_T_ref(z_ref, z, T):
    """
    Get temperature data at a specific depth.
    :param z_ref: Reference depth
    :param T: Temperature data
    :return: Temperature data at the reference depth
    """
    dz = np.abs(z - z_ref)
    iz = dz == np.min(dz)
    T_ref = T[iz, :]
    return T_ref

def get_T_profile(z_avg, z1, T1, z2, T2):
    """
    Get complete temperature profile from two sources of data.
    :param z_avg: Depths for new T_avg profile
    :param z1: Depths of reference profile
    :param T1: Temperatures of reference profile
    :param z2: Depths of profile to complete reference
    :param T2: Temperatures of profile to complete reference
    :return: Temperature profile of complete water column with
             depth z_avg
    """
    T_avg = np.interp(z_avg, z1, T1)
    iz = z_avg >= max(z1)
    T_avg[iz] = np.nan
    T_mltp_t = np.interp(z_avg, z2, T2)
    inan = np.isnan(T_avg)
    T_min = T_avg[~inan][-1]
    iT = T_mltp_t < T_min
    mask = inan & iT
    T_avg[mask] = T_mltp_t[mask]

    # To interpolate data and fill nan gaps for all water
    # column with depths z_avg.
    iT = ~np.isnan(T_avg)
    T_avg = np.interp(z_avg, z_avg[iT], T_avg[iT])
    return T_avg

def T_vol(T, vol):
    """
    Get volumentric average temperature from lake
    temperature data and lake volumen from
    hypsographic curves
    :param T: Temperature profile data
    :param vol: volume of layer per depth
    :return: Volumetric average temperature of water column
    :NOTE: Volume per layer and temperature profile must be
        in same order, descending from surface to bottom
    """
    Tvol_avg = np.nansum(T * vol) / np.nansum(vol)
    return Tvol_avg

def merge_profiles(z_out,
    z_shallow, rho_shallow,
    z_deep, rho_deep,
    extend="linear"  # "hold" or "linear"
):
    # Make a common target grid 0..zmax
    msk = ~np.isnan(rho_shallow)
    z_join = z_shallow[msk]
    z_join = np.max(z_join)

    # Interpolate both profiles onto the target grid where they exist
    rho_s = np.interp(
        z_out,
        np.asarray(z_shallow), np.asarray(rho_shallow),
        left=np.nan, right=np.nan
    )
    rho_d = np.interp(
        z_out,
        np.asarray(z_deep), np.asarray(rho_deep),
        left=np.nan, right=np.nan
    )

    # Start with deep profile
    rho_out = rho_d.copy()

    # Overwrite upper part with shallow profile (0..z_join)
    upper = z_out <= z_join
    rho_out[upper] = rho_s[upper]

    # Extend beyond max deep depth if needed (e.g., 460 -> 500)
    z_deep = np.asarray(z_deep)
    rho_deep = np.asarray(rho_deep)
    msk = ~np.isnan(rho_deep)
    zmax_deep = z_deep[msk]
    zmax_deep = np.nanmax(zmax_deep)

    beyond = z_out > zmax_deep
    if np.any(beyond):
        if extend == 'hold':
            msk = z_deep == zmax_deep
            rho_out[beyond] = rho_deep[msk][:]  # last value
        elif extend == 'linear':
            # Fit slope using last two deep points
            idx = np.argsort(z_deep)
            z_sorted = z_deep[idx]
            r_sorted = rho_deep[idx]
            z1, z2 = z_sorted[-2], z_sorted[-1]
            r1, r2 = r_sorted[-2], r_sorted[-1]
            slope = (r2 - r1) / (z2 - z1)
            rho_out[beyond] = r2 + slope * (z_out[beyond] - z2)
        else:
            raise ValueError("extend must be 'hold' or 'linear'")
    
    return rho_out

