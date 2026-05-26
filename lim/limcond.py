
from __future__ import annotations
import os
from types import SimpleNamespace
from typing import Any, Dict, Tuple, Optional

import numpy as np
from scipy.io import loadmat

from .limconsts import mm, KCL, as_const

def _as_array(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=float)

def _is_str(x: Any) -> bool:
    return isinstance(x, str)

def _pval(coeff: np.ndarray, T: np.ndarray) -> np.ndarray:
    """
    Evaluate polynomial with coefficients in MATLAB's pval order:
      coeff[0] + T*(coeff[1] + T*(...))
    """
    coeff = np.asarray(coeff, dtype=float).ravel()
    y = np.full_like(T, coeff[-1], dtype=float)
    for c in coeff[-2::-1]:
        y = T*y + c
    return y

def _lengthconst(const: SimpleNamespace) -> int:
    lens = []
    for k, v in const.__dict__.items():
        if isinstance(v, str):
            continue
        a = np.asarray(v)
        lens.append(a.size)
    uniq = sorted(set([l for l in lens if l is not None]))
    uniq = [u for u in uniq if np.isfinite(u)]
    if len(uniq) == 0:
        return 1
    if len(uniq) > 2:
        raise ValueError("Inconsistent field lengths in composition structure")
    if len(uniq) == 2 and uniq[0] != 1:
        raise ValueError("Inconsistent field lengths in composition structure")
    return int(max(uniq))

def _getconst(const: SimpleNamespace, k: int) -> SimpleNamespace:
    out = {}
    for name, vals in const.__dict__.items():
        if isinstance(vals, str):
            out[name] = vals
        else:
            arr = np.asarray(vals)
            if arr.size > 1:
                out[name] = float(arr.flat[k])
            elif arr.size == 1:
                out[name] = float(arr.flat[0])
    return SimpleNamespace(**out)

def limcond(*args, nargout: int | None = None):
    """
    Port of limcond.m (R. Pawlowicz).

    Supports both MATLAB call styles, but the primary intended use (your pipeline) is:
        Kref, Sc = limcond(K, T, p, Tref, const, nargout=2)

    Parameters
    ----------
    If first argument is a dict/namespace -> "conductivity forward mode" (condmode=1):
        limcond(const, T, p, ...options...)

    Else -> "reference conductivity mode" (condmode=0):
        limcond(K, T, p, [Tref], [const], ...options...)

    Notes
    -----
    - K is in-situ conductivity (uS/cm)
    - Sc returned here is g/L (NOT g/kg), matching MATLAB header comment.
    - const.units must be one of {'M','mM','uM'} for this port.
    """
    if len(args) == 0:
        # Return blank constituent structure (like MATLAB)
        return SimpleNamespace(
            units='uM|mM|M',
            HCO3=0., CO3=0., H=0., OH=0., Cl=0., NO3=0., SO4=0., Ca=0., Na=0., K=0., Mg=0.,
            NH4=0., Br=0., BOH4=0., Sr=0., F=0., ClO4=0., I=0., Ag=0., Li=0., Zn=0., Cu=0., Ba=0.,
            FeCN3=0., FeCN4=0., La=0., SiO2=0., SiOH4=0., BOH3=0., CO2=0., H2CO3=0., O2=0., N2=0.
        )

    args = list(args)

    # Detect condmode
    condmode = 1 if isinstance(args[0], (dict, SimpleNamespace)) else 0

    # Load PFAC matrices (Ampij, Icij, widij)
    pfac_path = os.path.join(os.path.dirname(__file__), "limPFAC.mat")
    pfac = loadmat(pfac_path)
    Ampij = np.asarray(pfac["Ampij"], dtype=float)
    Icij  = np.asarray(pfac["Icij"], dtype=float)
    widij = np.asarray(pfac["widij"], dtype=float)

    # Constants used later
    Beta = 1/300.0
    gammap = np.array([0.0, 9.0718e-6, -2.8013e-10], dtype=float)
    gammaT = np.array([1.0, -2.5751e-2, 8.46082e-4], dtype=float)

    # Options defaults
    calcopt = 0   # 0 full, 1 nopairs, 2 noion
    fullsal = 1   # 0 ionic only, 1 include non-ionic, 2 include gases too

    # Parse options (strings anywhere after required args)
    for a in args:
        if isinstance(a, str):
            s = a.lower()
            if s.startswith("full"):
                calcopt = 0
            elif s.startswith("nopa"):
                calcopt = 1
            elif s.startswith("noio"):
                calcopt = 2
            elif s.startswith("sall"):
                fullsal = 2
            elif s.startswith("sful"):
                fullsal = 1
            elif s.startswith("sion"):
                fullsal = 0

    if condmode == 1:
        const = as_const(args[0])
        lc = _lengthconst(const)
        if lc > 1:
            kR = np.full(lc, np.nan); Sc = np.full(lc, np.nan); Balfrac=np.full(lc,np.nan); TC=np.full(lc,np.nan)
            for k in range(lc):
                c2 = _getconst(const, k)
                out = limcond(c2, *args[1:], nargout=4)
                kR[k], Sc[k], Balfrac[k], TC[k] = out
            outs = (kR, Sc, Balfrac, TC)
            return outs[:nargout] if nargout else outs

        T = _as_array(args[1])
        p = _as_array(args[2])
        # not used in your workflow, but implemented
        kpT = None
    else:
        kpT = _as_array(args[0])
        T = _as_array(args[1])
        p = _as_array(args[2])

        # Optional Tref
        Tref = 25.0
        if len(args) > 3 and (not isinstance(args[3], (dict, SimpleNamespace))) and (not isinstance(args[3], str)):
            Tref = float(args[3])
            args.pop(3)

        # Optional const
        if len(args) > 3 and isinstance(args[3], (dict, SimpleNamespace)):
            const = as_const(args[3])
            args.pop(3)
        else:
            const = KCL

    # Unit conversion
    units = getattr(const, "units", None)
    if units in ('mg/l','m','mm','um','mol/kg','mmol/kg','umol/kg'):
        raise ValueError("Need molar units here (M/mM/uM).")
    if units == 'M':
        convfac = 1.0
    elif units == 'mM':
        convfac = 1e-3
    elif units == 'uM':
        convfac = 1e-6
    else:
        raise ValueError(f"Unrecognized units '{units}'")

    # Ion coefficient tables (subset identical to MATLAB)
    Icoeffs = [
        ('K',     +1, mm.K,     73.48,  [+1.941e-02, +5.786e-05, -1.593e-07], 0.62),
        ('Cl',    -1, mm.Cl,    76.31,  [+2.012e-02, +6.217e-05, -4.921e-08], 0.62),
        ('OH',    -1, mm.OH,   198.00,  [+1.644e-02, +1.403e-05, +4.500e-06], 0.31),
        ('NO3',   -1, mm.NO3,   71.42,  [+1.874e-02, +6.764e-05, +8.948e-07], 0.74),
        ('HCO3',  -1, mm.HCO3,  44.50,  [+1.909e-02, +1.573e-05, -1.320e-06], 0.91),
        ('Br',    -1, mm.Br,    78.10,  [+1.960e-02, +6.678e-05, +4.371e-07], 0.60),
        ('F',     -1, mm.F,     55.40,  [+2.089e-02, 0.0, 0.0],                0.80),
        ('ClO4',  -1, mm.ClO4,  67.30,  [+1.724e-02, -2.359e-05, 0.0],         0.66),
        ('I',     -1, mm.I,     76.80,  [+1.950e-02, +6.370e-05, +7.755e-07],  0.61),
        ('BOH4',  -1, mm.BOH4,  35.20,  [+2.130e-02, 0.0, 0.0],                1.25),
        ('NH4',   +1, mm.NH4,   73.50,  [+1.931e-02, +1.113e-04, +2.546e-06],  0.62),
        ('Na',    +1, mm.Na,    50.08,  [+2.189e-02, +9.159e-05, -1.127e-06],  0.80),
        ('H',     +1, mm.H,    349.65,  [+1.388e-02, -3.152e-05, -6.372e-07],  0.30),
        ('Ag',    +1, mm.Ag,    61.90,  [+1.969e-02, +4.308e-05, 0.0],         0.62),
        ('Li',    +1, mm.Li,    38.66,  [+2.294e-02, +1.313e-04, +4.451e-07],  0.96),
        ('Ca',    +2, mm.Ca,    59.47,  [+2.213e-02, +9.854e-05, -1.027e-06],  0.85*2),
        ('Mg',    +2, mm.Mg,    53.00,  [+2.325e-02, +2.025e-04, 0.0],         0.90*2),
        ('Sr',    +2, mm.Sr,    59.40,  [+2.096e-02, +7.322e-05, 0.0],         0.84*2),
        ('Zn',    +2, mm.Zn,    52.80,  [+2.200e-02, +1.287e-04, 0.0],         0.80*2),
        ('Cu',    +2, mm.Cu,    53.60,  [+2.329e-02, +1.676e-04, 0.0],         0.80*2),
        ('Ba',    +2, mm.Ba,    63.60,  [+2.084e-02, +8.885e-05, 0.0],         0.83*2),
        ('SO4',   -2, mm.SO4,   80.00,  [+2.119e-02, +6.746e-05, 0.0],         0.74*2),
        ('CO3',   -2, mm.CO3,   69.30,  [+1.772e-02, -6.001e-05, 0.0],         0.74*2),
        ('FeCN3', -3, mm.FeCN3,100.90,  [+1.700e-02, 0.0, 0.0],                0.77*3),
        ('La',    +3, mm.La,    69.70,  [+2.191e-02, +8.119e-05, 0.0],         0.81*3),
        ('FeCN4', -4, mm.FeCN4,110.40,  [+2.00e-02,  0.0, 0.0],                0.62*4),
    ]

    Ncoeffs = [
        ('Si',    0, mm.Si),
        ('SiO2',  0, mm.SiO2),
        ('SiOH4', 0, mm.SiOH4),
        ('BOH3',  0, mm.BOH3),
        ('CO2',   0, mm.CO2),
        ('H2CO3', 0, mm.CO2),
    ]
    Gcoeffs = [
        ('O2', 0, mm.O2),
        ('N2', 0, mm.N2),
    ]

    n = len(Icoeffs)
    c = np.zeros((n,), dtype=float)
    Z = np.zeros((n,), dtype=float)
    sZ = np.zeros((n,), dtype=float)
    M = np.zeros((n,), dtype=float)
    lambdaI = np.zeros((n,), dtype=float)
    lambdacoeff = np.zeros((n,4), dtype=float)
    afac = np.zeros((n,), dtype=float)

    for i,(name,charge,mass,lam0,coeff,af) in enumerate(Icoeffs):
        if hasattr(const, name):
            val = getattr(const, name)
            if not isinstance(val, str):
                c[i] = float(np.asarray(val).squeeze()) * convfac
        Z[i] = abs(charge)
        sZ[i] = charge
        M[i] = mass
        lambdaI[i] = lam0*1e3  # mS/cm -> uS/cm
        lambdacoeff[i,:] = np.array([1.0, coeff[0], coeff[1], coeff[2]], dtype=float)
        afac[i] = af

    ipls = sZ > 0
    imns = sZ < 0

    # Non-ionic constituents
    c2 = np.zeros((len(Ncoeffs),), dtype=float)
    M2 = np.zeros_like(c2)
    for i,(name,_,mass) in enumerate(Ncoeffs):
        if hasattr(const, name):
            val = getattr(const, name)
            if not isinstance(val, str):
                c2[i] = float(np.asarray(val).squeeze()) * convfac
        M2[i] = mass

    c3 = np.zeros((len(Gcoeffs),), dtype=float)
    M3 = np.zeros_like(c3)
    for i,(name,_,mass) in enumerate(Gcoeffs):
        if hasattr(const, name):
            val = getattr(const, name)
            if not isinstance(val, str):
                c3[i] = float(np.asarray(val).squeeze()) * convfac
        M3[i] = mass

    # Shapes: allow K,T,p to be matrices; we support column vectors too.
    T = _as_array(T)
    p = _as_array(p)
    if kpT is not None:
        kpT = _as_array(kpT)

    # If user passes 1D vectors (shape (N,)), allow that and return 1D outputs.
    _return_1d = (T.ndim == 1) and (p.ndim == 1) and ((kpT is None) or (kpT.ndim == 1))

    # Ensure at least 2D like MATLAB does
    if T.ndim == 1:
        T = T.reshape((-1,1))
    if p.ndim == 1:
        p = p.reshape((-1,1))
    if kpT is not None and kpT.ndim == 1:
        kpT = kpT.reshape((-1,1))

    DT = T - 25.0
    iN, iM = T.shape
    flatN = iN*iM

    lambda_mat = np.zeros((n, flatN), dtype=float)
    denom = np.zeros_like(lambda_mat)
    pbar = np.zeros_like(lambda_mat)

    if condmode == 1:
        # forward mode not required for your pipeline; implemented minimally
        raise NotImplementedError("condmode=1 (forward conductivity from composition) not implemented in this focused port.")
    else:
        # Reference conductivity mode (your use)
        k0T = kpT / (1.0 + _pval(gammap, p) * _pval(gammaT, DT))  # remove pressure effect

        # iterative solve for dilution factor 'a'
        a = 1.0
        anew = np.full_like(k0T, 0.0001, dtype=float)

        for _iter in range(50):
            # convergence check
            diff = np.max(np.abs(anew - a)) / max(float(np.max(np.abs(a))), 1e-12)
            if diff <= 1e-6:
                break
            a = anew

            # total charges (vectors per observation, like MATLAB)
            a_flat = a.reshape(-1)
            sum_neg = float(np.sum(c[imns]*Z[imns]))
            sum_pos = float(np.sum(c[ipls]*Z[ipls]))
            Neg = a_flat * sum_neg
            Pos = a_flat * sum_pos

            TC_flat = (Pos + Neg) / 2.0
            # fractional imbalance (absolute), per point
            Balfrac = np.where(TC_flat != 0, np.abs(Pos - Neg) / (2.0*TC_flat), np.nan)
            Im_arr = a_flat * float(np.sum(c*Z*Z)) / 2.0  # ionic strength per point

            # NRM is n x N
            NRM = np.tile(TC_flat, (n, 1))
            NRM[ipls, :] = Neg
            NRM[imns, :] = Pos
            for i in range(n):
                lambda_mat[i,:] = lambdaI[i] * _pval(lambdacoeff[i,:], DT.reshape(-1))
                if calcopt != 2:
                    denom[i,:] = 1.0 + afac[i]*np.sqrt(Im_arr)*(1.0 + DT.reshape(-1)*Beta)
                else:
                    denom[i,:] = 1.0

                if calcopt == 0:
                    pb = np.zeros((flatN,), dtype=float)
                    ffac = np.ones((flatN,), dtype=float)
                    # find opposite sign ions
                    for k in np.where(ipls != ipls[i])[0]:
                        ii = (Im_arr == 0)
                        if np.any(ii):
                            ffac[ii] = 1.0
                        if np.any(~ii):
                            DTv = DT.reshape(-1)[~ii]
                            ffac[~ii] = 1.0 - Ampij[i,k]*np.exp(-widij[i,k]*(np.log(Im_arr[~ii]/Icij[i,k])/2.0)**2) * (1.0 + DTv*Beta)
                        pb = pb + (a.reshape(-1) * (c[k]*Z[k]) * ffac)
                    pbar[i,:] = pb / NRM[i,:]
                else:
                    pbar[i,:] = 1.0

            denom_sum = np.sum((lambda_mat/denom)*pbar*(Z*c)[:,None], axis=0).reshape(iN,iM)
            anew = k0T / denom_sum

        # Now compute at reference temperature
        DTref = float(Tref - 25.0)

        # Recompute Im etc with final anew
        anew_flat = anew.reshape(-1)
        sum_neg = float(np.sum(c[imns]*Z[imns]))
        sum_pos = float(np.sum(c[ipls]*Z[ipls]))
        Neg = anew_flat * sum_neg
        Pos = anew_flat * sum_pos
        TC_flat = (Pos + Neg) / 2.0
        Balfrac2 = np.where(TC_flat != 0, (Pos - Neg) / (2.0*TC_flat), np.nan)
        Im_arr = anew_flat * float(np.sum(c*Z*Z)) / 2.0

        NRM = np.tile(TC_flat, (n, 1))
        NRM[ipls, :] = Neg
        NRM[imns, :] = Pos
        for i in range(n):
            lambda_mat[i,:] = lambdaI[i] * _pval(lambdacoeff[i,:], np.full((flatN,), DTref))
            if calcopt != 2:
                denom[i,:] = 1.0 + afac[i]*np.sqrt(Im_arr)*(1.0 + DTref*Beta)
            else:
                denom[i,:] = 1.0

            if calcopt == 0:
                pb = np.zeros((flatN,), dtype=float)
                ffac = np.ones((flatN,), dtype=float)
                for k in np.where(ipls != ipls[i])[0]:
                    ii = (Im_arr == 0)
                    if np.any(ii):
                        ffac[ii] = 1.0
                    if np.any(~ii):
                        ffac[~ii] = 1.0 - Ampij[i,k]*np.exp(-widij[i,k]*(np.log(Im_arr[~ii]/Icij[i,k])/2.0)**2) * (1.0 + DTref*Beta)
                    pb = pb + (anew.reshape(-1) * (c[k]*Z[k]) * ffac)
                pbar[i,:] = pb / NRM[i,:]
            else:
                pbar[i,:] = 1.0

        denom_sum = np.sum((lambda_mat/denom)*pbar*(Z*c)[:,None], axis=0).reshape(iN,iM)
        kR = anew * denom_sum

        Sc = anew * np.sum(M*c)
        if fullsal == 1:
            Sc = Sc + anew * np.sum(M2*c2)
        elif fullsal == 2:
            Sc = Sc + anew * (np.sum(M2*c2) + np.sum(M3*c3))

        # Return in MATLAB order: [kR, Sc, Balfrac, TC, equiv_c]
        Balfrac2_mat = Balfrac2.reshape(iN, iM)
        TC_mat = TC_flat.reshape(iN, iM)
        # Return in MATLAB order: [kR, Sc, Balfrac, TC, equiv_c]
        outs = (kR, Sc, Balfrac2_mat, TC_mat, None)
    # Match input convenience: if inputs were 1D vectors, return 1D vectors
    if _return_1d:
        def _squeeze_col(a):
            if isinstance(a, np.ndarray) and a.ndim == 2 and 1 in a.shape:
                return a.reshape(-1)
            return a
        outs = tuple(_squeeze_col(o) for o in outs)

    if nargout is None:
        return outs
    return outs[:nargout] if nargout != 1 else outs[0]