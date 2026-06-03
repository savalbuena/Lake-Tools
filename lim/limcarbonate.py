
from __future__ import annotations
import numpy as np
from types import SimpleNamespace
from typing import Any, Dict, Tuple, Optional

from .limstate import lim_state

def _as_array(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=float)

def _is_str(x: Any) -> bool:
    return isinstance(x, str)

def lim_carbonate(*args, nargout: int | None = None):
    """
    Port of limcarbonate.m (R. Pawlowicz).

    This implementation follows the MATLAB logic closely and supports the main 'TYPE' modes:
      'TA', 'DIC', 'SAT', 'AD', 'HH', and 'K' (constants).

    Units:
      - If TYPE ends with 'KG' (e.g., 'TAKG'): inputs are per kg solution (mol/kg).
      - Otherwise inputs are per liter solution (molar), internally converted to mol/kg using limstate().
    """
    if len(args) == 0:
        raise TypeError("limcarbonate requires arguments")

    typ = 'TA'
    args = list(args)
    if _is_str(args[0]):
        typ = args[0].upper()
        args = args[1:]

    # Default arguments
    T = 25.0
    S = 0.0
    bor = 0.0
    Pdbar = 0.0
    T2 = None
    Pdbar2 = 0.0

    # Special testing: return constants
    if typ == 'K':
        if len(args) < 2:
            raise TypeError("limcarbonate('K', T, S, [P]) requires T and S")
        T = args[0]
        S = args[1]
        Pdbar = args[2] if len(args) >= 3 else 0.0
        Kw,K1,K2,K0,KB,KSi,KSO4,KHF,KHS,pHT2pHF,pHT2pHS = Kconst(T, S, Pdbar)
        water_type = 'Seawater' if np.any(np.asarray(S) < 0) else 'freshwater'
        Carb = {
            'S': np.abs(_as_array(S)),
            'T': _as_array(T),
            'type': water_type,
            'Kw': Kw, 'pKw': -np.log10(Kw),
            'K0': K0, 'pK0': -np.log10(K0),
            'K1': K1, 'pK1': -np.log10(K1),
            'K2': K2, 'pK2': -np.log10(K2),
            'KB': KB, 'pKB': -np.log10(KB),
            'KSi': KSi, 'pKSi': -np.log10(KSi),
            'KSO4': KSO4, 'pKSO4': -np.log10(KSO4),
            'KHF': KHF, 'pKHF': -np.log10(KHF),
            'KHS': KHS, 'pKHS': -np.log10(KHS),
            'pHT2pHF': pHT2pHF,
            'pHT2pHS': pHT2pHS,
        }
        return _maybe_nargout(Carb, nargout)

    # Parse optional T,S,P,T2,P2 (after first 2 required args for each TYPE)
    if len(args) >= 3:
        T = args[2]
    if len(args) >= 4:
        S = args[3]
        bor = 0.0
        if isinstance(S, dict) or isinstance(S, SimpleNamespace):
            Sd = S.__dict__ if isinstance(S, SimpleNamespace) else S
            bor = Sd.get('B', 0.0)
            S = Sd.get('S', 0.0)
    if len(args) >= 5:
        Pdbar = args[4]
    if len(args) >= 6:
        T2 = args[5]
    if len(args) >= 7:
        Pdbar2 = args[6]

    # Broadcast sizing to match MATLAB behavior
    X = _as_array(args[0])
    Y = _as_array(args[1]) if len(args) > 1 else None
    Sarr = _as_array(S)
    Tarr = _as_array(T)
    nsz = np.broadcast(X, X, Sarr, Tarr).shape
    if np.size(Sarr) == 1 and np.prod(nsz) > 1:
        Sarr = np.full(nsz, float(Sarr))
    if np.size(Tarr) == 1 and np.prod(nsz) > 1:
        Tarr = np.full(nsz, float(Tarr))

    # Units logic
    if typ.endswith('KG'):
        rho = np.ones_like(Sarr, dtype=float)
        units = 'mol/kg'
    else:
        rho = lim_state(Sarr, Tarr, Pdbar)
        units = 'M'

    # Chemistry constants at (T,S,P)
    Kw,K1,K2,K0,KB,KSi,KSO4,KHF,KHS,pHT2pHF,pHT2pHS = Kconst(Tarr, Sarr, Pdbar)

    # Default borate
    if np.all(np.asarray(bor) == 0):
        bor_est = -Sarr/1.80655*(0.0052004+0.0162211)*1e-3
        bor_est = np.where(Sarr >= 0, 0.0, bor_est)
        bor_use = bor_est
    else:
        bor_use = _as_array(bor)

    # Main switch
    typ_base = typ  # keep original
    if typ_base.endswith('KG'):
        typ_base = typ_base[:-2]  # drop 'KG' suffix for dispatch

    Carb = {}

    if typ_base == 'TA':
        TA = X / rho
        pHt = _as_array(args[1]) * np.ones_like(rho)

        H = 10.0**(-pHt)
        OH = Kw / H
        B4 = bor_use*KB/(H+KB)
        B3 = bor_use - B4

        HCO3 = (TA - B4 - OH + H) / (1.0 + 2.0*K2/H)
        CO3 = HCO3 * K2 / H
        H2CO3 = HCO3 * H / K1

    elif typ_base == 'DIC':
        DIC = X / rho
        pHt = _as_array(args[1]) * np.ones_like(rho)

        H = 10.0**(-pHt)
        OH = Kw / H
        B4 = bor_use*KB/(H+KB)
        B3 = bor_use - B4

        HCO3 = DIC / (H/K1 + 1.0 + K2/H)
        CO3 = HCO3 * K2 / H
        H2CO3 = HCO3 * H / K1

    elif typ_base == 'SAT':
        TA = X / rho
        fug = _as_array(args[1])
        H2CO3 = K0 * fug

        H = np.full_like(K1, np.nan, dtype=float)
        # Solve polynomial per element (as MATLAB does)
        TAf = np.ravel(TA); K1f=np.ravel(K1); K2f=np.ravel(K2); Kwf=np.ravel(Kw); KBf=np.ravel(KB); borf=np.ravel(bor_use); H2f=np.ravel(H2CO3)
        outH = np.full_like(TAf, np.nan, dtype=float)
        for i in range(outH.size):
            if np.isnan([TAf[i],K1f[i],K2f[i],Kwf[i],H2f[i]]).any():
                continue
            if borf[i] == 0:
                C = np.array([1.0, TAf[i], -(K1f[i]*H2f[i] + Kwf[i]), -2.0*K1f[i]*K2f[i]*H2f[i]])
            else:
                C = np.array([
                    1.0,
                    TAf[i] + KBf[i],
                    (KBf[i]*(TAf[i]-borf[i]) - (K1f[i]*H2f[i] + Kwf[i])),
                    -(KBf[i]*(K1f[i]*H2f[i] + Kwf[i]) + 2.0*K1f[i]*K2f[i]*H2f[i]),
                    -KBf[i]*2.0*K1f[i]*K2f[i]*H2f[i]
                ])
            roots = np.roots(C)
            real_roots = roots[np.isclose(roots.imag, 0)]
            if real_roots.size == 0:
                continue
            val = np.max(real_roots.real)
            outH[i] = val if val > 0 else np.nan
        H = outH.reshape(TA.shape)
        pHt = -np.log10(H)

        HCO3 = K1*H2CO3/H
        CO3 = K2*HCO3/H
        OH = Kw/H
        B4 = bor_use*KB/(H+KB)
        B3 = bor_use - B4

    elif typ_base == 'AD':
        TA = X / rho
        DIC = _as_array(args[1]) / rho

        TAf = np.ravel(TA); DICf=np.ravel(DIC); K1f=np.ravel(K1); K2f=np.ravel(K2); Kwf=np.ravel(Kw); KBf=np.ravel(KB); borf=np.ravel(bor_use)
        outH = np.full_like(TAf, np.nan, dtype=float)
        for i in range(outH.size):
            if np.isnan([TAf[i],DICf[i],K1f[i],K2f[i],Kwf[i]]).any():
                continue
            if borf[i] == 0:
                C = np.array([
                    1.0/K1f[i],
                    (TAf[i]/K1f[i] + 1.0),
                    (-Kwf[i]/K1f[i] + TAf[i] + K2f[i] - DICf[i]),
                    (-Kwf[i] + K2f[i]*TAf[i] - 2.0*K2f[i]*DICf[i]),
                    -K2f[i]*Kwf[i]
                ])
            else:
                C = np.array([
                    -1.0,
                    -TAf[i]-KBf[i]-K1f[i],
                    DICf[i]*K1f[i]-TAf[i]*(KBf[i]+K1f[i])+KBf[i]*borf[i]+Kwf[i]-KBf[i]*K1f[i]-K1f[i]*K2f[i],
                    (DICf[i]*(KBf[i]*K1f[i]+2*K1f[i]*K2f[i]) - TAf[i]*(KBf[i]*K1f[i]+K1f[i]*K2f[i]) + KBf[i]*borf[i]*K1f[i]
                     + (Kwf[i]*KBf[i]+Kwf[i]*K1f[i]-KBf[i]*K1f[i]*K2f[i])),
                    (2.*DICf[i]*KBf[i]*K1f[i]*K2f[i] - TAf[i]*KBf[i]*K1f[i]*K2f[i] + KBf[i]*borf[i]*K1f[i]*K2f[i]
                     + (Kwf[i]*KBf[i]*K1f[i] + Kwf[i]*K1f[i]*K2f[i])),
                    Kwf[i]*KBf[i]*K1f[i]*K2f[i]
                ])
            roots = np.roots(C)
            real_roots = roots[np.isclose(roots.imag, 0)]
            if real_roots.size == 0:
                continue
            val = np.max(real_roots.real)
            outH[i] = val if val > 0 else np.nan
        H = outH.reshape(TA.shape)
        pHt = -np.log10(H)

        OH = Kw/H
        H2CO3 = DIC / (1.0 + K1/H + K1*K2/(H*H))
        HCO3 = H2CO3*K1/H
        CO3 = HCO3*K2/H
        B4 = bor_use*KB/(H+KB)
        B3 = bor_use - B4

    elif typ_base == 'HH':
        HCO3 = X / rho
        CO3 = _as_array(args[1]) / rho
        H = K2*HCO3/CO3
        pHt = -np.log10(H)
        OH = Kw/H
        H2CO3 = HCO3*H/K1
        B4 = bor_use*KB/(H+KB)
        B3 = bor_use - B4

    else:
        raise ValueError(f"Unsupported TYPE '{typ}' in this Python port.")

    # Closed-system move to (T2,P2) if provided
    if T2 is not None:
        DIC = H2CO3 + HCO3 + CO3
        Kw2,K12,K22,K02,KB2,KSi2,_,_,_,pHT2pHF2,pHT2pHS2 = Kconst(_as_array(T2), Sarr, Pdbar2)
        # Solve at new T2 keeping DIC and TA fixed using 'AD' approach (TA and DIC)
        # We implement the same polynomial method by re-calling limcarbonate('AD', TA*rho, DIC*rho, T2, S, P2)
        # Convert back to original units convention:
        TA_input = ( (HCO3 + 2*CO3 + (bor_use*KB/(H+KB)) + OH - H) )  # mol/kg
        # TA_input and DIC are mol/kg; scale by rho if original units were mol/L.
        if not typ.endswith('KG'):
            TA_call = TA_input * rho
            DIC_call = DIC * rho
        else:
            TA_call = TA_input
            DIC_call = DIC
        Carb2 = lim_carbonate('AD' + ('KG' if typ.endswith('KG') else ''), TA_call, DIC_call, T2, Sarr, Pdbar2)
        return _maybe_nargout(Carb2, nargout)

    # Assemble output structure (dict)
    # MATLAB stores concentrations in mol/l solution if inputs were molar.
    # Here we keep the same convention.
    Carb.update({
        'type': ('Seawater' if np.any(Sarr < 0) else 'freshwater'),
        'units': units,
        'T': Tarr,
        'S': np.abs(Sarr),
        'pH_total': pHt,
        'pH_free': pHt + pHT2pHF,
        'pH_SWS': pHt + pHT2pHS,
        'H2CO3': H2CO3 * rho,   # back to per L if molar inputs
        'HCO3':  HCO3  * rho,
        'CO3':   CO3   * rho,
        'H':     (10.0**(-(pHt + pHT2pHF))) * rho,  # 'free' scale as in MATLAB output line
        'OH':    OH * rho,
        'BOH4':  (bor_use*KB/(10.0**(-pHt)+KB)) * rho,
        'BOH3':  (bor_use - bor_use*KB/(10.0**(-pHt)+KB)) * rho,
        'DIC':   (H2CO3 + HCO3 + CO3) * rho,
        'TA':    (HCO3 + 2*CO3 + (bor_use*KB/(10.0**(-pHt)+KB)) + OH - 10.0**(-pHt)) * rho,
        'fCO2':  (H2CO3 / K0),  # fugacity that would match
    })

    return _maybe_nargout(Carb, nargout)

def _maybe_nargout(obj, nargout: int | None):
    if nargout is None:
        return obj
    if isinstance(obj, tuple):
        return obj[:nargout]
    # for dict / scalar, MATLAB would return one output only
    return obj

def Kconst(T, Sin, Pdbar=0.0):
    """
    Port of nested Kconst() in limcarbonate.m.
    Returns (Kw,K1,K2,K0,KB,KSi,KSO4,KHF,KHS,pHT2pHF,pHT2pHS)
    """
    T = _as_array(T)
    Sin = _as_array(Sin)
    Pdbar = _as_array(Pdbar)

    T, Sin, Pdbar = np.broadcast_arrays(T, Sin, Pdbar)

    TK = T + 273.15
    lTK = np.log(TK)
    T100 = TK/100.0

    # Determine case
    finite = np.isfinite(Sin)
    if np.all((Sin[finite] >= 0)):
        limcase = True
        S = Sin/1.00488
    elif np.all((Sin[finite] <= 0)):
        limcase = False
        S = -Sin
    else:
        raise ValueError("Kconst: cannot mix S>=0 and S<=0 in same call (ignoring NaNs).")

    S05 = np.sqrt(S)
    S15 = S*S05
    S2 = S*S

    # K0 Weiss 1974 (mol/kg/atm)
    K0 = np.exp(-60.2409 + 9345.17/TK + 23.3585*np.log(T100)
                + (0.023517 - 0.023656*T100 + 0.0047036*T100**2)*S)

    # Kw (mol/kg solution)
    Kw = np.exp(148.9652 - 13847.26/TK - 23.6521*lTK
                + (-5.977 + 118.67/TK + 1.0495*lTK)*S05
                - 0.01615*S)

    # Boric acid (total scale)
    KB = np.exp(( -8966.90 - 2890.53*S05 - 77.942*S + 1.728*S15 - 0.0996*S2)/TK
                + (148.0248 + 137.1942*S05 + 1.62142*S)
                + (-24.4344 - 25.085*S05 - 0.2474*S)*lTK
                + 0.053105*S05*TK)

    # Ionic strength
    I = 19.924*S/(1000 - 1.005*S)
    sI = np.sqrt(I)

    # Silicic acid
    KSi = np.exp(117.385 - 8904.2/TK - 19.334*lTK
                 + (-458.79/TK + 3.5913)*sI
                 + (188.74/TK - 1.5998)*I
                 + (-12.1652/TK + 0.07871)*I**2
                 + np.log(1 - 0.001005*S))

    # Sulfate (free scale)
    KSO4 = np.exp(-4276.1/TK + 141.328 - 23.093*lTK
                  + (-13856/TK + 324.57 - 47.986*lTK)*sI
                  + (35474/TK - 771.54 + 114.723*lTK)*I
                  + (-2698/TK)*I*sI
                  + (1776/TK)*I**2
                  + np.log(1 - 0.001005*S))

    # HF (free scale)
    KHF = np.exp(1590.2/TK - 12.641 + 1.525*sI + np.log(1 - 0.001005*S))

    # K1/K2
    if limcase:
        K1 = np.exp(290.9097 - 14554.21/TK - 45.0575*lTK
                    + (-228.39774 + 9714.36839/TK + 34.485796*lTK)*S05
                    + (54.20871 - 2310.48919/TK - 8.19515*lTK)*S
                    + (-3.969101 + 170.22169/TK + 0.603627*lTK)*S15
                    - 0.00258768*S2 + np.log(1 - 0.001005*S))
        K2 = np.exp(207.6548 - 11843.79/TK - 33.6485*lTK
                    + (-167.69908 + 6551.35253/TK + 25.928788*lTK)*S05
                    + (39.75854 - 1566.13883/TK - 6.171951*lTK)*S
                    + (-2.892532 + 116.270079/TK + 0.45788501*lTK)*S15
                    - 0.00613142*S2 + np.log(1 - 0.001005*S))
    else:
        K1 = 10.0**-(-61.2172 + 3633.86/TK + 9.67770*lTK + S*(-0.011555 + S*(0.0001152)))
        K2 = 10.0**-(25.9290 + 471.78/TK - 3.16967*lTK + S*(-0.01781 + S*(0.0001122)))

    # H2S
    KHS = np.exp(225.838 - 13275.3/TK - 34.6435*lTK + (0.3449)*S05 + (-0.0274)*S)

    # Pressure corrections (Millero 1983 approach, as in MATLAB)
    P = Pdbar/10.0
    R = 83.1451
    T2 = T*T
    nS05 = np.sqrt(S/35.0)

    # Helper to apply exp factor
    def _apply(K, dV0, dk0, dVs, dks):
        dV = dV0 + nS05*(dVs - dV0)
        dK = dk0 + nS05*(dks - dk0)
        return K * np.exp(-dV/(R*TK)*P + dK/(2*R*TK)*P**2)

    # Water
    Kw = _apply(Kw,
                -25.60 + 0.2324*T - 3.6246e-3*T2,
                (-7.33 + 0.1368*T - 1.233e-3*T2)/1000.0,
                -20.02 + 0.1119*T - 1.409e-3*T2,
                (-5.13 + 0.0794*T)/1000.0)

    # Boric acid
    KB = _apply(KB,
                -38.75 + 0.1991*T - 2.813e-3*T2,
                (-10.13 + 0.2396*T - 3.961e-3*T2)/1000.0,
                -28.56 + 0.1211*T - 0.321e-3*T2,
                (-3.0 + 0.0427*T)/1000.0)

    # Silicic acid (same as borate)
    KSi = KSi * np.exp((np.log(KB) - np.log(KB)))  # no-op placeholder; actual MATLAB applies same factor
    KSi = _apply(KSi,
                 -38.75 + 0.1991*T - 2.813e-3*T2,
                 (-10.13 + 0.2396*T - 3.961e-3*T2)/1000.0,
                 -28.56 + 0.1211*T - 0.321e-3*T2,
                 (-3.0 + 0.0427*T)/1000.0)

    # Bisulfate
    KSO4 = _apply(KSO4,
                  -23.10 + 0.0990*T - 1.560e-3*T2,
                  (-5.39 + 0.093*T - 1.896e-3*T2)/1000.0,
                  -18.03 + 0.0466*T + 0.316e-3*T2,
                  (-4.53 + 0.0900*T)/1000.0)

    # HF
    KHF = _apply(KHF,
                 -12.81 - 0.0028*T - 1.254e-3*T2,
                 (-5.57 + 0.1310*T - 1.160e-3*T2)/1000.0,
                 -9.78 - 0.0090*T - 0.942e-3*T2,
                 (-3.91 + 0.054*T)/1000.0)

    # K1
    K1 = _apply(K1,
                -30.54 + 0.1849*T - 2.3366e-3*T2,
                (-6.22 + 0.1368*T - 1.233e-3*T2)/1000.0,
                -25.59 + 0.1976*T - 3.111e-3*T2,
                (-3.90 + 0.1700*T)/1000.0)

    # K2
    K2 = _apply(K2,
                -29.81 + 0.1150*T - 1.816e-3*T2,
                (-5.74 + 0.093*T - 1.896e-3*T2)/1000.0,
                -17.32 + 0.1758*T - 2.647e-3*T2,
                (-4.87 + 0.0900*T)/1000.0)

    # H2S (apply analogous correction to KHS)
    KHS = _apply(KHS,
                 -14.80 + 0.0020*T - 0.400e-3*T2,
                 (-(5.57-1.02) + 0.1310*T - 1.160e-3*T2)/1000.0,  # matches MATLAB "added 1.02" comment intent
                 -11.07 - 0.009*T - 0.942e-3*T2,
                 (-2.89 + 0.054*T)/1000.0)

    # pH scale conversions
    if limcase:
        pHT2pHF = np.zeros_like(S)
        pHT2pHS = np.zeros_like(S)
    else:
        SO4 = 1.4573830e-3*S/1.80655
        pHT2pHF = np.log10(1.0 + SO4/KSO4)
        F = 0.0035266e-3*S/1.80655
        pHT2pHS = pHT2pHF - np.log10(1.0 + SO4/KSO4 + F/KHF)

    return Kw,K1,K2,K0,KB,KSi,KSO4,KHF,KHS,pHT2pHF,pHT2pHS
