
from __future__ import annotations
import numpy as np
from typing import Any

def _as_array(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _as_vec(x: Any) -> np.ndarray:
    """Convert scalars/arrays to float ndarray, squeezing 2D column/row vectors to 1D."""
    a = _as_array(x)
    if a.ndim == 2 and 1 in a.shape:
        return a.reshape(-1)
    return a


def lim_state(S, T, p=0, *, nargout: int | None = None):
    """
    Port of limstate.m (R. Pawlowicz).

    Returns density in kg/L solution.

    - If S >= 0 (ignoring NaNs): limnological formulation (Chen & Millero 1986)
    - If S <= 0 (ignoring NaNs): seawater formulation (Millero et al. 1980)
    """
    S = _as_vec(S)
    T = _as_vec(T)
    p = _as_vec(p)
    S, T, p = np.broadcast_arrays(S, T, p)

    P = p / 10.0  # dbar -> bar

    finite = np.isfinite(S)
    if np.all(S[finite] >= 0):
        # Limnological equation
        rho0 = (
            0.9998395
            + T*(6.7914e-5
            + T*(-9.0894e-6
            + T*(1.0171e-7
            + T*(-1.2846e-9
            + T*(1.1592e-11
            + T*(-5.0125e-14))))))
            + S*(8.181e-4 + T*(-3.85e-6 + T*(4.96e-8)))
        )
        K = (
            19652.17
            + T*(148.113 + T*(-2.293 + T*(1.256e-2 + T*(-4.18e-5))))
            + P*(3.2726 + T*(-2.147e-4 + T*(1.128e-4)))
            + S*(53.238 + T*(-0.313) + P*(5.728e-3))
        )
        rho = rho0 / (1.0 - P / K)

    elif np.all(S[finite] <= 0):
        # Seawater equation, input is negative PSS-78
        Spos = -S
        S15 = Spos**1.5
        S2 = Spos**2

        rhow = (
            999.842594
            + T*(6.793952e-2
            + T*(-9.095290e-3
            + T*(1.001685e-4
            + T*(-1.120083e-6
            + T*(6.536332e-9)))))
        )

        rho0 = (
            rhow
            + Spos*(8.24493e-1
            + T*(-4.0899e-3
            + T*(7.6438e-5
            + T*(-8.2467e-7
            + T*(5.3875e-9)))))
            + S15*(-5.72466e-3 + T*(1.0227e-4 + T*(-1.6546e-6)))
            + S2*(4.8314e-4)
        )

        Kw = (
            19652.21
            + T*(148.4206
            + T*(-2.327105
            + T*(1.360477e-2
            + T*(-5.155288e-5))))
        )

        Kw0 = (
            Kw
            + Spos*(54.6746 + T*(-0.603459 + T*(1.09987e-2 + T*(-6.1670e-5))))
            + S15*(7.944e-2 + T*(1.6483e-2 + T*(-5.3009e-4)))
        )

        K = (
            Kw0
            + P*(3.239908
            + T*(1.43713e-3 + T*(1.16092e-4 + T*(-5.77905e-7)))
            + Spos*(2.2838e-3 + T*(-1.0981e-5 + T*(-1.6078e-6)))
            + S15*(1.91075e-4))
            + P*P*(8.50935e-5
            + T*(-6.12293e-6 + T*(5.2787e-8))
            + Spos*(-9.9348e-7 + T*(2.0816e-8 + T*(9.1697e-10))))
        )

        rho = rho0 / (1.0 - P / K)
        rho = rho / 1000.0  # kg/m^3 -> kg/L

    else:
        raise ValueError("limstate: S must be all >=0 or all <=0 (ignoring NaNs), matching MATLAB behavior.")

    if nargout is None:
        return rho
    outs = (rho,)
    return outs[:nargout] if nargout != 1 else rho
