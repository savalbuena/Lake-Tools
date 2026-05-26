
"""
Constants ported from limconsts.m (R. Pawlowicz).
"""

from __future__ import annotations
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Dict

def _ns(d: Dict[str, Any]) -> SimpleNamespace:
    return SimpleNamespace(**d)

# Atomic weights (2009 values from the MATLAB file)
mm = _ns({
    'H':   1.00794,
    'Li':  6.941,
    'B':  10.811,
    'C':  12.0107,
    'N':  14.0067,
    'O':  15.9994,
    'F':  18.9984032,
    'Na': 22.98976928,
    'Mg': 24.3050,
    'Si': 28.0855,
    'P':  30.9737,
    'S':  32.065,
    'Cl': 35.453,
    'Ar': 39.948,
    'K':  39.0983,
    'Ca': 40.078,
    'Mn': 54.938045,
    'Fe': 55.845,
    'Cu': 63.546,
    'Zn': 65.409,
    'Br': 79.904,
    'Sr': 87.62,
    'Ag': 107.8682,
    'I':  126.90447,
    'Ba': 137.327,
    'La': 138.90547,
    'As': 74.92160,
    'Al': 26.9815386,
    'Sb': 121.760,
    'Sn': 118.710,
})

# Subsidiary masses for complex molecules
mm.SO4   = mm.S + mm.O*4
mm.HCO3  = mm.H + mm.C + mm.O*3
mm.CO3   = mm.C + mm.O*3
mm.H2CO3 = mm.C + mm.O*2
mm.CO2   = mm.C + mm.O*2
mm.BOH3  = mm.B + (mm.O + mm.H)*3
mm.BOH4  = mm.B + (mm.O + mm.H)*4
mm.OH    = mm.O + mm.H
mm.NO3   = mm.N + mm.O*3
mm.NH4   = mm.N + mm.H*4
mm.SiO2  = mm.Si + mm.O*2
mm.SiOH4 = mm.Si + (mm.O + mm.H)*4
mm.O2    = mm.O*2
mm.N2    = mm.N*2
mm.ClO4  = mm.Cl + mm.O*4
mm.FeCN3 = mm.Fe + 6*(mm.C + mm.N)
mm.FeCN4 = mm.Fe + 6*(mm.C + mm.N)
mm.H2S   = mm.S + mm.H*2
mm.HS    = mm.S + mm.H

# Charges on ions
zval = _ns({
    'H':1,'Li':1,'B':0,'C':0,'N':0,'O':0,'F':-1,'Na':1,'Mg':2,'Si':0,
    'P':0,'S':-2,'Cl':-1,'Ar':0,'K':1,'Ca':2,'Mn':2,'Cu':2,'Zn':2,'Br':-1,'Sr':2,
    'Ag':1,'I':-1,'Ba':2,'La':3,
    'SO4':-2,'HCO3':-1,'CO3':-2,'H2CO3':0,'CO2':0,'BOH3':0,'BOH4':-1,'OH':-1,'NO3':-1,'NH4':1,
    'SiO2':0,'SiOH4':0,'O2':0,'N2':0,'ClO4':-1,'FeCN3':-3,'FeCN4':-4,'H2S':0,'HS':-1
})

# Blank ion structure for natural waters
cstruct = _ns({
    'units':'uM|mM|M',
    'HCO3':0., 'CO3':0., 'H':0., 'OH':0., 'Cl':0., 'NO3':0., 'SO4':0., 'Ca':0., 'Na':0., 'K':0., 'Mg':0.,
    'NH4':0., 'Br':0., 'BOH4':0., 'Sr':0., 'F':0.,
})

# KCL stock solution
KCL = _ns({**cstruct.__dict__})
KCL.units = 'uM'
KCL.K  = 1000.
KCL.Cl = 1000.
KCL.H  = 0.1
KCL.OH = 0.1

def as_const(obj: Any) -> SimpleNamespace:
    """
    Convert a dict or SimpleNamespace into a SimpleNamespace, leaving already-namespace objects unchanged.
    """
    if isinstance(obj, SimpleNamespace):
        return obj
    if isinstance(obj, dict):
        return SimpleNamespace(**obj)
    # Allow simple objects with __dict__
    if hasattr(obj, "__dict__"):
        return SimpleNamespace(**obj.__dict__)
    raise TypeError("const must be a dict-like object or SimpleNamespace")
