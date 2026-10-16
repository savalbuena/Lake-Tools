
# lim (Python port of MATLAB limnological chemistry / conductivity functions)

This folder contains a pure-Python (NumPy/SciPy) port of the MATLAB functions for `lim` created by Rich Palowics. Manuscript available at https://doi.org/10.4319/lom.2008.6.489 

- `limcarbonate`  (carbonate system helper)
- `limcond`       (conductivity/salinity from in-situ conductivity)
- `limstate`      (limnological equation of state / density)

The goal is to preserve **function names** and the overall **call signatures** you used from MATLAB.  
In Python you don't have `nargout`, so these functions accept an optional `nargout=` keyword that
mimics MATLAB output slicing.

## Install / import (local)

```bash
# from the directory containing this README
python -c "import sys; sys.path.insert(0,'.'); import lim; print(lim.limstate(0,10,1000))"
```

Or add the parent folder to `PYTHONPATH`.

## Usage on a simple python script

```python
import numpy as np
from lim import limcarbonate, limcond, limstate

# 1) Ionic composition (molar units)
CQ = {
    'units': 'M',
    'alk': 9.18e-4,
    'pH': 7.3,
    'Cl': 0.57e-4,
    'Si': 3.68e-4,
    'SO4': 0.21e-4,
    'Ca': 2.10e-4,
    'Na': 2.61e-4,
    'K': 0.41e-4,
    'Mg': 1.03e-4
}

T_ref = 25.0

# 2) Carbonate speciation (if you need it)
carb = limcarbonate('TA', CQ['alk'], CQ['pH'], T_ref)

# 3) Create 1D column vectors (N x 1)
cond_mat = cond_avg[:, i].reshape(-1, 1)
T_mat    = T_avg[:, i].reshape(-1, 1)
p_mat    = pressure.reshape(-1, 1)

# 4) limcond (MATLAB-style nargout)
K25, ScLTP = limcond(cond_mat, T_mat, p_mat, T_ref, CQ, nargout=2)

Sc_gkg = ScLTP.flatten() / 1000.0
rho_lim = limstate(Sc_gkg.reshape(-1,1), T_mat, p_mat)
```

## Notes

- `limcond` reads `limPFAC.mat` (included) for the ion-pair interaction fit parameters.
- This port targets your main workflow (`condmode=0`, i.e., reference conductivity from measured conductivity).
