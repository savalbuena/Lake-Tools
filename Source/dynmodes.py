# dynmodes.py calculates ocean dynamic vertical modes

# Usage dynmodes(Nsq, p, nmodes)
# taking a column vector of Brunt-Vaisala values (Nsq) at different pressures (p) and calculating some number of dynamic modes (nmodes).
# Note: The input pressures need not be uniformly spaced, and the deepest pressure is assumed to be the bottom.

# The program is based on the developement in section 6.13 of Atmosphere-Ocean Dynamics by Adrian Gill, Academic Press, 1982.
# The variables describing flow in the ocean or atmosphere are separated into vertical (depending only on height or depth), horizontal, and time. For example, the vertical velocity (w) is represented as a sum of different vertical structures as

# $$ w = \sum_{n=0}^\infty \hat h_n(z) \hat w_n (x,y,t). $$

# There is a related set of vertical functions ($\hat p_n(z)$) which represent the vertical structure of horizontal velocity and pressure anomaly. In the program below, wmodes is a matrix containing the $\hat h_n(z)$ functions and pmodes are the $\hat p_n(z)$ functions. wmodes(:,n) or pmodes(:,n) is the nth vertical structure at the depths (p) provided as input.

# Each set of functions is mutually orthogonal with appropriate weighting. The governing equation for the vertical structures (6.11.18 in Gill) is
#  $$ d^2/dz^2 \hat h_n + (N^2(z) /c_e) \hat h_n = 0, n = 1, 2, ... , $$

# where $N^2(z)$ is the Brunt-Vaiasala freqency squared and $c_e$ is the eigenvalue, which is the phase speed of the internal gravity wave with the indicated vertical structure. The boundary conditions are that $\hat h$ vanishes at the surface and bottom (the rigid lid assumption). The solution for n=0 is a special case described in Gill which is very weakly affected by the density structure.

# The vertical velocity structure (wmodes) is calculated at the depths (p) for which $N^2$ is provided on input. The differential equation above is converted to a set of algebraic equations using second order central finite differences. The linear algebraic equation is $A * wmodes + \lambda B * wmodes = 0$, where B contains $N^2$ along its diagonal. The routine eig(A,B) calculate the eigenvalues $\lambda$ which are the reciprocals of the wave speed $c_e$.

# USAGE: [wmodes,pmodes,ce]=dynmodes(Nsq,p,nmodes);

# Inputs: Nsq = column vector of Brunt-Vaisala buoyancy frequency (s^-2)
#         p = column vector of pressure (decibars)
#         nmodes = number of vertical modes to calculate

# Outputs: wmodes = vertical velocity structure
#          pmodes = horizontal velocity structure
#          ce = modal speed (m/s)

#  developed by John Klinck. July, 1999
# Edited by Sergio Valbuena July 2020.

# NOTE: the wmodes result has not been validated yet.

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

# --------------- Function to create the bathymetry file ----------------------------

def dynmodes(*args):
    nargin = len(args)

    if nargin < 1:
        nplot = 3
        p = np.arange(0,410,10)
        z = -p
        n = len(p)
        Nsq = np.zeros(n)
        Nsq[0:n] = 1e-6
        nmodes = 3
    else:
        Nsq = args[0]
        p = args[1]
        nmodes = args[2]
        z = np.zeros(len(p))


    rho0 = 1028
    if len(Nsq) != len(p):
        print('Error, the vectors Nsq and p do not have the same length')
        exit()
    else:
        n = len(Nsq)

    if p[0] > 0:
        # Add surface pressure with top Nsq value
        z = np.zeros(len(Nsq)+1)
        z[1:n+1] = -p[0:n]
        N2 = np.zeros(len(Nsq)+1)
        N2[0] = Nsq[0]
        N2[1:n+1] = Nsq[0:n]
        nz = n + 1
    else:
        z = -p
        N2 = Nsq
        nz = n
    # Calculte depths and spacing
    dz = np.zeros(nz-1)
    dz[:] = np.nan
    zm = np.zeros(nz-1)
    zm[:] = np.nan
    dz[0:nz-1] = z[0:nz-1] - z[1:nz]
    # midpoint depth
    zm[0:nz-1] = z[0:nz-1] - 0.5*dz[0:nz-1]
    # midpoint spacing
    dzm = np.zeros(nz)
    dzm[1:nz-1] = zm[0:nz-2] - zm[1:nz-1]
    dzm[0] = dzm[2]
    dzm[-1] = dzm[-2]

    # Get dynamic modes
    A = np.zeros((nz,nz))
    B = np.zeros((nz,nz))
    for i in range(1,nz-1):
        A[i,i] = 1/(dz[i-1]*dzm[i]) + 1/(dz[i]*dzm[i])
        A[i,i-1] = -1/(dz[i-1]*dzm[i])
        A[i,i+1] = -1/(dz[i]*dzm[i])

    for i in range(0,nz):
        B[i,i] = N2[i]

    # Setting boundary conditios
    A[0,0] = -1
    A[nz-1,0] = -1
    # Obtaining the eigenvalues and eigenvectors
    [e, wmodes] = la.eig(A,B)
    idxreal = np.isreal(e)
    e = e[idxreal]
    e = e.real
    wmodes = wmodes[:,idxreal]
    ind = e >= 1e-10
    e = e[ind]
    wmodes = wmodes[:,ind]
    e1 = np.sort(e)
    ind = np.argsort(e)
    wmodes[:,:] = wmodes[:,ind]

    nm = len(e)
    ce = 1/np.sqrt(e1)

    pmodes = np.zeros(np.shape(wmodes))
    for i in range(0,nm):
         pr = np.diff(wmodes[:,i])
         pr[0:nz-2] = pr[0:nz-2]/dz[0:nz-2]
         pr = pr*rho0*ce[i]*ce[i]
         pmodes[1:nz-1,i] = 0.5*(pr[1:nz-1]+pr[0:nz-2])
         pmodes[0,i] = pr[0]
         pmodes[nz-1,i] = pr[nz-2]

    return wmodes, pmodes, ce
