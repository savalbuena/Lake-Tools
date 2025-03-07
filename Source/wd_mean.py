# Calculate mean of cardinal direction (wind direction) following circular statistics decribed in Berens, P. (2009). CircStat: A MATLAB toolbox for circular statistics. Journal of Statistical Software, 31(10), 1–21, doi:10.1002/wics.10

# Input arguments:
# wd - 1D vector of wind direction values
# units - 'deg' or 'rad' (default is 'deg')
# nanvar - 'omitnan' to use nanmean averaging; else NaN's in data create
# NaN output (default is to generate NaN output when NaN's exist)

# Use of function
# wd_mean(wd,units,nanvar)
# Returns the scalar variable meanwd, average wind direction in deg

import math
import numpy as np

pi = math.pi


def wd_mean(wd, *args):
    nargin = len(args)
    if nargin <= 1:  # set default units and nanvar approach if unspecified
        units = 'deg'
        nanvar = 'omitnan'
    elif nargin <= 2:  # set default nanvar approach if unspecified
        nanvar = 'omitnan'
        units = args[0]

    if units == 'rad':
        wd = wd * (180 / pi)

    u = np.sin(np.radians(wd))
    v = np.cos(np.radians(wd))

    if nanvar == 'omitnan':
        mean_u = np.nanmean(u)
        mean_v = np.nanmean(v)
    else:
        mean_u = np.mean(u)
        mean_v = np.mean(v)

    mean_r = np.sqrt(mean_u**2 + mean_v**2)
    meanwd = np.arctan2(mean_u, mean_v) * (180 / pi)

    if meanwd < 0:
        meanwd = 360 + meanwd
    elif meanwd == 360:
        meanwd = 0
    elif meanwd > 360:
        meanwd = meanwd - 360

    return meanwd
