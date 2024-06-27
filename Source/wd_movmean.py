# Calculate moving mean of cardinal direction (wind direction) following circular
# statistics decribed in Berens, P. (2009). CircStat: A MATLAB toolbox for
# circular statistics. Journal of Statistical Software, 31(10), 1–21, doi:10.1002/wics.10

# Input arguments:
# wd - 1D vector of wind direction values
# k - if one value, is index size of centered moving average function
#       if two values, gives iBack and iForward described below:
#       iBack - index steps back from value i for moving mean (can be 0)
#       iForward - index steps forward from value i for moving mean (can be 0)
# units - 'deg' or 'rad' (default is 'deg')
# nanvar - 'omitnan' to use nanmean averaging; else NaN's in data create
# NaN output (default is to generate NaN output when NaN's exist)

#Returns the 1D vector variable meanwd, average moving wind direction

import numpy as np
import pandas as pd
import math
pi = math.pi

def wd_movmean(wd,k,units,nanvar):
    if units == 'rad':
        wd = wd*(180/pi)

    u = np.sin(np.radians(wd))
    v = np.cos(np.radians(wd))

    if nanvar == 'omitnan':
        mean_u = pd.Series(u).rolling(window=k, min_periods=1).mean().values
        mean_v = pd.Series(v).rolling(window=k, min_periods=1).mean().values
    else:
        mean_u = pd.Series(u).rolling(window=k).mean().values
        mean_v = pd.Series(v).rolling(window=k).mean().values

    mean_r = np.sqrt(mean_u**2 + mean_v**2)
    meanwd = np.arctan2(mean_u,mean_v)*(180/pi)
    meanwd[meanwd<0] = 360 + meanwd[meanwd<0]
    meanwd[meanwd==360] = 0
    meanwd[meanwd>360] = meanwd[meanwd>360] - 360

    return meanwd
