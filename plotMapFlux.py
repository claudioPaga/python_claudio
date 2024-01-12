#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 14:08:30 2021

@author: cp232
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 10:45:30 2019

Inputs - Flux maps from WM processing.

Outputs - 2D maps plots + Median cumulative flux along AC, to check for CTI signatured   

Summary - compares maps of damaged WMRs generated using different trap species/bg values, plotting the flux difference.

Updates - Proc updated many times as new flux maps were generated after processing WMRs with new CDM libraries with 3rd trap species added.
          Last update on 19/12/2019, after processing WRMs with a 3rd trap species with density of 0.5 traps/pixels.

@author: cp232
"""

import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

#
#------------------------------
# Damaged flux, BG = 0
#------------------------------
#

# Read in array of output flux, columns are x, y, flux
tdata = ascii.read('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux5fr3T-100/node_flux_frn0.txt')

X = tdata['Xbinned']
Y = tdata['Ybinned']
Z = tdata['Flux']
# Remove NANs
filter = np.logical_not(np.isnan(Z))
x = X
y = Y
z = Z
x = x[filter]
y = y[filter]
z = z[filter]

xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
print('Total flux = ', np.nanmean(z2))
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
meanAC = np.zeros(z2.shape[1])
for i in range(z2.shape[1]):
    meanAC[i] = z2[:,i].mean()

meanAL = np.zeros(z2.shape[0])
for i in range(z2.shape[0]):
    meanAL[i] = z2[i,:].mean()
    