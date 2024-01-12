#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:35:10 2019

Inputs - Residual maps from WM processing.
    1 - Obs - Model residuals
    2 - Obs - CDM residuals

Outputs - 2D maps plots    

@author: cp232


"""

# 

import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

# ******************************
# Plot1 - Obs-Model residuals

# Read in array of residuals 
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsModelResG14_15.acs')


X = tdata['alIdx']
Y = tdata['acIdx']
Z = tdata['flux']
z = Z
# Remove NANs
filter = np.logical_not(np.isnan(z))
z = z[filter]
x = X
y = Y
x = x[filter]
y = y[filter]

xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

# Set width and height of image (it might be in cm) and the resolution
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Residuals (obs-model)')

color_map = plt.imshow(z2, origin = "lower", extent = [-9, 9, -6, 6])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_modelG14_15.png")

#
# ******************************
# Plot2 - Obs-Model residuals
#

# Read in array of residuals 
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsDamagedResG14_15.acs')


X = tdata['alIdx']
Y = tdata['acIdx']
Z = tdata['flux']
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

# Set width and height of image (it might be in cm) and the resolution
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Residuals (obs-cdm)')

color_map = plt.imshow(z2, origin = "lower", extent = [-9, 9, -6, 6])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_cdmG14_15.png")



