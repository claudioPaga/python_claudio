#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:35:04 2019

@author: cp232
"""

import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
import pandas as pd
import numpy as np
from scipy.interpolate import griddata


tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/deltaMuLinFitsG19.acs')
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/FpaCdmParametersObmt3260.asc')

X = tdata['row']
Y = tdata['strip']
# Z = tdata['totElectronsReleasedCI50K']
Z = tdata['chargeVolumeCoeff']

x = X
y = Y
z = Z

xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

w = 8
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_title('Total charge release in post-scan trails')

color_map = plt.imshow(z2, origin = "lower", extent = [0.5, 7.5, 0.5, 9.5])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/fpa_map_release.png")


Z = tdata['ctiIndex']
z = Z
xyz = {'x': x, 'y': y, 'z': z}


# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

w = 8
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)

ax = fig.gca()
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_title('CTI coefficient, G>19')
#ax.get_clim()
color_map = plt.imshow(z2, origin = "lower", extent = [0.5, 7.5, 0.5, 9.5])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/fpa_map_cti_cfs.png")

     
# Plot of median BG distribution over focal plane

Z = tdata['bgMedian'] # The BG median value for each device
z = Z
xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

w = 8
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)

ax = fig.gca()
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_title('Median CCD Background')
#ax.get_clim()
color_map = plt.imshow(z2, origin = "lower", extent = [0.5, 7.5, 0.5, 9.5])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/fpa_map_bg_cfs.png")

