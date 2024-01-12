#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 08:33:51 2019

@author: cp232
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import astropy
from astropy.io import ascii
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/deltaMuLinFitsG19.acs')

X = tdata['row']
Y = tdata['strip']
#Z = tdata['totElectronsReleasedCI50K']
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)
fig.colorbar(surf, shrink=2.0, aspect=5)

plt.title('Original Code')
# ~~~~ MODIFICATION TO EXAMPLE BEGINS HERE ~~~~ #
import pandas as pd
from scipy.interpolate import griddata
# create 1D-arrays from the 2D-arrays
X = tdata['row']
Y = tdata['strip']
Z = tdata['totElectronsReleasedCI50K']
x = X
y = Y
z = Z
#x = X.reshape(1600)
#y = Y.reshape(1600)
#z = Z.reshape(1600)
xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 

# re-create the 2D-arrays
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

fig = plt.figure(figsize = [16,12])
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x2, y2, z2, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False )
ax.set_zlim(1100, 1600)
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_zlabel('Released e-')

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title('Total charge release in post-scan trails')
# ~~~~ MODIFICATION TO EXAMPLE ENDS HERE ~~~~ #
plt.savefig('/Users/cp232/eclipse-workspace/CTILab/outputs/release.png')
plt.show()
