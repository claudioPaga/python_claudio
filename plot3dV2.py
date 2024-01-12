#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 08:58:03 2019

@author: cp232
"""
from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import Image
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import astropy
from astropy.io import ascii
import pandas as pd
from scipy.interpolate import griddata


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d


tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/deltaMuLinFitsG19.acs')


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


fig = plt.figure(figsize = [16,12])
ax = Axes3D(fig)
ax.view_init(30, 135)
surf = ax.plot_trisurf(x, y, z, cmap=cm.coolwarm, linewidth=0.1)
ax.set_zlim(1100, 1600)
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_zlabel('Released e-')

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('/Users/cp232/eclipse-workspace/CTILab/outputs/releaseV2.png')
plt.show()


Z = tdata['ctiIndex']
z = Z
xyz = {'x': x, 'y': y, 'z': z}

# put the data into a pandas DataFrame (this is what my data looks like)
df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 
fig = plt.figure(figsize = [16,12])
ax = Axes3D(fig)
ax.view_init(30, 135)
surf = ax.plot_trisurf(x, y, z, cmap=cm.coolwarm, linewidth=0.1)
ax.set_zlim(0,1)
ax.set_xlabel('Row Index')
ax.set_ylabel('Strip Index')
ax.set_zlabel('CTI coefficient, G>19')
plt.title('Total charge release in post-scan trails')

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('/Users/cp232/eclipse-workspace/CTILab/outputs/ctiIndexV2.png')
plt.show()





