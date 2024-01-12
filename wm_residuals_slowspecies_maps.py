#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10th 11:35:10 2019

Inputs - Residual maps from WM processing.
    1 - Obs - Model residuals
    2 - Obs - CDM residuals
    3 - Obs - CDM residuals with 3rd trap species.

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
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsModelResG18_19.acs')


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
plt.clf()
meanAC = np.zeros(z2.shape[1])
for i in range(z2.shape[1]):
    meanAC[i] = z2[:,i].mean()
plt.plot(range(z2.shape[1]),meanAC)
plt.xlabel("AL (4x subsampled)")
zeroline72 = np.zeros(72)
plt.plot(range(z2.shape[1]),zeroline72,"r")
plt.show()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AL_flux_residualsWM_obs_modelG18_19.png")

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
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_modelG18_19.png")

#
# ******************************
# Plot2 - Obs-Damaged residuals
#

# Read in array of residuals 
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsDamagedResG18_19.acs')


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
print('Total flux = ', np.nanmean(z2))
plt.clf()
meanAC = np.zeros(z2.shape[1])
for i in range(z2.shape[1]):
    meanAC[i] = z2[:,i].mean()
plt.plot(range(z2.shape[1]),meanAC)
plt.xlabel("AL (4x subsampled)")
zeroline72 = np.zeros(72)
plt.plot(range(z2.shape[1]),zeroline72,"r")
#plt.show()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AL_flux_residualsWM_obs_cdmG18_19.png")
totalmeanAC = np.sum(meanAC[1:72])
print(totalmeanAC)

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
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_cdmG18_19.png")

#
# ******************************
# Plot3 - Obs-Damaged residuals with 3rd slow trap species.
#

# Read in array of residuals 
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091_cdm_bg50_tauE-3._G_mag_18_to_19_AC_1000_to_2000.acs')


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
print('Total flux = ', np.nanmean(z2))
plt.clf()
meanAC = np.zeros(z2.shape[1])
for i in range(z2.shape[1]):
    meanAC[i] = z2[:,i].mean()
plt.plot(range(z2.shape[1]),meanAC)
plt.xlabel("AL (4x subsampled)")
zeroline72 = np.zeros(72)
plt.plot(range(z2.shape[1]),zeroline72,"r")
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AL_flux_residualsWM_obs_cdm_slowSpecies_bg50_tauE-3._G_mag_18_to_19_AC_1000_to_2000.png")
totalmeanAC = np.sum(meanAC[1:72])
print(totalmeanAC)
#plt.show()


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
ax.set_title('Residuals (obs-cdm) with 3rd trap species')

color_map = plt.imshow(z2, origin = "lower", extent = [-9, 9, -6, 6])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_cdm_slowSpeciesG18_19.png")

#
# ******************************
# Plot4 - Overplot of mean flux along AL for 2 and 3rd trap species, to check what difference 3rd trap species makes.
#

# Read in array of residuals with 2 trap species
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsDamagedResG18_19.acs')

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
print('Total flux = ', np.nanmean(z2))
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
meanAC = np.zeros(z2.shape[1])
for i in range(z2.shape[1]):
    meanAC[i] = z2[:,i].mean()
plt.plot(range(z2.shape[1]),meanAC)
plt.xlabel("AL (4x subsampled)")
zeroline72 = np.zeros(72)
plt.plot(range(z2.shape[1]),zeroline72,"r")

meanAL = np.zeros(z2.shape[0])
for i in range(z2.shape[0]):
    meanAL[i] = z2[i,:].mean()
    

# Read in array of residuals with slow species
tdata = ascii.read('/Users/cp232/eclipse-workspace/CTILab/outputs/overSampledArrayObsDamagedResG18_19slowSpecies.acs')

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
z2_3traps = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
print('Total flux = ', np.nanmean(z2_3traps))
meanAC_3traps = np.zeros(z2_3traps.shape[1])
for i in range(z2_3traps.shape[1]):
    meanAC_3traps[i] = z2_3traps[:,i].mean()
plt.plot(range(z2_3traps.shape[1]),meanAC_3traps)
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AL_flux_residualsWM_obs_cdm_2and3species_overplotG18_19.png")
totalmeanAC_3traps = np.sum(meanAC_3traps[1:72])
print(totalmeanAC_3traps)

# Plot of mean flux along AC
meanAL_3traps = np.zeros(z2_3traps.shape[0])
for i in range(z2_3traps.shape[0]):
    meanAL_3traps[i] = z2_3traps[i,:].mean()
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
plt.plot(range(z2.shape[0]),meanAL)
plt.xlabel("AC (4x subsampled)")
zeroline72 = np.zeros(np.size(meanAL))
plt.plot(range(z2.shape[0]),zeroline72,"r")
plt.plot(range(z2_3traps.shape[0]),meanAL_3traps)
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AC_flux_residualsWM_obs_cdm_2and3species_overplotG18_19.png")  

plt.clf()
plt.plot(range(z2.shape[0]),meanAL-meanAL_3traps)
plt.xlabel("AC (4x subsampled)")
plt.ylabel("OBS-CDM (2 traps) - OBS-CDM (3 traps)")
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/AC_flux_residualsWM_obs_cdm_2and3species_differenceG18_19.png")


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
ax.set_title('Residuals (obs-cdm) difference with 2 and 3rd trap species')

color_map = plt.imshow(z2-z2_3traps, origin = "lower", extent = [-9, 9, -6, 6])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/eclipse-workspace/CTILab/outputs/residualsWM_obs_cdm_2_3_species_differenceG18_19.png")
