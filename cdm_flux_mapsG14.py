#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:45:15 2019

Inputs - Flux maps from WM processing.

Outputs - 2D maps plots + Median cumulative flux along AC, to check for CTI signatured   

Summary - compares maps of damaged WMRs generated using different trap species/bg values, plotting the flux difference.

Updates - The code is copied over from cdm_flux_maps.py, to work with input G=14 files

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

# Read in array of damaged flux values from 1st map.
# 3 trap species,  BG = 0 (or low BG, eg bg = 10)

tdata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091_cdm_bg0_rho05_sigma500k_tauE-3._G_mag_14_to_15_AC_1000_to_2000_damaged.acs')

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

meanAL = np.zeros(z2.shape[0])
for i in range(z2.shape[0]):
    meanAL[i] = z2[i,:].mean()
    
#
#------------------------------
# Damaged flux, BG = 100
#------------------------------
#

# Read in array of damaged flux with BG = 100
tdata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091_cdm_bg100_rho05_sigma500k_tauE-3._G_mag_14_to_15_AC_1000_to_2000_damaged.acs')

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
z2_bg100 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
print('Total flux = ', np.nanmean(z2_bg100))
meanAC_bg100 = np.zeros(z2_bg100.shape[1])
for i in range(z2_bg100.shape[1]):
    meanAC_bg100[i] = z2_bg100[:,i].mean()
    
meanAL_bg100 = np.zeros(z2_bg100.shape[0])
for i in range(z2_bg100.shape[0]):
    meanAL_bg100[i] = z2_bg100[i,:].mean()
    
#
#------------------------------
# Damaged flux, 2 TRAP SPECIES, BG added to WRMS
#------------------------------
#

# Read in array of damaged flux with 2 TRAP SPECIES
tdata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091_bckAdded_cdm_LowSignalConstraint._G_mag_14_to_15_AC_1000_to_2000_damaged.acs')

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
z2_2traps = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
print('Total flux = ', np.nanmean(z2_2traps))
meanAC_2traps = np.zeros(z2_2traps.shape[1])
for i in range(z2_2traps.shape[1]):
    meanAC_2traps[i] = z2_2traps[:,i].mean()
    
   
    
#
#------------------------------
# Observed flux
#------------------------------
#

# Read in array of Model Flux
observedfluxdata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091._G_mag_14_to_15_AC_1000_to_2000_observed.acs')    

observedfluxX = observedfluxdata['alIdx']
observedfluxY = observedfluxdata['acIdx']
observedfluxZ = observedfluxdata['flux']
# Remove NANs
filter = np.logical_not(np.isnan(observedfluxZ))
observedfluxx = observedfluxX
observedfluxy = observedfluxY
observedfluxz = observedfluxZ
observedfluxx = observedfluxx[filter]
observedfluxy = observedfluxy[filter]
observedfluxz = observedfluxz[filter]    
observedfluxxyz = {'x': observedfluxx, 'y': observedfluxy, 'z': observedfluxz}
# put the data into a pandas DataFrame (this is what my data looks like)
observedfluxdf = pd.DataFrame(observedfluxxyz, index=range(len(observedfluxxyz['x']))) 
# re-create the 2D-arrays
observedfluxx1 = np.linspace(observedfluxdf['x'].min(), observedfluxdf['x'].max(), len(observedfluxdf['x'].unique()))
observedfluxy1 = np.linspace(observedfluxdf['y'].min(), observedfluxdf['y'].max(), len(observedfluxdf['y'].unique()))
observedfluxx2, observedfluxy2 = np.meshgrid(observedfluxx1, observedfluxy1)
z2_observedflux = griddata((observedfluxdf['x'], observedfluxdf['y']), observedfluxdf['z'], (observedfluxx2, observedfluxy2), method='cubic')
print('Total flux = ', np.nanmean(z2_observedflux))
meanAC_observedflux = np.zeros(z2_observedflux.shape[1])
for i in range(z2_observedflux.shape[1]):
    meanAC_observedflux[i] = z2_observedflux[:,i].mean()
 # Plot of mean flux along AC
meanAL_observedflux = np.zeros(z2_observedflux.shape[0])
for i in range(z2_observedflux.shape[0]):
    meanAL_observedflux[i] = z2_observedflux[i,:].mean()       
        
    
#
#------------------------------
# Model flux, BG = 0
#------------------------------
#

# Read in array of Model Flux
modeldata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091._G_mag_14_to_15_AC_1000_to_2000_model.acs')
#modeldata = ascii.read('/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/FOV1_ROW4_AF5_NOGATE_WINCLASS0_AC_WindowModelRecord_51091._G_mag_13_to_14_AC_1000_to_2000_model.acs')

modelX = modeldata['alIdx']
modelY = modeldata['acIdx']
modelZ = modeldata['flux']
# Remove NANs
filter = np.logical_not(np.isnan(modelZ))
modelx = modelX
modely = modelY
modelz = modelZ
modelx = modelx[filter]
modely = modely[filter]
modelz = modelz[filter]

modelxyz = {'x': modelx, 'y': modely, 'z': modelz}

# put the data into a pandas DataFrame (this is what my data looks like)
modeldf = pd.DataFrame(modelxyz, index=range(len(modelxyz['x']))) 

# re-create the 2D-arrays
modelx1 = np.linspace(modeldf['x'].min(), modeldf['x'].max(), len(modeldf['x'].unique()))
modely1 = np.linspace(modeldf['y'].min(), modeldf['y'].max(), len(modeldf['y'].unique()))
modelx2, modely2 = np.meshgrid(modelx1, modely1)
z2_model = griddata((modeldf['x'], modeldf['y']), modeldf['z'], (modelx2, modely2), method='cubic')
print('Total flux = ', np.nanmean(z2_model))
meanAC_model = np.zeros(z2_model.shape[1])
for i in range(z2_model.shape[1]):
    meanAC_model[i] = z2_model[:,i].mean()
 # Plot of mean flux along AC
meanAL_model = np.zeros(z2_model.shape[0])
for i in range(z2_model.shape[0]):
    meanAL_model[i] = z2_model[i,:].mean()       
    
    
# Plots the LSF/2-traps/3-traps BG=0/3-traps BG=100 profiles along AL (non really interesting for my analysis, just for sanity)   
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
plt.plot(range(z2.shape[1]),meanAC, "-b", label = 'BG = 0') #CDM BG = 0
plt.xlabel("AL (4x subsampled)")
plt.ylabel("Mean damagedFlux")
plt.plot(range(z2_bg100.shape[1]),meanAC_bg100, "-r", label = 'BG = 100') #CDM BG = 100
plt.plot(range(z2_model.shape[1]),meanAC_model, "-y", label = 'PSF') #PSF
plt.plot(range(z2_2traps.shape[1]),meanAC_2traps, "-g", label = '2 TRAPS') #CDM, 2 TRAPS
plt.legend(loc = "upper left")
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AL_damagedFlux_2traps_model_bg0_bg100_overplotG14_15.png")
totalmeanAC_bg100 = np.sum(meanAC_bg100[1:72])
print(totalmeanAC_bg100)

# Plot of mean flux along AC
meanAL_bg100 = np.zeros(z2_bg100.shape[0])
for i in range(z2_bg100.shape[0]):
    meanAL_bg100[i] = z2_bg100[i,:].mean()

    
meanAL_2traps = np.zeros(z2_2traps.shape[0])
for i in range(z2_2traps.shape[0]):
    meanAL_2traps[i] = z2_2traps[i,:].mean()
    
    
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
plt.plot(range(z2_model.shape[0]),meanAL_model, "-y", label = 'PSF, G=14')
plt.xlabel("AC (4x subsampled)")
plt.ylabel("Mean Flux")
zeroline72 = np.zeros(np.size(meanAL))
#plt.plot(range(z2.shape[0]),zeroline72,"r")
plt.plot(range(z2_bg100.shape[0]),meanAL_bg100, "-r", label = '3 TRAPS, BG = 100')
plt.plot(range(z2.shape[0]),meanAL, "-b", label = '3 TRAPS, BG = 0')
plt.plot(range(z2_2traps.shape[0]),meanAL_2traps, "-g", label = '2 TRAPS')
plt.legend(loc = "upper left")
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AC_damagedFlux_2traps_model_bg0_bg100_rho05_sigma500k_tauE-3_overplotG14_15.png")


# Write the profiles out
ac = np.linspace(1,48,48)
tableACprofiles = {'AC': ac, 'Model': meanAL_model, 'DamagedBG0': meanAL, 'DamagedBG100': meanAL_bg100, 'Damaged2Traps': meanAL_2traps, 'Observed': meanAL_observedflux}
ascii.write(tableACprofiles, '/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/ACFluxProfilesG14_15.dat', formats={'AC': '%i', 'Model': '%.5f', 'DamagedBG0': '%.5f', 'DamagedBG100': '%.5f', 'Damaged2Traps': '%.5f', 'Observed': '%.5f'})


#OverPlot the median flux AC profiles of model/CDM(bg100)/CDM(BG0)/Observed data
plt.clf()
w = 12
h = 8
d = 70
fig = plt.figure(figsize=(w, h), dpi=d)
plt.plot(range(z2_model.shape[0]),meanAL_model, "-y", label = 'PSF, G=14')
plt.xlabel("AC (4x subsampled)")
plt.ylabel("Mean Flux")
zeroline72 = np.zeros(np.size(meanAL))
#plt.plot(range(z2.shape[0]),zeroline72,"r")
plt.plot(range(z2_bg100.shape[0]),meanAL_bg100, "-r", label = '3 TRAPS, BG = 100')
plt.plot(range(z2.shape[0]),meanAL, "-b", label = '3 TRAPS, BG = 0')
plt.plot(range(z2_observedflux.shape[0]),meanAL_observedflux, "-g", label = 'Observed')
plt.legend(loc = "upper left")
#plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AC_damagedFlux_bg10_bg100_overplotG17_18.png")  
#plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AC_damagedFlux_bg0_bg100_overplotG13_14.png")  
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AC_damagedFlux_observed_model_bg_0_bg100_rho05_sigma500k_tauE-3_overplotG14_15.png")




# Plot of LSF difference between BG=0 and BG=100
plt.clf()
fig = plt.figure(figsize=(w, h), dpi=d)
plt.plot(range(z2.shape[0]),meanAL-meanAL_bg100)
plt.xlabel("AC (4x subsampled)")
plt.ylabel("Mean damaged flux difference (BG=0 - BG=100)")
plt.legend(loc = "upper left")
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/AC_damagedFlux_bg0_bg100_rho05_sigma500k_tauE-3_differenceG14_15.png")


#Plot the damage flux difference between BG=0 and BG=100

# Set width and height of image (it might be in cm) and the resolution
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x, y, z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Damaged flux difference (BG = 0 - BG = 100)')
color_map = plt.imshow(z2-z2_bg100, origin = "lower", extent = [-9, 9, -6, 6])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/damagedFlux_bg0_bg100_rho05_sigma500k_tauE-3_differenceG14_15.png")


#Plot the flux difference between PSF and damaged flux when BG=0

# Set width and height of image (it might be in cm) and the resolution
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Flux difference, PSF - damaged(BG = 0)')

color_map = plt.imshow(z2_model-z2, origin = "lower", extent = [-9, 9, -6, 6], clim=[-140, 300])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/psf_damagedFlux_difference_bg0_rho05_sigma500k_tauE-3_G14_15.png")

print('Mean flux difference, WRM processed with BG=0 = ', np.nanmean(z2_model-z2))
print('Total flux difference, WRM processed with BG=0 = ', np.nansum(z2_model-z2))

#Plot the flux difference between PSF and damaged flux when BG=100

# Set width and height of image (it might be in cm) and the resolution
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Flux difference, PSF - damaged(BG = 100)')

color_map = plt.imshow(z2_model - z2_bg100, origin = "lower", extent = [-9, 9, -6, 6], clim=[-140, 300])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/psf_damagedFlux_difference_bg100_rho05_sigma500k_tauE-3_G14_15.png")

print('Mean flux difference, WRM processed with BG=100 = ', np.nanmean(z2_model-z2_bg100))
print('Total flux difference, WRM processed with BG=100 = ', np.nansum(z2_model-z2_bg100))

# Plot the difference between the 2-traps vs 3 traps (low BG) case

# Set width and height of image (it might be in cm) and the resolution
fig = plt.figure(figsize=(w, h), dpi=d)
ima = [x,
     y,
     z]
ax = fig.gca()
ax.set_xlabel('AL')
ax.set_ylabel('AC')
ax.set_title('Flux difference, 3-2 traps')

color_map = plt.imshow(z2 - z2_2traps, origin = "lower", extent = [-9, 9, -6, 6], clim=[-100, 10])
color_map.set_cmap("Blues")
plt.colorbar()
plt.savefig("/Users/cp232/Gaia/data/WindowModelRecords/data/medianFluxMaps/Flux_difference_2traps_3traps_rho05_sigma500k_tauE-3_G14_15.png")


print('Mean flux difference, WRM processed with BG=100 = ', np.nanmean(z2 - z2_2traps))
print('Total flux difference, WRM processed with BG=100 = ', np.nansum(z2 - z2_2traps))
