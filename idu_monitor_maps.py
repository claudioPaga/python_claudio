#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 29 Jan 2020

Inputs - Flux maps from IDU monitorning checks..

Outputs - 2D maps plots of IPD, AGIS and difference wrt the assigned window

Summary - compares maps of damaged WMRs generated using different trap species/bg values, plotting the flux difference.

Updates - Proc borrowed from cdm_flux_map.py

@author: cp232
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
from astropy.io import ascii
from scipy.interpolate import griddata

#
#------------------------------
# Read in IDP binned map
#------------------------------
#

# Search ipd filenames 


binnedMapFiles = glob.glob('/Users/cp232/Gaia/data/IPD_MON/OUTPUT_TEST/SubSampleLocDistribution/debug/revs5000_5100/binnedMap*WinClass0*.acs')

for filenameCounter, name in enumerate(binnedMapFiles):
    
    # Set output plot filename for maps
    namesplit = name.split("binned")
    endNameSplit = namesplit[1].split(".")
    outMapPlotName = namesplit[0]+endNameSplit[0]+'.png'
    
    outACHistoPlotName = namesplit[0]+ endNameSplit[0]+'_AChisto.png'
    outALHistoPlotName = namesplit[0]+ endNameSplit[0]+'_ALhisto.png'
    
    # Read in array of counts in the Window
    tdata = ascii.read(name)
    x = tdata['al']
    y = tdata['ac']
    ipd = tdata['ipdN']
    agis = tdata['agisN']
    delta = tdata['deltaN']

    # Prepare the grid for the plot
    xyz = {'x': x, 'y': y, 'z': ipd}
    
    # put the data into a pandas DataFrame (this is what my data looks like)
    df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 
    
    # re-create the 2D-arrays
    x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
    y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
    x2, y2 = np.meshgrid(x1, y1)
    zipd = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
    
    # Calculate mean counts along AC/AL
    meanAL_ipd = np.zeros(zipd.shape[1])
    for i in range(zipd.shape[1]):
        meanAL_ipd[i] = zipd[:,i].mean()
    meanAC_ipd = np.zeros(zipd.shape[0])
    for i in range(zipd.shape[0]):
        meanAC_ipd[i] = zipd[i,:].mean()       
    

    #------------------------------
    # AGIS
    #------------------------------
    
    # Prepare the grid for the plot
    xyz = {'x': x, 'y': y, 'z': agis}
    
    # put the data into a pandas DataFrame (this is what my data looks like)
    df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 
    
    # re-create the 2D-arrays
    x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
    y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
    x2, y2 = np.meshgrid(x1, y1)
    zagis = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
     
    # Calculate mean counts along AC/AL
    meanAL_agis = np.zeros(zagis.shape[1])
    for i in range(zagis.shape[1]):
        meanAL_agis[i] = zagis[:,i].mean()
    meanAC_agis = np.zeros(zagis.shape[0])
    for i in range(zagis.shape[0]):
        meanAC_agis[i] = zagis[i,:].mean()       
    
    
    #------------------------------
    # Delta IDP-AGIS binned map
    #------------------------------
    
    
    # Prepare the grid for the plot
    xyz = {'x': x, 'y': y, 'z': delta}
    
    # put the data into a pandas DataFrame (this is what my data looks like)
    df = pd.DataFrame(xyz, index=range(len(xyz['x']))) 
    
    # re-create the 2D-arrays
    x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
    y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
    x2, y2 = np.meshgrid(x1, y1)
    zdelta = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
    
     # Calculate mean counts along AC/AL
    meanAL_delta = np.zeros(zdelta.shape[1])
    for i in range(zdelta.shape[1]):
        meanAL_delta[i] = zdelta[:,i].mean()
    meanAC_delta = np.zeros(zdelta.shape[0])
    for i in range(zdelta.shape[0]):
        meanAC_delta[i] = zdelta[i,:].mean()       
    
    
    # Plot of all the images in a multiplot, stacked horizontally
    
    w = 20  
    h = 8
    d = 70
    fig = plt.figure(figsize=(w, h), dpi=d)
    ima = [x,
         y,
         delta]
    plt.clf()
    plt.subplot(1, 3, 1) # = 1 column, 3 rows, 1 = the first position
    ax = fig.gca()
    ax.set_xlabel('AL')
    ax.set_ylabel('AC')
    color_map = plt.imshow(zipd, origin = "lower", extent = [0, 19, 0, 19])
    color_map.set_cmap("Blues")
    plt.xlim(4, 15)
    plt.ylim(4, 15)
    plt.title("IPD map")
    plt.colorbar(orientation="horizontal")
    plt.subplot(1, 3, 2) # = 1 column, 3 rows, 2 = the second position
    ax = fig.gca()
    color_map = plt.imshow(zagis, origin = "lower", extent = [0, 19, 0, 19])
    color_map.set_cmap("Blues")
    plt.xlim(4, 15)
    plt.ylim(4, 15)
    ax = fig.gca()
    ax.set_xlabel('AL')
    ax.set_ylabel('AC')
    plt.title("AGIS map")
    plt.colorbar(orientation="horizontal")
    plt.subplot(1, 3, 3) # = 1 column, 3 rows, 3 = the third position
    color_map = plt.imshow(zdelta, origin = "lower", extent = [0, 19, 0, 19])
    color_map.set_cmap("Blues")
    plt.xlim(4, 15)
    plt.ylim(4, 15)
    ax = fig.gca()
    ax.set_xlabel('AL')
    ax.set_ylabel('AC')
    plt.title("IPD-AGIS map")
    plt.colorbar(orientation="horizontal")    
    plt.savefig(outMapPlotName)
    plt.close()
    
    
    # Plot of Cumulative histo along AC   
    
    plt.clf()
    w = 12
    h = 8
    d = 70
    fig = plt.figure(figsize=(w, h), dpi=d)
    plt.plot(range(zipd.shape[0]),meanAC_ipd, "-y", label = 'IPD')
    plt.xlabel("AC ")
    plt.ylabel("Mean n. of Sources")
    zeroline72 = np.zeros(np.size(meanAC_ipd))
    #plt.plot(range(z2.shape[0]),zeroline72,"r")
    plt.plot(range(zagis.shape[0]),meanAC_agis, "-r", label = 'AGIS')
    plt.legend(loc = "upper left")
    plt.savefig(outACHistoPlotName)
    plt.close()
    
    # Plot of Cumulative histo along AL    
    
    plt.clf()
    w = 12
    h = 8
    d = 70
    fig = plt.figure(figsize=(w, h), dpi=d)
    plt.plot(range(zipd.shape[0]),meanAL_ipd, "-y", label = 'IPD')
    plt.xlabel("AL ")
    plt.ylabel("Mean n. of Sources")
    zeroline72 = np.zeros(np.size(meanAL_ipd))
    #plt.plot(range(z2.shape[0]),zeroline72,"r")
    plt.plot(range(zagis.shape[0]),meanAL_agis, "-r", label = 'AGIS')
    plt.legend(loc = "upper left")
    plt.savefig(outALHistoPlotName)
    plt.close()
    
    
    

# Plot of IPD map

#plt.clf()
#w = 8
#h = 8
#d = 70
#fig = plt.figure(figsize=(w, h), dpi=d)
#ima = [x,
#     y,
#     ipd]
#ax = fig.gca()
#ax.set_xlabel('AL')
#ax.set_ylabel('AC')
#ax.set_title('IPD Counts')
#
#color_map = plt.imshow(zipd, origin = "lower", extent = [0, 19, 0, 19])
#color_map.set_cmap("Blues")
#plt.xlim(5, 15)
#plt.ylim(5, 15)
#plt.colorbar()
#plt.savefig("/Users/cp232/Gaia/data/IPD_MON/OUTPUT_TEST/SubSampleLocDistribution/debug/ipdMap_AF3_ROW0_FoV0_WinClass0.png")


# Plot oa AGIS map

#w = 8
#h = 8
#d = 70
#fig = plt.figure(figsize=(w, h), dpi=d)
#ima = [x,
#     y,
#     agis]
#ax = fig.gca()
#plt.xlim(5, 15)
#plt.ylim(5, 15)
#ax.set_xlabel('AL')
#ax.set_ylabel('AC')
#ax.set_title('AGIS Counts')
#
##color_map = plt.imshow(z2, origin = "lower", extent = [0, 20, 0, 20], clim=[0, 100])
#color_map = plt.imshow(zagis, origin = "lower", extent = [0, 19, 0, 19])
#color_map.set_cmap("Blues")
#plt.colorbar()
#plt.savefig("/Users/cp232/Gaia/data/IPD_MON/OUTPUT_TEST/SubSampleLocDistribution/debug/agisMap_AF3_ROW0_FoV0_WinClass0.png")


# Plot of difference (IPD-AGIS)

#w = 8
#h = 8
#d = 70
#fig = plt.figure(figsize=(w, h), dpi=d)
#ima = [x,
#     y,
#     delta]
#ax = fig.gca()
#ax.set_xlabel('AL')
#ax.set_ylabel('AC')
#ax.set_title('IPD-AGIS Counts')
#plt.xlim(5, 15)
#plt.ylim(5, 15)
#
##color_map = plt.imshow(z2, origin = "lower", extent = [0, 20, 0, 20], clim=[0, 100])
#color_map = plt.imshow(zdelta, origin = "lower", extent = [0, 19, 0, 19])
#color_map.set_cmap("Blues")
#plt.colorbar()
#plt.savefig("/Users/cp232/Gaia/data/IPD_MON/OUTPUT_TEST/SubSampleLocDistribution/debug/deltaMap_AF3_ROW0_FoV0_WinClass0.png")
