#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:32:14 2023

Analysis of Signal, BgSignal distributons of outliers with large AcCorr values in DIPD42 test run (revs 4400-4600) generated 
by correction libraries derived during March/April 2023

The goal is to verify that outliers are at the extreme of the Signal and BgSignal distribution

Use a debug file from ROW2 as input file

Example

ctiMonitorOutliers('/Users/cp232/Gaia/CTIMonitor/EmpiricalSCTICorrections/DebugIPD42Revs4400_4600/Debug_AcLocationBiasByMuAndSignalRow2.asc')

@author: cp232
"""

import math
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
from itertools import compress

def ctiMonitorOutliers(debugFile):
    """Plots outliers Signal, BgSignal distribution
    """
    # Read in table, converting fields to integer
    data = np.loadtxt(debugFile) 
    
    # Select OUTliers with large AcCorr values 
    
    acCorr = data[:,13].ravel()
    largeCorrBool = abs(acCorr) > 2
    largeCorrIndex = list(compress(range(len(largeCorrBool)), largeCorrBool))
    
    ccdStrip = data[:,1].astype(int).ravel()
    fov = data[:,2].astype(int).ravel()
    gate = data[:,3].astype(int).ravel()
    gateL = data[:,4].astype(int).ravel()
    texp = data[:,5].ravel()
    gmag = data[:,6].ravel()
    flux = data[:,7].ravel()
    signal = data[:,8].ravel()
    mu = data[:,9].ravel()
    bgsignal = data[:,10].ravel()
    acSmear = data[:,11].ravel()
    acOffset = data[:,12].ravel()
    

    # Setup output figure
    w = 12
    h = 8
    d = 70
    fig = plt.figure(figsize=(w, h), dpi=d)
 
    # Loop over outliers, selecting their CalUnit and displaying the Signal v BgSignal

    for i in largeCorrIndex:
        ccdStripOL = data[i,1].astype(int).ravel()[0]
        fovOL = data[i,2].astype(int).ravel()[0]
        gateOL = data[i,3].astype(int).ravel()[0]
        
        # CCD and node selection
        
        strip_index = ccdStrip == ccdStripOL
        fov_index = fov == fovOL
        gate_index = gate == gateOL

        strip_fov_index = np.logical_and(strip_index, fov_index)
        calU_index = np.logical_and(strip_fov_index, gate_index)        
        calU_index_OL = np.logical_and(calU_index, largeCorrBool)   
    
        stripCU = ccdStrip[calU_index]
        fovCU = fov[calU_index]
        gateCU = gate[calU_index]
        signalCU = signal[calU_index]
        bgsignalCU = bgsignal[calU_index]
        acCorrCU = acCorr[calU_index]       

        signalCU_OL = signal[calU_index_OL]
        bgsignalCU_OL = bgsignal[calU_index_OL]


         
        # Use object-oriented interface to produce a simple scatter plot with plt    
        fig, ax = plt.subplots(figsize=(12, 6))
        outPlotTitle = 'CU Strip'+str(ccdStripOL)+' FoV'+str(fovOL)+' Gate'+str(gateOL)
        ax.plot(signalCU, bgsignalCU, '+', color='blue')
        ax.plot(signalCU_OL, bgsignalCU_OL, 'o', color='red')
        ax.set_title(outPlotTitle)
        ax.set_xlabel("Signal")
        ax.set_ylabel("BgSignal")
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.show()        
        
        outPlotName = 'ol_strip'+str(ccdStripOL)+'_fov'+str(fovOL)+'_gate'+str(gateOL)+'.png'
        
        fig.savefig(outPlotName)
        plt.close()

        
ctiMonitorOutliers('/Users/cp232/Gaia/CTIMonitor/EmpiricalSCTICorrections/DebugIPD42Revs4400_4600/Debug_AcLocationBiasByMuAndSignalRow2.asc')                
#ctiMonitorOutliers('/Users/cp232/Gaia/CTIMonitor/EmpiricalSCTICorrections/DebugIPD42Revs4400_4600/test100.asc')             
        
        
        
        

