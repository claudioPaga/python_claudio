#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 11:08:32 2022

@author: cp232
"""


import numpy as np
import matplotlib.pyplot as plt

def diffuse_bg_record(inputfluxfilestore, inputfluxfilenode):
    # Set output plot filenames for the 2 types of plots
    namesplit = inputfluxfilestore.split("store")
    endNameSplit = namesplit[1].split(".")
    ns2 = namesplit[0].split('/')
    filetit = ns2[len(ns2)-2]
    outPlotName = namesplit[0]+'store'+endNameSplit[0]+'losses.png'
    outCTIName = namesplit[0]+'store'+endNameSplit[0]+'losses.txt'

    # Read in the table as a numpy array
    datain = np.loadtxt(inputfluxfilestore, skiprows=1)
    store_transfers = datain[:,0]
    flux = datain[:,1]
    
    xrays = flux[(flux > 80) & (flux < 2000)]
    input_e = np.array([150.0, 500.0, 1600.0, 150.0, 500.0, 1600.0, 150.0, 500.0, 1600.0])
    total_tr = np.array([100, 120, 140, 2000, 2020, 2040, 3500, 3520, 3540]) + 719
    print('Input energies')
    print(input_e)
    print('Total transfers')
    print(total_tr)

    losses = input_e - xrays
    cti = 1. - (1.-losses/xrays)**(1./total_tr)
    print('Losses')
    print(losses)
    # Record to file
    outputarr = np.array([input_e, total_tr, losses, cti])
    outputarrt = outputarr.T
    np.savetxt(outCTIName, outputarrt, header = "Input_flux  N_p_transfers  Losses   CTI", fmt='%i, %i, %1.2f, %1.2e')   
    
    # Simple scatter plot showing losses
    fig, axs = plt.subplots(2, figsize=(12, 6))
    
    fig.suptitle('CTI '+filetit)
    axs[0].plot(input_e, losses, 'o', color='blue')
    axs[1].plot(total_tr, losses, 'o', color= 'red')
    axs[0].set_xlabel("Input X-ray Flux")
    axs[0].set_ylabel("Losses")
    
    axs[1].set_xlabel("Transfers parallel")
    axs[1].set_ylabel("Losses")
        
    #ax.set_ylim([0, 5])

    fig.savefig(outPlotName)
    
    # Evaluate losses in serial register
    storey = datain[:,0]
    storeflux = datain[:,1]
    
    data = np.loadtxt(inputfluxfilenode, skiprows=1)
    serialy = data[:,1]
    serialx = data[:,0]
    nodeflux = data[:,2]
    
    store_xrays  = storeflux[(storeflux > 50) & (storeflux < 2000)]
    store_xrays_y = storey[(storeflux > 50) & (storeflux < 2000)]

    print('Store input flux ', store_xrays)
    for index, store_flux_value in enumerate(store_xrays):
        print('Xray n and store input energy ', index, store_flux_value)
        serial_xrays_flux = nodeflux[(nodeflux > 50.) & (nodeflux < 2000.) & (serialy == store_xrays_y[index])] 
        print('Flux in serial register:', serial_xrays_flux)
        serial_xrays_x = serialx[(nodeflux > 50.) & (nodeflux < 2000.) & (serialy == store_xrays_y[index])] 
        print('SerialX location:', serial_xrays_x)
        losses_serial = store_flux_value -serial_xrays_flux
        print('Losses', losses_serial)
        print('Next...')
    
    