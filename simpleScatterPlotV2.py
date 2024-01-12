#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 14:50:27 2022

@author: cp232
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 12:30:27 2022

Reads in a table consisting of [TransferStoreN, Flux]
reporting the flux at the bottomw of the SXI store section during phaseB of an image readout.
It's the binned lines flux that goes into the serial register after the CDM simulation of a diffuse BG + X-rays at specific energies and 
@author: cp232
"""

import numpy as np
import matplotlib.pyplot as plt

def simpleScatterPlotV2(inputfluxfile):
    # Set output plot filenames for the 2 types of plots
    namesplit = inputfluxfile.split("store")
    endNameSplit = namesplit[1].split(".")
    outPlotName = namesplit[0]+'store'+endNameSplit[0]+'plot.png'
    #outScatterName = namesplit[0]+'store'+endNameSplit[0]+'scatter.png'

    # Read in the table as a numpy array
    datain = np.loadtxt(inputfluxfile, skiprows=1)
    store_transfers = datain[:,0]
    flux = datain[:,1]
    # Simple scatter plot with plt
    plt.figure(figsize=(12,6))    
    plt.plot(store_transfers, flux, 'o', color = 'red')
    plt.xlabel("Store Transfers in Phase B (Store readout) ")
    plt.ylabel("Flux [e-/pixel]")
    plt.ylim([0, 5])
    plt.savefig(outPlotName)

#    plt.plot(store_transfers, flux, '-ok', color = 'red')
#    plt.ylim(0, 5)
#    plt.show()
#    plt.savefig(outPlotName)
#    plt.close()