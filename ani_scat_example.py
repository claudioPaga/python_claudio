#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:47:55 2021

@author: cp232
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def main():
    numframes = 4314
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/ima_evo_frn0_tr0_10.txt')
    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/ima_evo_frn0.txt')
    transfers = np.arange(3790)
    flux_fr0 = flux[0:3790]    

    fig = plt.figure()
    ax = plt.axes(ylim=(0.01, 0.11), xlim=(0, 3790), title = 'Flux = 1e-', xlabel = 'DETY', ylabel = 'Flux')
    scat = plt.scatter(transfers, flux_fr0, s=1)

#    scat = plt.scatter(x, y, c=c, s=20)

    #ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes),
    #                             fargs=(color_data, scat))
    
    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes),
                                  fargs=(transfers, flux, ax, scat))
     # Save animation to file
    ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux1_frames2/image_evo_animation_sc1.gif', writer='imagemagick', fps=60)

    plt.show()

def update_plot(i, transfers, flux, ax, scat):
   # scat.set_offsets(data[i])
    flux_plot =  flux[3791*i:3791*i+3790]
    ax.set_title('Flux = 1e-, TransferN='+ str(i))
    data = [transfers, flux_plot]
    b = np.transpose(data)
    scat.set_offsets(b)
    return scat,

main()