#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:47:55 2021

Summary - Reads in a text file consisting of flux values along a column during readout. the column is 3791 pixels long, therefore flux values
          after exposure are in flux[0:3790], after the first transfer are in flux[3791:7582] etc. This flux change is what is animated in a scatter plot

@author: cp232
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def main():
    # Number of frames whose flux is sotred in input file is 719*6
    numframes = 4314
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux01_frames2/ima_evo_frn0_tr0_10.txt')
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux01_frames2/ima_evo_frn0.txt')
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux5_frames3/ima_evo_frn0.txt') 
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux1_frames2_T-100C/ima_evo_frn0.txt') 
#    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux5fr3T-120/ima_evo_frn0.txt') 
    flux = np.loadtxt('/Users/cp232/SMILE/TestRadHardness/diffuse_bg/test_runs_pixel/flux5fr3T-80/ima_evo_frn2.txt') 
    transfers = np.arange(3790)
    flux_fr0 = flux[0:3790]    

    fig = plt.figure()
    # Set the axis of the plot
    # ax = plt.axes(ylim=(0.01, 0.11), xlim=(0, 3790), title = 'Flux = 1e-', xlabel = 'DETY', ylabel = 'Flux')
    ax = plt.axes(ylim=(4.5, 5.1), xlim=(0, 3790), title = 'Flux = 5e-, T=-80C', xlabel = 'DETY', ylabel = 'Flux')
    #ax = plt.axes(ylim=(0.6, 1.1), xlim=(0, 3790), title = 'Flux = 1e-, T=-100C', xlabel = 'DETY', ylabel = 'Flux')
    # Initial scatter plot for flux values at beginning of readout
    scat = plt.scatter(transfers, flux_fr0, s=1)

    # Animation    
    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes),
                                  fargs=(transfers, flux, ax, scat))
     # Save animation to file
    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux01_frames2/image_evo_animation_scatter.gif', writer='imagemagick', fps=60)
    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux5_frames3/image_evo_animation_scatter.gif', writer='imagemagick', fps=300)
    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux1_frames2_T-100C/image_evo_animation_scatter.gif', writer='imagemagick', fps=300)
    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux5fr3T-120/image_evo_fr0_animation_scatter.gif', writer='imagemagick', fps=300)
    ani.save('/Users/cp232/SMILE/TestRadHardness/diffuse_bg/test_runs_pixel/flux5fr3T-80/image_evo_fr2_animation_scatter.gif', writer='imagemagick', fps=300)
    plt.show()
    
    # Quicker animation every 10 frames
    frequency = 10
    ani = animation.FuncAnimation(fig, update_plot, frames=range(0,numframes,frequency),
                                  fargs=(transfers, flux, ax, scat))
     # Save animation to file
    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux01_frames2/image_evo_animation_scatter_step10.gif', writer='imagemagick', fps=60)
    scat = plt.scatter(transfers, flux_fr0, s=1)     

    #ani.save('/Users/cp232/SMILE/TestHarness/diffuse_bg/test_runs_pixel/flux5_frames3/image_evo_animation_scatter_step10.gif', writer='imagemagick', fps=20)
    ani.save('/Users/cp232/SMILE/TestRadHardness/diffuse_bg/test_runs_pixel/flux5fr3T-80/image_evo_fr2_animation_scatter_step10.gif', writer='imagemagick', fps=20)
    plt.show()

def update_plot(i, transfers, flux, ax, scat):
   # Select the flux along the column during transfer i
    flux_plot =  flux[3791*i:3791*i+3790]
#    ax.set_title('Flux = 5e-, TransferN='+ str(i))
    ax.set_title('Flux = 5e-, T=-80C, TransferN='+ str(i))
    data = [transfers, flux_plot]
    b = np.transpose(data)
    scat.set_offsets(b)
    return scat,

main()