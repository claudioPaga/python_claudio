#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:10:32 2022

Summary - Reads in a text file consisting of trap filling values in the serial register during readout.


@author: cp232
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def main():
    # Number of frames whose flux is sotred in input file is 719
    nframes = 719
    #flux = np.loadtxt('/Users/cp232/IDLWorkspace/optical_ci/t1_serial_filled_frn0.txt') 
    flux = np.loadtxt('/Users/cp232/IDLWorkspace/optical_ci/t1_serial_filled_frn1.txt') 
    storeTr = flux[:,0]
    detX = flux[:,1]
    fill = flux[:,2]
    
    detXL = 2255
    detX0 = np.arange(detXL)
    fill0 = fill[0:detXL]
    

    fig = plt.figure()
    # Set the axis of the plot
    # ax = plt.axes(ylim=(0.01, 0.11), xlim=(0, 3790), title = 'Flux = 1e-', xlabel = 'DETY', ylabel = 'Flux')
    ax = plt.axes(ylim=(0.0, 0.2), xlim=(0, detXL), title = 'E-centre filling', xlabel = 'DETX', ylabel = 'Filled traps')
    #ax = plt.axes(ylim=(0.6, 1.1), xlim=(0, 3790), title = 'Flux = 1e-, T=-100C', xlabel = 'DETY', ylabel = 'Flux')
    # Initial scatter plot for flux values at beginning of readout
    scat = plt.scatter(detX0, fill0, s=1)

    # Animation    
    ani = animation.FuncAnimation(fig, update_plot, frames=range(nframes-1),
                                  fargs=(detX0, fill, ax, scat))
     # Save animation to file
    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/image_evo_animation_scatter.gif', writer='imagemagick', fps=60)
    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux5_frames3/image_evo_animation_scatter.gif', writer='imagemagick', fps=300)
    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux1_frames2_T-100C/image_evo_animation_scatter.gif', writer='imagemagick', fps=300)
    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux5fr3T-120/image_evo_fr0_animation_scatter.gif', writer='imagemagick', fps=300)
    #ani.save('/Users/cp232/IDLWorkspace/optical_ci/fill1_serial_animation_scatter_frn0.gif', writer='imagemagick', fps=40)
    ani.save('/Users/cp232/IDLWorkspace/optical_ci/fill1_serial_animation_scatter_frn1.gif', writer='imagemagick', fps=40)
    
    plt.show()
    
    # Quicker animation every 10 frames
    frequency = 10
    ani = animation.FuncAnimation(fig, update_plot, frames=range(0,nframes-1,frequency),
                                  fargs=(detX0, fill, ax, scat))
     # Save animation to file
    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/image_evo_animation_scatter_step10.gif', writer='imagemagick', fps=60)
    scat = plt.scatter(detX0, fill0, s=1)     

    #ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux5_frames3/image_evo_animation_scatter_step10.gif', writer='imagemagick', fps=20)
    #ani.save('/Users/cp232/IDLWorkspace/optical_ci/fill1_serial_animation_scatter_step10_frn0.gif', writer='imagemagick', fps=20)
    ani.save('/Users/cp232/IDLWorkspace/optical_ci/fill1_serial_animation_scatter_step10_frn1.gif', writer='imagemagick', fps=20)
    
    plt.show()

def update_plot(i, detX0, fill, ax, scat):
   # Select the flux along the column during transfer i
    fill_plot =  fill[2255*i:2255*i+2255]
#    ax.set_title('Flux = 5e-, TransferN='+ str(i))
    ax.set_title('E-centre filling Frame 1- TransferN='+ str(i))
    data = [detX0, fill_plot]
    b = np.transpose(data)
    scat.set_offsets(b)
    return scat,

main()