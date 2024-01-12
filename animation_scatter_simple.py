#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:28:23 2021

@author: cp232
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def main():
    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/ima_evo_frn0_tr0_10.txt')
    transfers = np.arange(3790)
    flux_fr0 = flux[0:3790] 
    fig, ax = plt.subplots()


    fig = plt.figure()
    ax = plt.axes(ylim=(0.01, 0.11), xlim=(0, 3790), title = 'Flux = 1e-', xlabel = 'DETY', ylabel = 'Flux')
    numframes = 10
    numpoints = 10
    color_data = np.random.random((numframes, numpoints))
    x, y, c = np.random.random((3, numpoints))

    scat = plt.scatter(transfers, flux_fr0, s=100)
  #  scat = plt.scatter(x, y, s=100)

    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes),
                                  fargs=(color_data, scat))
    # Save animation to file
    ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux1_frames2/image_evo_animation.gif', writer='imagemagick', fps=60)
    plt.show()

def update_plot(i, data, scat):
    scat.set_offsets(data[i])
    return scat,

main()

