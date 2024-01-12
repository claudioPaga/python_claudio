#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:47:23 2021

@author: cp232

A simple example of an animated plot
Summary - Reads in a text file consisting of flux values along a column during readout. the column is 3791 pixels long, therefore flux values
          after exposure are in flux[0:3790], after the first transfer are in flux[3791:7582] etc. This flux change is what is animated
@author: cp232

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux01_frames2/ima_evo_frn0.txt')


transfers = np.arange(3790)
fig, ax = plt.subplots()


fig = plt.figure()
ax = plt.axes(ylim=(0.01, 0.11), xlim=(0, 3790), title = 'Flux = 1e-', xlabel = 'DETY', ylabel = 'Flux')
flux_fr0 = flux[0:3790]  
# The animation will expect a sequence of artists, thus the trailing comma.  
line, = ax.plot(transfers, flux_fr0)

def animate(i):
    xdata = transfers
    ydata = flux[3791*i:3791*i+3790]
    ax.set_title('Flux = 1e-, TransferN='+ str(i))
    line.set_data(xdata, ydata)
    return line,

#Init only required for blitting to give a clean slate.
def init():
    line.set_data([], [])
    return line,

# Method to create animation, calling animate for a number of times equals to frames
ani = animation.FuncAnimation(fig, animate, frames=10, init_func=init,
    interval=20, blit=True)

# Save animation to file
ani.save('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux1_frames2/image_evo_animation.gif', writer='imagemagick', fps=60)


plt.show()