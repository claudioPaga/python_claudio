#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:30:50 2019

Plotting, labels, legend, log, multiplots, histograms

@author: cp232
"""

def practice17plots():

    import numpy as np
    import matplotlib.pyplot as plt
    print("Test")
    x = np.arange(-5, 5, 0.1)
    y = np.sin(x)   
    y2 = np.cos(x)
    plt.plot(x, y);
# Note - ; avoids printing out the matplotlib.pyplot object type

    plt.plot(x,y, 'sr-', markersize = 1, linewidth = 2, label = 'Sin(x)')
    plt.plot(x,y2, 'ob-', markersize = 1, linewidth = 2, label = 'Cos(x)')
    plt.xlabel("X")
    plt.ylabel("\sigma") # NOTE - myplot knows latex!!!
    plt.legend(loc = "upper right")
    plt.axis([-5, 5, -1.1, 1.1]) # NOTE plt.axis([xmin, xmax, ymin, ymax])
    plt.savefig("plotExamplePractice17.png")
    y3 = x**2
    plt.show()
    plt.loglog(x,y3,"bo-")
    plt.subplot(231)
    plt.plot(x,y, 'sr-', markersize = 1, linewidth = 2, label = 'Sin(x)')
    plt.subplot(232)
    plt.plot(x,y2, 'ob-', markersize = 1, linewidth = 2, label = 'Cos(x)')
    plt.subplot(233)
    plt.loglog(x,y3,"bo-")
    plt.subplot(234)
    plt.loglog(x,y3,"bo-")
    
    
    x = np.random.normal(size=10000)
    plt.hist(x)
