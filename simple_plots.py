#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 16:41:31 2022

CP, Plotting examples

Simple plots, defined in 2 functions, called as below:

plot_in_subplots(x,y1,y2,y3)
plot_same_frame(x, y1, y2, y3)   

The first function plots 3 subplots, the second 3 curves in the same plot.

@author: cp232
"""

import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-5, 5, 50)
y1 = x
y2 = x ** 2 - 4
y3 = x ** 3 + x**2 - 4*x


def plot_in_subplots(x_var, y_1, y_2, y_3):
    
    fig = plt.figure()
    fig.suptitle('3 Graphs separated', fontsize=16)

    fig.set_size_inches(15,5)

    ax = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    ax.plot(x_var, y_1)
    ax.grid()
    ax2.plot(x_var, y_2)
    ax2.grid()
    ax3.plot(x_var, y_3)
    ax3.grid()
    plt.show()
    
    
def plot_same_frame(x_var, y_1, y_2, y_3):
    plt.figure()
    
    plt.plot(x_var, y_1, 'black')
    plt.plot(x_var, y_2, 'r')
    plt.plot(x_var, y_3, 'b--')
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    
    plt.title('Plot title')
    
    plt.grid()
    plt.show()
    
