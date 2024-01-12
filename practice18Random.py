#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:05:15 2019

@author: cp232
"""

def practice18Random():
    import random
    import numpy as np
    import matplotlib.pyplot as plt

    x = random.choice(range(100))
    print(x)
    # Will pick a sample of UNIQUE values (without replacement) from the input list (range(100))
    x = random.sample(range(100), 10)
    plt.figure()
    plt.hist(x)
    plt.show
    
    # Simulate random walk using Numpy
    
    # Starting point, a numpy array of dimension 2x1
    x0 = np.array([[0],[0]])
    # 1000000 Random steps of lenght from the normal distribution in 2 dimensions
    delta_steps = np.random.normal(0,1, [2, 1000000])
    # Accumulate the steps.
    # The input is a 2D array, I need to specify over which axis I want to sum (1 = sum the columns)
    delta_steps = np.cumsum(delta_steps, 1)
    # I concatenate the steps with the starting point of (0,0)
    X = np.concatenate((x0, delta_steps), axis=1 )
    plt.figure()
    plt.plot(X[0], X[1], "yo-")
    plt.savefig("randomWalk.png")
    