#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:39:04 2019

@author: cp232
"""

def moving_window_average(x, n_neighbors=1):
    import statistics
    n = len(x)
    width = n_neighbors*2 + 1
    x = [x[0]]*n_neighbors + x + [x[-1]]*n_neighbors
    print(width)
    print(x)
    # To complete the function,
    # return a list of the mean of values from i to i+width for all values i from 0 to n-1.
    averages = []
    for i in range(1, n+1):
        meanValue = statistics.mean(x[i-n_neighbors: i + n_neighbors+1])
        averages.append(meanValue)
    return averages

    