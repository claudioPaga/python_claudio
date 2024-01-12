#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:30:51 2019

Examples of how to create numpy arrays and other numpy stuff

@author: cp232
"""

def practice16numpyArrays():
    
    import numpy as np
    import glob
    
    # Find all text files in this dir
    files_stats = glob.glob('*.txt')
    
    #Load data from table as a numpy array
    flux = np.loadtxt('/Users/cp232/SMILE/TestHarness/optical_loading/test_runs_pixel/flux1_frames2_T-100C/ima_evo_frn0.txt') 
    
    # First, create arrays of zeros
    array1dZeros = np.zeros(5)
    print(array1dZeros)
    # This will create 2D array (a matrix), with 3 = rows, 2 = columns
    matrix2dZeros = np.zeros((3,2))
    print(matrix2dZeros)
    
    # Now create and populate arrays
    
    arrayExample = np.array([2,3,3,4])
    print(arrayExample)
    array2dExample = np.array([[2,3,3], [1,1,1]])
    print(array2dExample)
    
    x = np.array([2, 4, 32])
    print(x[1]) 
    #Slicing
    print(x[0:2])
    slicedx = x[0:2]
    print('x = ', x)
    print('slicedx = ', slicedx)
    slicedx[0] = 1
    print('modified slicedx =', slicedx)
    print('x = ', x)
    print('Warning, slicing does not create a new array, it is just a slice of the original one')
    print('Changing the slice has changed the original one as well')
    index = [0,1]
    clonedx = x[index]
    print('clonedx =',clonedx)
    clonedx[0]=10
    print('Modified clonedx = ', clonedx)
    print('x = ', x)
    print('Warning! Using the index to define clonedx I cloned the original array. When I modify clonedx the original array is not modified')
    x2d = np.array([[1,3,5,2], [3, 4, 5, 1]])
    # This will print the first row
    print(x2d[0])
    # This will print the second row
    print(x2d[1,:])
    # This willprint the 2nd column
    print(x2d[:,1])
    # This will print the 1st and second column
    print(x2d[:,0:2])
    # Reverse
    print(x2d[::-1, ::-1])
    # Mask
    mask = arrayExample > 3
    print(arrayExample[mask])
    y = np.array([3,2, 1])
    x = x + y
    while x[0] < 100:
        x = x + y
        print(x)
        
    # Diffrence between slicing an array and indexing it
    # Slice = returns a view of the array
    # Index = returns a COPY 
    
    # Example of slicing, return a view, modifying the view will also modify the original !!!!!
    a = np.array([4, 5, 3, 2])
    b = a[0:2]
    b[0] = 1
    print(b)
    print(a)
    # Now 2D slicing
    a = np.array([[1.1, 2, 3], [1.5, 2.2, 3.1]])
    aSlicedFirstRows = a[0,:]
    print(aSlicedFirstRows)
    a = np.array([[1.1, 2, 3, 5, 4], [1.5, 2.2, 3.1, 4.4, 4.2]])
    aS = a[0:2:1, 1:4:2]
    print(aS) # Selected rows 0 and 1 (2 not included) and columns 1 and 3 (2 is the step)
    a[0,:] = [1, 2, 3, 3, 3] # Replaced row 0
    
    
    # Example of indexing, returns a copy, modifying the copy will NOT modify the original !!! (Very confusing)
    a = np.array([4, 5, 3, 2])
    index = [0,1]
    b = a[index]
    b[0] = 1
    print(b)
    print(a)
    
    # Example of concatenation
    a = np.array([[1.1, 2, 3], [1.5, 2.2, 3.1]])
    b = np.array([[0.1, 0.2, 0.3]])
    c = np.concatenate((a, b), axis=0)   
    # Create equally-spaced array
    x = np.linspace(0, 10, 11) #start, stop and number of elements
    print(x)
    x = np.log10(100)
    print(x)
    x = np.log(100)
    print(x)
    x = np.log2(100)
    print(x)
    
    #Create an ordered list of 20 elements, 0 to 19
    a = np.arange(20)
    # Check values
    b = np.all(a > 5) #Will return false
    b = np.any(a > 5) #will retrun true
    
    # Random numbers and plot:
    s = np.random.uniform(0, 1, 1000)
    import matplotlib.pyplot as plt
    count, bins, ignored = plt.hist(s, 15, density=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()
    
    
    rng = np.random.default_rng(42) # This to guarantee reproducibility
    rng.random(10)
    print(rng.integers(1,10,100)) # Will generate 100 integers between 1 and 9 
    
    
    
    
    
    
    
    