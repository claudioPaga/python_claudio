#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

Functions can be used as elements of a list.
Functions can be applied to every element of a list.

In general, this is formalised in Python using maps, a higher order procedure (HOP)

@author: cp232
"""

def applyFunctionToListExample(L, f):
    """
    Inputs:
        L is a list of numbers
        f is a function
    Outputs:
        Returns the list of f(L)
    """
    
    for counter in range(len(L)):
        L[counter] = f(L[counter])
        
    return(L)
    
    
def applyListOfFunctionsToValue(L, value):
    """
    Inputs:
        L is a list of function
        value is a number
    Outputs:
        Returns the list of L(value)
        
    """
    Lout = []
    for f in L:
        print(f(value))
        Lout.append(f(value))
    return(Lout)    
    

def setAbsMap(L):
     """
     Summary: this will do exactly as above but using the "map" HOP
     Inputs:
        L is a list of numbers
     Outputs:
         Returns the list of abs(L)      
     """
     Lout = []
     for element in map(abs, L):
         Lout.append(element)
         print(element)
     return(Lout)
     
     
def setMinMap(L1, L2):
     """
     Summary: Maps can be used to apply functions on multiple sets of input values
        L1, L2 are list of numbers is a list of function
       
     Outputs:
         Returns the list of min(L1, L2)      
     """
     Lout = []
     # NOTE - map produces an "interable", something that can be iterated over.
     # We need something like a "for loop" to access the elements of map.
     for element in map(min, L1, L2):
         Lout.append(element)
         print(element)
     return(Lout)
     
        
        
        
    
    
    
    
  