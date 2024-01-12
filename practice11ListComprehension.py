#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

List comprehension is a fast way of computing a list of operations (or better, an operation on a list)
It's elegant and it's fast

@author: cp232
"""

def squareListComprehension():
    
    """
    Returns a message in case of any exception

    @author: cp232

    """
    
    numbers = range(10)
    
    # Simple for loop to get squares
    
    squaresFor = []
    
    for n in numbers:
        square = n**2
        squaresFor.append(square)
        
    print(squaresFor)


    # List comprehnsion
    
    squaresComp = [n**2 for n in numbers]
    
    print(squaresComp)
    
    return(squaresComp)
     
 
        
     
     