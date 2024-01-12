#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 08:55:59 2019

Simple function to illustrate what happens to parameters passed to a function in python.
 

@author: cp232
"""

def functionParametersExample(par1Integer, par2List):
    
    """
    Input:
        par1Integer is an integer variable
        par2List is a List
    Output: 
        outTuple tuple, returning input parameters after modification if function 

    @author: cp232

    """
    global c
    par1Integer +=1
    c += 1
    # NOTE - a will be modified locally, but not globally! (the same happens in Java I believe)
    # If I try to run this:
    # a = 5
    # b = [1, 2, 3, 4]
    # rtuple = functionParametersExample(a, b)
    # rtuple
    # (6, [2, 2, 3, 4])
    # a
    # 5
    # b 
    #[2, 2, 3, 4]

    # NOTE - Doing this will modify b, because b is passed as an argument as reference.
    # Modifying it in the function will modify the values of the list, so b will be modified outside the function as well.
    par2List[0] = par2List[0]+1
    
    outTuple = (par1Integer, par2List, c)
    import string
    dic = {"aa":1, "bb":2}
    
    return(outTuple)    
    