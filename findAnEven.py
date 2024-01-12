#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:04:01 2019

@author: cp232
"""

def findAnEven(L):
    """ Assumes L is a list of integers
        Returns the first even number in L
        Raises ValueError if L does not contain an even number"""
    foundEven = False
    for value in L:
        if value % 2 == 0:
            foundEven = True
            print(value)
            return(value)
    if foundEven == False:
        raise ValueError('There is no even number in the list')
            
    