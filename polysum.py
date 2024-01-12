#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:47:34 2019

@author: cp232
"""
def polysum(n, s):
    '''
    n: int, n>2, number of sides of a regolar polygon
    s: side length
    
    returns: Sum of the area and square of the perimeter, rounded to 4 decimal places
    '''
    
    import math 
    
    area = 0.25*n*s**2/math.tan(math.pi/n)
    perimeter = n*s
    
    return round(area + perimeter ** 2, 4)