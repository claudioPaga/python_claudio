#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:34:07 2019

Example of how scope rules work.
In increment definition, the n parameter is used as the variable to refer to the object n from the main practice14scope() variable
the n varialble will be incremented in increment() but not once back in practice14scope()
@author: cp232
"""

def increment(n): 
    n += 1 
    print(n) 
    #This will print 2


def practice14scope():
    n = 1 
    increment(n) 
    print(n) 
    # This will print 1, as n has not been updated in practice14scope