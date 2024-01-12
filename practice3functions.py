#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:41:06 2019

@author: cp232
"""

def positive(x):
    
    """
    Returns true if number is positive
    """
    
    flag = False
    if (x>0):
        print(x, " is positive!")
        flag = True
    return flag        
# NOTE - If no return statement, python will return None value
# NOTE - When calling this, I can either use
# positive(10)
# positive(x = 10)   


def positive2(x, flag = False):
    
    """
    Returns true if number is positive
    By default, the flag is False, so I don't need the initial assignment
    """

    if (x>0):
        print(x, " is positive!")
        flag = True
    return flag        

