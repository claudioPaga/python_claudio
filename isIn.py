#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 14:41:34 2019

@author: cp232
"""

def isIn(char, aStr):
    '''
    char: a single character
    aStr: an alphabetized string
    
    returns: True if char is in aStr; False otherwise
    '''
    # Your code here
    print(aStr)
    print(len(aStr))
    if len(aStr) == 0:
        return False;
    if len(aStr) == 1:
        print(char, aStr)
        if char == aStr:
            return True
        else:
            return False
    if len(aStr) > 1:
       middle =  len(aStr) // 2
       if char == aStr[middle]:
           return True
       else: 
           if char<aStr[middle]:
               return isIn(char, aStr[0:middle])  
           else:
               return isIn(char, aStr[middle+1:len(aStr)])
    