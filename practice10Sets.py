#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

A few examples of set operations.

Sets are collections of UNIQUE values
Sets in python have "normal" sets operations, union, etc..

@author: cp232
"""

def setOperations():
    
    """
    Returns a message in case of any exception

    @author: cp232

    """
    
    ids = set([1, 3, 5, 14, 33, 7])
    
    # other definitio of a set
    
    ids2 = {4, 5, 555}
    
    females = set([1, 14])
    
    males = ids - females
    
    everybody = females | males
    
    print(everybody)
    
    both = females & males
    
    print(both)
    
  
     
 
        
     
     