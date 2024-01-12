#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:20:18 2019

@author: cp232
"""

def biggest(aDict):
    '''
    aDict: A dictionary, where all the values are lists.

    returns: The key with the largest number of values associated with it
    '''
    # Your Code Here
    lengthValues = {}
    for elements in aDict:
        value = aDict[elements]
        lengthValues[elements] = len(value)

    length = lengthValues.values()
    mostFrequent = max(length)

    keysFrequent = None
    for keys in lengthValues:
        if lengthValues[keys] == mostFrequent:
            keysFrequent = keys
            
    return keysFrequent