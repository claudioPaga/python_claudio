#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

Example of a dictionary

@author: cp232
"""

def studentsDictionaryPopulate(Lnames, Lgrades, Lage):
    """
    Inputs:
        Lgrades, Lnames Lage are lists  
    Outputs:
        A dictionary
    """

    studentsDictionary = {}
    
    # Populate the dictionary
    # Key will be the name of the student
    # Values will be a tuple, consisting of the grades and age of the student
    
    for counter in range(len(Lnames)): 
        value = (Lgrades[counter], Lage[counter])
        studentsDictionary[Lnames[counter]] = value
    
    # Print keys and values
    
    keys = studentsDictionary.keys()
    print(keys)
    for element in keys:
        print(element)
    values = studentsDictionary.values()
    print(values)
    for element in values:
        print(element)
    
    # Print out all the keys in the dictionary and the values of each key
    
    for students in studentsDictionary:
        print(students)
        print(studentsDictionary[students])
        
        
        
    
    return(studentsDictionary)    
        
    
    

        
        
        
    
    
    
    
  