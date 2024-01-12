#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

Shows difference between copies (aliases) and clones of lists.

Aliases: Point to the same memory location, changing one changes the other
Clones: I am making a clone, a new memory location is creating. Changing the original will not change the clone.

@author: cp232
"""

def aliasesClones():
    
 
    
    listExample = [4, 3, 2, 1, 4, 3, 5, 4, 3]
    
    listCopy = listExample
    
    listExample.append(12)
    print('ListExample= ', listExample)
    print('ListCopy= ',listCopy)
    
    #NOTE the notation to create a CLONE!!
    listClone = listExample[:]
    listExample.append(111)
    print('ListExample= ', listExample)
    print('ListClone= ', listClone)
    
    # Example with sort and sorted
    
    #sorted() will NOT mutate the original array
    listExampleSorted = sorted(listExample)
    print('ListExample= ', listExample)
    print('ListExampleSorted= ', listExampleSorted)
    
    # NOTE - sort() will mutate the original array
    listExample.sort()
    print('ListExample after .sort() = ', listExample)
    
    
    
    
  