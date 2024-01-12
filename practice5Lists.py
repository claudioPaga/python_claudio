#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

@author: cp232
"""

def listForLoop():
    
    # Lists are not really .
    # For loops with lists
    
    arrayExample = [4, 3, 2, 1, 4, 3, 5, 4, 3]
    
    total = 0
    
    for i in range(len(arrayExample)):
        total += arrayExample[i]
        
    print(total)

    
    total = 0
    
    for elementOfArray in arrayExample:
        # NOTE!!!  elementOfArray is not the index, it's the element of the array!
        total += elementOfArray
        
    print(total)
    
    
def listOperations(a):

    """
    Created on Thu Aug 15 17:10:05 2019

    Input = a list
    Output = the updated list after operations have been applied

    @author: cp232
    
    
    """ 

    print(a)
    a.append(12)
    print(a)
    a.append(14)
    b = a
    print(b)
    b.extend([12, 11])
    # With lists + is not the sum of the elements of the list, it's the concatenation!!!
    c = a + b
    
    # Other operations
    c.pop() # Removed last element
    c.remove(12) # Will remove the 1st element in the list that is equal to 12
    del(c[0]) # Will delete the first element of c
    print(c)
    
    # Sorting
    
    d = sorted(c) # NOTE, c is NOT sorted!!
    print(c)
    print(d)
    
    c.sort()
    print(c)
    
    c.reverse()
    print(c)
    
    f = 'addio!!'
    print(f)
    g = list(f)
    print(g)
    h = f.split('i')
    print(h)
    joinString = ''.join(g)
    print(joinString)
    joi = '__'.join(g)
    print(joi)
    return(c)
    