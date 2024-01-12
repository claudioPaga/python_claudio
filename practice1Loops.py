#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:41:06 2019

@author: cp232
"""

n = 0
while n<5:
    print(n)
    n += 1

for n in range(10):
    print(n)

# range(min, max, step)
for n in range(10, 40, 5):
    print(n)
    
# Break out of loops
 
sum = 0    
for n in range(10, 40, 5):
    sum +=5
    print(sum)
    if sum == 30:
        break

       
    
     
    