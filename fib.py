#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:14:43 2019

Calculates fibonacci number using recursive function + call to dictionary for efficiencty
n - number of generations
d - dictionary with base cases

Needs to be setup dictionary with the 2 base fibonacci cases (first generation = 1 birth, 2nd generation = 2 births)
    d = {1:1, 2:2}
    
returns fibonacci number

@author: cp232
"""

def fib(n,d):
    #Base cases first
    if n in d:
        return d[n]
    else:
        #Recursive calculation
        answer = fib(n-1, d) + fib(n-2, d)
        #Update dictionary with new lookup case
        d[n] = answer
        return answer
    