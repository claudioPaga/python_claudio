#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:41:06 2019

@author: cp232
"""

value = int(input("Please enter a number larger than 1 and I will find its square root: "))

low = 1
high = value
guess = (high+low)*0.5
tolerance = 0.001

print("initial guess = ", guess)
while abs(guess**2  -  value) >= tolerance:
    if guess**2 > value:
        # Guess was too large, I go lower
        high = guess
        guess = (high+low)*0.5
        print("Guess was too large, new guess = ", guess)
    else:
        # Guess was too small, I go higher
        low = guess
        guess = (high+low)*0.5
        print("Guess too small, new guess = ", guess)
           
if guess**2 == value:        
    print("Square root of ", value, " is ", guess)
else:
    print("I've found an approximate square of ", value, " being ", guess, " within a tolerance of ", tolerance)
    print(" Abs(value = ", value, "-", guess**2, ")<", tolerance)
       
    
     
    