#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:04:06 2019

@author: cp232
"""

def f(n):
   """
   n: integer, n >= 0.
   """
   if n <= 1:
      return n
   else:
      return n * f(n-1)