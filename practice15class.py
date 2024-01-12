#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:48:55 2019

Simple example of class definition and class methods.

The new class is called Mylist.
I define it using class (NOT def, it's not a function!!!)
It inherits from the List class

@author: cp232
"""


class MyList(list):
    # class --> indicates this is a new class i'm defining here
    # MyList --> Type of the class. Instances will be objects of the myList class
    # List --> MyList inherits from the list class, will inherit its construction method, so it does not
    # need its constructor method
    
    # Class with two methods
    # First method of the class
    # self refers to the instance of the object on which the method is used.
    # For example when calling the method as y.remove_min(), self refers the instance y of the MyList class
    # the Method does not return anything, it simply modifies the object
    def remove_min(self):
        self.remove(min(self))                     
            


def practice15class():
    # As I havent defined a constructor for MyList it will inherith the List constructor
    L = MyList([1, 2, 3, 4])
    print(L)
    L.remove_min()
    print(L)
    
