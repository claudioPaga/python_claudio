#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

@author: cp232
"""

def multiplication(a, b):
    
    """
    Input = 2 positive integers to multiply
    Return = a*b
    
    Created on Thu Aug 15 17:10:05 2019

    @author: cp232

    """

    #Returns the multiplication of a*b by iteration, given that a*b = a+a+a --- b times
    
    result = 0
    while b>0:
        result += a
        b -= 1
    return result    


def multiplicationIter(a, b):
    
    """
    Input = 2 positive integers to multiply
    Return = a*b
    
    Created on Thu Aug 15 17:10:05 2019

    @author: cp232

    """

    #Returns the multiplication of a*b by recursion.
    #This is based on the fact that a * b = a + a * (b-1) 
    #And a * (b-1) = a + a * (b-2)
    #and so on until b = 1, in that case a * 1 = a
    
    if b == 1:
        return a
    else:
        b -= 1
        print("b = ",b)
        print("a +  multiplicationIter(a, b) = ", a + multiplicationIter(a, b))
        return a + multiplicationIter(a, b)
     
def factorialRecursive(n):
    if n == 1:
        return 1
    else:
        print(n)
        return n * factorialRecursive(n-1)
    
    
def fibonacci(number):
     
    """
    Input = an integer number
    Return = Fibonacci of that number
    
    It's the number of rabbits after 5 years, each year the rabbits reproduce and it takes 1 year for the baby to be born
    
    The recursion has 2 base cases, with n = 1 and n = 2
    @author: cp232

    """
    
    if number == 1:
        return 1
        
    elif number == 2:
        return 2
    
    else:
        return fibonacci(number-1) + fibonacci(number-2)
    
    
def fibonacciEfficient(number, dictionary):
    """
    Input = an integer number
    Return = Fibonacci of that number
    
    It's the number of rabbits after 5 years, each year the rabbits reproduce and it takes 1 year for the baby to be born
    
    The recursion has 2 base cases, with n = 1 and n = 2
    
    The problem with the simple implementation is that it is slow, because it calculates many times the same thing.
    For example if n = 6 it needs to call fibonacci(5) and fibonacci(4), but fib(5) calls fib(4) and fib(3), and fib(4) calls fib(3) and fib(2).
    It keeps on calling stuff that it has already calculated.
    
    To avoid repetition we can use a dictionary, where the first time we calculate something we put it in the dictionary
    and the following times we just look it up in the dictionary.
    @author: cp232

    """
    
    
    if number in dictionary:
        # This will return the value of fibonacci already in the dictionary
        return dictionary[number]
    
    else:
        answer = fibonacciEfficient(number-1, dictionary) + fibonacciEfficient(number-2, dictionary)
        dictionary[number] = answer
        return answer
    
    