#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

A few examples with exceptions.

Goal of exceptions: - Program tries something, but in case of errors stops what it was doing and jumps to the par
of the code where I handle the exception.
There are different options on what to do after I handle the exception, I can do something if there were no exceptions (else: block)
or I can do something in any case (finally: block)

List of common errors:
    ValueError (Value is OK Type but not OK to do that operation, for example a division by zero)
    TypeError (trying to do something with incorrect type (abs(T)))
    NameError (variable name does not exist)
    IndexError (array index out of bound)
    SyntaxError
    AttributError (for classes)
    ZeroDivisionError
    
@author: cp232
"""

def generalException():
    
    """
    Returns a message in case of any exception

    @author: cp232

    """
    
    try:
        a = int(input("Please give me a number: "))
        b = int(input("Please give me another number: "))
        c = a/b
        print("All OK, the answer of the decision is ", c)
        # NOTE - With this structure the program will jump to except: for errors in any of the statements above
        # Code tries to go thru statements one line at a time but as soon as an error happens it jumps to except:
    except:
        print("The user gave incorrect input")
    else:
        print("Program finished")    
    
    #Code will always do this line, in case of exceptions or no exceptions    
    print("Program finished")    
    
def specificException():
    
    """
    Returns a message in case of a few specific exception

    """
    
    try:
        a = int(input("Please give me a number: "))
        b = int(input("Please give me another number: "))
        c = a/b
        print("All OK, the answer of the decision is ", c)
        # NOTE - With this structure the program will jump to except: for errors in any of the statements above
        # Code tries to go thru statements one line at a time but as soon as an error happens it jumps to except:
    except ValueError:
        print("User gave incorrect input Value")
    except TypeError:
        print("User gave incorrect type of input Value")    
    except ZeroDivisionError:
        print("Error, division by zero!")
    print("Program finished")    
       
    
def specificExceptionElseFinally():
    
    """
    Returns a message in case of a few specific exception

    """
    
    try:
        a = int(input("Please give me a number: "))
        b = int(input("Please give me another number: "))
        c = a/b
        print("All OK, the answer of the decision is ", c)
        # NOTE - With this structure the program will jump to except: for errors in any of the statements above
        # Code tries to go thru statements one line at a time but as soon as an error happens it jumps to except:
    except ValueError:
        print("User gave incorrect input Value")
    except TypeError:
        print("User gave incorrect type of input Value")    
    except ZeroDivisionError:
        print("Error, division by zero!")
    else:
        # Will do these operations only in case of no exception
        print("Well done user, the input values were OK")
    finally:
        print("I print this line in any case, exceptions or not")
        print("This finally block is mostly useful to close files that I've opened")
    print("Program finished")    
       

def checkOnInputValues():
     """
     Example of how to make sure that inputs are OK 
     """
     while True:
         #This is an infinite loop, will keep on running unless I break out of i
         try:
             a = int(input("Please enter a number: "))
             #The code will move on to the next line and break out of the loop if
             #no exceptions were raised up to now, in particular if the code was able to read in OK the input
             #value and to cast it to an int using int()
             break
             # The code will break out of the while loop and continue outside 
         except:
             print("You've entered a wrong input, try again!")
             
     print("You've entered an OK value")
        
     
def checkOnInputValuesWithAssert(L1, L2):
     """
     Input: Two Lists of the same length
     Example of check that input value is OK, if not just stop execution
     """
     assert len(L1)==len(L2), "Lists are of different lenght!!"
     Lconcat = L1+L2
     print(Lconcat)
     
 
        
     
     