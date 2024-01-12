#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:41:06 2019

@author: cp232
"""

# Read input from user + casting (to integer)
x = input("Dammi un numero please ")
xint = int(x)
print(xint)
y = input("Dammi un altro numero please ")
yint = int(y)
print(yint)

# Operations
# Division
# Strings concatenation, 
division = xint/yint
# Print using strings concatenation
print("Division " + str(division))
output = "Division " + str(division)
print(output)
# Print without concatenation, note, in this case I don't need to convert division to a string in the print statement
print("Division ", division)
# 
# Integer Division
intx = xint//yint
print("Integer Division " + str(intx))
#Reminder
reminderx = xint%yint
print("Reminder " + str(reminderx))
#Power
power = xint**yint
print("Power " + str(power))

#
# Operators and nested conditionals
#

if (x < y):
    print(str(x) + " is < than " + str(y))
else:
    if (x > y):
        print(str(x) + " is > than " + str(y))

if (x == y):
     print(str(x) + " is equal to " + str(y))   

if (x != y):
     print(str(x) + " is NOT equal to " + str(y))   
  

if (x > y):
    print(str(x) + " is greater than " + str(y))
elif (x == 2):
   print(str(x) + " is equal to " + str(y)) 
else:
   print(str(x) + " is < than " + str(y))
        
# String 
print(len(output))   
print(output[0])
print(output[1:2])#NOTE- 2 not included!!!
print(output[1:4:2])#NOTE- 2 is the step 

#Values as reference
x = 3.45
y = x
x = 5.0
print("x="+str(x))
print("y="+str(y))