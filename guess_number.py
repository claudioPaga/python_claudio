#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:57:18 2019

@author: cp232
"""

print("Please think of a number between 0 and 100!")
guess = 50
low = 0
high = 100
print("Is your secret number ", guess, "?")
answer= input("Enter 'h' to indicate the guess is too high. Enter 'l' to indicate the guess is too low. Enter 'c' to indicate I guessed correctly. ")
while answer != 'c':
    while answer != 'l' and answer != 'h':
        print("Sorry, I did not understand your input.")
        print("Is your secret number ", guess, "?")
        answer= input("Enter 'h' to indicate the guess is too high. Enter 'l' to indicate the guess is too low. Enter 'c' to indicate I guessed correctly. ")
    if answer == 'l':
        low = guess
        guess = (high+low)//2
    else:
        high = guess
        guess = (high+low)//2
    print("Is your secret number ", guess, "?")
    answer= input("Enter 'h' to indicate the guess is too high. Enter 'l' to indicate the guess is too low. Enter 'c' to indicate I guessed correctly. ")

print("Game over. Your secret number was: ", guess)




