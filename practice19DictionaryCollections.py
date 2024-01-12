#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:10:05 2019

Dictionary and collections

Goal - Count frequency of words in a text.
It can be done using a dictionary or using the collection Counter classs

The two implementations are below

@author: cp232
"""

def countWords(inputText):
    """
    Inputs:
        inputText, some text
    Outputs:
        A dictionary with keys = words, values = freqency of the key word in a text
    """

    wordDictionary = {}
    
    # Populate the dictionary
    # Key will be the name of the student
    # Values will be a tuple, consisting of the grades and age of the student
    
    # Step 1 - all lower case using the lower() method of the String class
    
    text = inputText.lower()
    
    # Step 2 - I want to remove most common punctuation
    
    punct = [".", ",","/",":",";","[","]","£","$","%","^","&","*","(",")","-","_","+","=","{","}","[","]","?",">","<","~",'"']
    
    # The replace() method will replace each occurance of a punctuation with a blank
    for bad in punct:
        text = text.replace(bad, "")
    print(text)
    # Step 3 - Count words in text, put them in dictionary
    # First, use split() to single out words in the text, then increase the counter of the dictionary if 
    # word already present in dictionary or initialize that key to 1.
    
    for word in text.split(" "):
        if word in wordDictionary:
            wordDictionary[word] +=1
        else:
           wordDictionary[word] = 1 
           
    return wordDictionary       



def countWordsCollections(inputText):
    """
    Inputs:
        inputText, some text
    Outputs:
        A dictionary with keys = words, values = freqency of the key word in a text
    """

    from collections import Counter
    wordDictionary = {}
    
    # Populate the dictionary
    # Key will be the name of the student
    # Values will be a tuple, consisting of the grades and age of the student
    
    # Step 1 - all lower case using the lower() method of the String class
    
    text = inputText.lower()
    
    # Step 2 - I want to remove most common punctuation
    
    punct = [".", ",","/",":",";","[","]","£","$","%","^","&","*","(",")","-","_","+","=","{","}","[","]","?",">","<","~",'"']
    
    # The replace() method will replace each occurance of a punctuation with a blank
    for bad in punct:
        text = text.replace(bad, "")
    print(text)
    # Step 3 - Count words in text, using the counter object of the collections class
    # it will take text.split(" "): as input
    
    wordDictionary = Counter(text.split(" "))
           
    return wordDictionary       



    
   
        
        
   
        
    
    

        
        
        
    
    
    
    
  