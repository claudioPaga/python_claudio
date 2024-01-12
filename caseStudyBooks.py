#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:53:42 2019

@author: cp232
"""

import os
import pandas as pd
import numpy as np
from collections import Counter

def count_words_fast(text):
    text = text.lower()
    skips = [".", ",", ";", ":", "'", '"', "\n", "!", "?", "(", ")"]
    for ch in skips:
        text = text.replace(ch, "")
    word_counts = Counter(text.split(" "))
    return word_counts

def read_book(title_path):
    text   = pd.read_csv(title_path, sep = "\n", engine='python', encoding="utf8")
    text = text.to_string(index = False)
    return text

def word_stats(word_counts):
    num_unique = len(word_counts)
    counts = word_counts.values()
    return (num_unique, counts)


hamlets = read_book('/Users/cp232/python/Harvard_tutorial/asset-v1_HarvardX+PH526x+2T2019+type@asset+block@hamlets.csv')
#language, text = hamlets.iloc[0]

counted_text = count_words_fast(hamlets)
data = pd.DataFrame(columns = ["word", "count"])
index = 1
#for values in counted_text.keys():
#    data.loc[index] = values, counted_text[values]
#    print(values, counted_text[values])
#    index +=1
#    print(index)
    
data

words =  counted_text.keys()
counts = counted_text.values()   
raw_data = {"word": words, "counts": counts}

data = pd.DataFrame(raw_data)

#data = pd.DataFrame(counted_text.keys(), counted_text.values(), columns = ["word", "count"])

