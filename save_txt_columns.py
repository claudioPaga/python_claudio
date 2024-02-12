#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

CP, Save columns in txt file

Example of how to save an ascii table with columns

@author: cp232
"""
    
import numpy as np
from astropy.io import ascii
from astropy.table import Table


def output_ascii_cols(col1, col2, col3):

    data = Table()
    data['X'] = col1
    data['Y'] = col2
    data['Z'] = col3
    ascii.write(data, 'example_ascii_col_file.dat', overwrite=True)


x = np.linspace(-5, 5, 50)
y = x ** 2 - 4
z = x ** 3 + x**2 - 4*x

output_ascii_cols(x,y,z)




    
