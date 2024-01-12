#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:13:13 2019

@author: cp232
"""

def healpix_map_creation():
    import numpy as np
    import healpy as hp
    from astropy.io import ascii
    NSIDE = 32
    NPIX = hp.nside2npix(NSIDE)
    print(NPIX)
    print(360*180)
    m = np.arange(NPIX)
    hp.mollview(m, title="Mollview image RING")
    hp.graticule()
    
    pointings = hp.pix2ang(NSIDE, m, lonlat = True)
    print(pointings) 
    
    table = {'HealPixels': m, 'l': pointings[0], 'b': pointings[1]}
    
    ascii.write(table, 'healGrid.dat', formats={'l': '%.2f', 'b': '%.2f'})
    
    
    NSIDE = 64
    NPIX = hp.nside2npix(NSIDE)
    m = np.arange(NPIX)
    hp.mollview(m, title="Mollview image RING")
    hp.graticule()
    
    pointings = hp.pix2ang(NSIDE, m, lonlat = True)
    print(pointings) 
    
    table = {'HealPixels': m, 'l': pointings[0], 'b': pointings[1]}
    
    ascii.write(table, 'healGrid64.dat', formats={'l': '%.2f', 'b': '%.2f'})


NSIDE = 64
NPIX = hp.nside2npix(NSIDE)
m = np.arange(NPIX)
hp.mollview(m, title="Mollview image RING")
hp.graticule()

pointings = hp.pix2ang(NSIDE, m, lonlat = True)
print(pointings) 

table = {'HealPixels': m, 'l': pointings[0], 'b': pointings[1]}

ascii.write(table, 'healGrid64.dat', formats={'l': '%.2f', 'b': '%.2f'})

