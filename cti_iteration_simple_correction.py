#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 12:36:05 2022

Summary - Simple script to evaluate correctness of iterative CTI correction algorithm.
 
Capture probability expressed as
Pc = 1 - exp(-cti_alpha * En(KeV)^(1-cti_beta))

Losses expressed as
L = n_transfers * Pc * trap_density_per_pixel

Example:
cti_iteration_simple_correction(0.2, 0.2, 0.2, 0.2, 0.1, 0.1)

@author: cp232
"""

import numpy as np
import math

def cti_iteration_simple_correction(scti_alpha, scti_beta, pcti_alpha, pcti_beta, s_trap_density_pixel, p_trap_density_pixel):

 # Generate 100 random input energies between 0 and 10 keV and their X,Y locations
  energies_input_Kev = np.random.uniform(0, 10, 10)
  s_transfers = np.random.randint(1, 1250, 10)
  p_transfers = np.random.randint(1, 3790, 10)

  store_transfers = 719
  thresholdKev = 0.001
  iter_max = 10
  
  for i, energy_input in enumerate(energies_input_Kev):
     counter = 0
     offsetKev = 10.
     en_iterationKev = energy_input
     enMeasuredKeV = energy_input
     damagedEnPKev = 0.
     damagedEnSKev = 0.
          
     store_transfers = 719                           
    
     while (counter < iter_max) and (offsetKev > thresholdKev):
        
        capture_prob = 1. - math.exp(-pcti_alpha * en_iterationKev**(1.-pcti_beta))
        #print('Pc_parallel = ', capture_prob)
        losses = (p_transfers [i]+ store_transfers) * capture_prob * p_trap_density_pixel * 3.65/1000.
        damagedEnPKev =  (en_iterationKev - losses)
        cti_p = 1.-(damagedEnPKev/en_iterationKev)**(1./(p_transfers[i] + store_transfers))
            
        serial_capture_prob = 1. - math.exp(-scti_alpha * damagedEnPKev**(1.-scti_beta))
        #print('Pc_serial = ', serial_capture_prob)
        serial_losses = s_transfers[i] * serial_capture_prob * s_trap_density_pixel * 3.65/1000.
    
        damagedEnSKev =  damagedEnPKev - serial_losses
        cti_s = 1.-(damagedEnSKev/damagedEnPKev)**(1./s_transfers[i])
        offsetKev = (enMeasuredKeV - damagedEnSKev)
        #print('enMeasuredKeV, en_iterationKeV, losses, serial_losses, damagedEnPKev, damagedEnSKev, offsetKeV (Measured-Damaged), New Correction guess')
        #print(enMeasuredKeV*1000, en_iterationKev*1000, losses*1000, serial_losses*1000, damagedEnPKev*1000, damagedEnSKev*1000, offsetKev*1000, (en_iterationKev + offsetKev)*1000)

        en_iterationKev = en_iterationKev + offsetKev
        counter = counter + 1
        
     print('Measured en, Corrected, DamagedP, DamagedS, X, Y, CTIp, CTIs')
     print(energies_input_Kev[i]*1000, en_iterationKev*1000, damagedEnPKev*1000, damagedEnSKev*1000, s_transfers[i], p_transfers[i], cti_p, cti_s)
