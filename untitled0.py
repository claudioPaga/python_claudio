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

def cti_iteration_simple_correction(scti_alpha, scti_beta, pcti_alpha, pcti_beta, s_trap_density_pixel, p_trap_density_pixel)

  ;;; Generate 100 random input energies between 0 and 10 keV and their X,Y locations
  energies_input_Kev = randomu(10, 100)*10.
  s_transfers = fix(randomu(1, 100)*1250+1)
  p_transfers = fix(randomu(1, 100)*3790+1)

  store_transfers = 719
  thresholdKev = 0.001
  iter_max = 10
  
  for i = 0, n_elements(energies_input_Kev ) - 1 do begin
     counter = 0
     losses_ev = 0.
     offsetkev = 10.
     en_iterationKev = energies_input_Kev[i]
     enMeasuredKeV = energies_input_Kev[i]
     damagedEnPKev = 0.
     damagedEnSKev = 0.
     while ((counter le iter_max) and (offsetKev gt thresholdKev)) do begin
        
        capture_prob = 1. - exp(-pcti_alpha * en_iterationKev^(1-pcti_beta))
        print, 'Pc_parallel = ', capture_prob
        ;;; Losses in KeV
        losses = (p_transfers[i] + store_transfers) * capture_prob * p_trap_density_pixel * 3.65/1000.
        damagedEnPKev =  (en_iterationKev - losses)

        cti_p = 1.-(damagedEnPKev/en_iterationKev)^(1./(p_transfers[i] + store_transfers))
            
        serial_capture_prob = 1. - exp(-scti_alpha * damagedEnPKev^(1-scti_beta))
        print, 'Pc_serial = ', serial_capture_prob
        ;;; Losses in KeV
        serial_losses = s_transfers[i] * serial_capture_prob * s_trap_density_pixel * 3.65/1000.
        damagedEnSKev =  damagedEnPKev - serial_losses

        cti_s = 1.-(damagedEnSKev/damagedEnPKev)^(1./s_transfers[i])
        
        offsetKev = (enMeasuredKeV - damagedEnSKev)
        counter = counter + 1
        print, 'enMeasuredKeV, en_iterationKeV, losses, serial_losses, damagedEnPKev, damagedEnSKev, offsetKeV (Measured-Damaged), New Correction guess'
        print, enMeasuredKeV*1000, en_iterationKeV*1000, losses*1000, serial_losses*1000, damagedEnPKev*1000, damagedEnSKev*1000, offsetKev*1000, (en_iterationKev + offsetKev)*1000

        en_iterationKeV = en_iterationKev + offsetKev
        
        
        stop

     endwhile

     print, 'Measured en, Corrected, Damaged, X, Y, CTIp, CTIs'
     print, energies_input_Kev[i]*1000, en_iterationKeV*1000, damagedEnPKev*1000, damagedEnSKev*1000, s_transfers[i], p_transfers[i], cti_p, cti_s

     stop
  endfor


end
