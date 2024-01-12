#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 09:24:46 2022

@author: cp232

CATEGORY: CTI, CTIMonitor

INPUTS:
table_corrections = Table including CalUnits info and CTI offsets as a function of AC

EXAMPLE:
ctimonitor_corr_calunits_stats('/Users/cp232/Gaia/CTIMonitor/EmpiricalSCTICorrections/AcLocationBiasByMuAndSignal-IPD32-4P-NoW_Rev1500_5400_I50_poly4.asc')    
    

HISTORY:
    - CP, 1 Sept 2022
    Initial version

"""

import numpy as np


def ctimonitor_corr_calunits_stats(table_corrections):
    
    # Set output  filenames 
    namesplit = table_corrections.split("poly4")
    endNameSplit = namesplit[1].split(".")
    outStatsName = namesplit[0]+'poly4'+endNameSplit[0]+'_CalUStats.txt'

    # Read in table, converting fields to integer
    data = np.loadtxt(table_corrections) 
    
    # Select high stats bins only
    
    numInBin = data[:,13].astype(int)
    highStatIndex = np.where(numInBin > 100)
    
    # .astype(int) converts numpy array to integer
    # data[highStatIndex,0].astype(int) will be a [1xN] array, I need to flatten it to a 1D array
    # ravel() will do the 1D flattening
    
    ccdRow = data[highStatIndex,0].astype(int).ravel()
    ccdStrip = data[highStatIndex,1].astype(int).ravel()
    wc = data[highStatIndex,2].astype(int).ravel()
    gate = data[highStatIndex,3].astype(int).ravel()
    fov = data[highStatIndex,4].astype(int).ravel()
    calUnitId = data[highStatIndex,5].astype(int).ravel()
    obmtRev = data[highStatIndex,6].ravel()
    obmtRevMin = data[highStatIndex,7].ravel()
    obmtRevMax = data[highStatIndex,8].ravel()
    acPos = data[highStatIndex,9].astype(int).ravel()
    signal = data[highStatIndex,10].ravel()
    clippedMeanAcLocBias = data[highStatIndex,11].ravel()
    acLocBiasSig = data[highStatIndex,12].ravel()
    
    
    gates = [0, 4, 7, 8, 9, 10, 11, 12]
    countPrintLines = 0
    print('CalU, CalU_Entries, MeanBias, Median, STD')
    with open(outStatsName, 'w') as f:
         f.write('#CalU CalU_Entries MeanCTIBias MedianCTIBias CtiStd')
         f.write('\n') 
         f.close()    
         
    #Go over ccdRow, Strip, FoV and Gate selection to select a specific CalUnit     
    for ccdRowN in range(1, 8):
        for ccdStripN in range(1, max(ccdStrip)+1):
            for gateN in gates:
                for fovN in range(1,3):
                    #calU_index = np.where(ccdRow == ccdRowN and ccdStrip == ccdStripN and gate == gateN and fov == fovN)
                    
                    calU_index = np.where(ccdRow == ccdRowN)
                    if (calU_index[0]).shape[0] > 0:
                        ccdStripSR = ccdStrip[calU_index]
                        gateSR = gate[calU_index]
                        fovSR = fov[calU_index]
                        clippedMeanAcLocBiasSR = clippedMeanAcLocBias[calU_index]
                        calUnitIdSR = calUnitId[calU_index]
                        
                        # Search over strip
                        calU_indexStrip = np.where(ccdStripSR == ccdStripN)
                        if (calU_indexStrip[0]).shape[0] > 0:
                            gateSS = gateSR[calU_indexStrip]
                            fovSS = fovSR[calU_indexStrip]
                            clippedMeanAcLocBiasSS = clippedMeanAcLocBiasSR[calU_indexStrip]
                            calUnitIdSS = calUnitIdSR[calU_indexStrip]
                            
                            # Search over gates
                            calU_indexGate = np.where(gateSS == gateN)
                            if (calU_indexGate[0]).shape[0] > 0:
                                fovSG = fovSS[calU_indexGate]
                                clippedMeanAcLocBiasSG = clippedMeanAcLocBiasSS[calU_indexGate]
                                calUnitIdSG = calUnitIdSS[calU_indexGate]
                                
                                # Search over FOV
                                calU_indexFov = np.where(fovSG == fovN)
                                if (calU_indexFov[0]).shape[0] > 0:
                                    clippedMeanAcLocBiasSF = clippedMeanAcLocBiasSG[calU_indexFov]
                                    calUnitIdSF = calUnitIdSG[calU_indexFov]
    
                                    # Calculate statistics
                                    meanCu = clippedMeanAcLocBiasSF.mean()
                                    medianCu = np.median(clippedMeanAcLocBiasSF)
                                    stDevCu = clippedMeanAcLocBiasSF.std()
                                    
                                        
                                    #print('CalU = ' + str(calUnitIdSF[0]) + ' CalUN ' + str((calU_indexFov[0]).shape[0]))
                                    print(str(calUnitIdSF[0]), str((calU_indexFov[0]).shape[0]), meanCu, medianCu, stDevCu) 
                                    out_txt = '' + str(calUnitIdSF[0]) + ' ' + str((calU_indexFov[0]).shape[0]) + ' ' + str(meanCu) + ' ' + str(medianCu) + ' ' + str(stDevCu)
                                        
                                    # Output stats                                     
                                    with open(outStatsName, 'a') as f:
                                        f.write(out_txt)
                                        f.write('\n') 
                                        f.close()    
                                        
    


