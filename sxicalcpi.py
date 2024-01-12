"""sxicalcpi module - performs the gain and CTI correction"""


def calcpi(evt_filename, gain_filename):
    """Calculates the PHA to PI gain/CTI correction

    Args:
        evtfile  (str): The input event file.
        gainfile (str): The input gain file (or CALDB).

    DRAFT VERSION, CP, 16 Sept 2022
    Implement an initial version of the iterative CTI corrections.

    NOTE: This basic version is a placeholder for things to be decided/developed in the future.
    It is intended as an example of a SXI tasks that can be used to understand what is needed and how to structure the software in the future.
    This draft version is based on 3 requirements:
    1 - Simple, fast, functional modelling of CTI losses.
    2 - Iterative procedure to determine corrected energy.
    3 - CTI modelling to retain a degree to physical representation of the traps, with parameters such as trap density and coefficients to model the capture probability.
    
    """

    import numpy as np
    import astropy
    import math
    
    print("calcpi: entering")

    print("calcpi:  evtfile =", evt_filename)
    print("calcpi: gainfile =", gain_filename)


    # Set placeholder values for iteration convergence threshold and number of iterations
    threshold = 1.0
    iter_max = 5
    
    # Read-in fits files, select relevant parameters

    evt = astropy.io.fits.open(evt_filename)
    evt_data = evt[1].data # Select data the 1st extension.
    
    pi = evt_data['PI']
    rawx = evt_data['RAWX'] # CCD coordinates, measured from the output node
    rawy = evt_data['RAWY']
    mode = evt[0].header['DATAMODE']

    # Scale CCD (X,Y) coordinates by the binning.
    binning = 1
    if mode == 'FT':
        binning = 1
    elif mode == 'FF':
        binning = 6
    rawx *= binning
    rawy *= binning

    gain =  astropy.io.fits.open(gain_filename)
    gain_coeff = gain[1].data # Select the first extension
    gain_caldb_time = gain_coeff['TIME']
    pcti_alpha_time = gain_coeff['P_ALPHA']
    pcti_beta_time = gain_coeff['P_BETA']
    scti_alpha_time = gain_coeff['S_ALPHA']
    scti_beta_time = gain_coeff['S_BETA']
    p_trap_density_time = gain_coeff['IMA_DENSITY']
    s_trap_density_time = gain_coeff['SERIAL_DENSITY']

    # Select the CTI coefficients matching the observing time using interpolation.
    obstime = evt[0].header['TSTART']
    pcti_alpha = np.interp(obstime, gain_caldb_time, pcti_alpha_time, left=None, right=None, period=None)
    pcti_beta = np.interp(obstime, gain_caldb_time, pcti_beta_time, left=None, right=None, period=None) 
    scti_alpha = np.interp(obstime, gain_caldb_time, scti_alpha_time, left=None, right=None, period=None)
    scti_beta = np.interp(obstime, gain_caldb_time, scti_beta_time, left=None, right=None, period=None)
    p_trap_density = np.interp(obstime, gain_caldb_time, p_trap_density_time, left=None, right=None, period=None)
    s_trap_density = np.interp(obstime, gain_caldb_time, s_trap_density_time, left=None, right=None, period=None)
           
    # Setup numpy array of corrected energies
    enCorrected = np.zeros(len(pi))

    # Derive CTI-corrected energies
    
    for i, enMeasured in enumerate(pi):
        enCorrected[i] = iterCorrection(enMeasured/1000., rawx[i], rawy[i], pcti_alpha, pcti_beta, scti_alpha, scti_beta, p_trap_density, s_trap_density, threshold, iter_max)
        
    print("calcpi: exiting")

    return enCorrected


def iterCorrection(en_measuredKev, s_transfers, p_transfers, pcti_alpha, pcti_beta, scti_alpha, scti_beta, p_trap_density_pixel, s_trap_density_pixel, threshold, iter_max):
    """
    Args:
    en_measuredKeV: measured energy, in KeV
    s_transfers: unbinned serial readout transfers
    p_transfers: unbinned parallel readout transfers
    pcti_alpha: first coeff of parallel capture probability function
    pcti_beta: second coeff of parallel capture probability function
    scti_alpha: first coeff of serial capture probability function
    scti_beta: second coeff of serial capture probability function
    p_trap_density_pixel: trap density per pixel in image and store section
    s_trap_density_pixel: trap density per pixel in serial register
    threshold: convegence threshold of iterative correction algorithm
    
    ITERATIVE CTI CORRECTION
    
    For details, see https://space.mit.edu/ACIS/cticor.pdf
    Townsley et al 2000, ApJ 534 L139 
    Refines the estimate of the corrected energy by converging on
    the damaged energy comparing the measured energy with the
    CTI-processed energy.
    The initial guess is the measured energy, to which a delta
    equal to the difference between measured and cti-damaged energy
    is added at each refined step.
    D = Damaged
    M = Measured
    C = Corrected
    O = Offsets, M-D
    T = Threshold  
 
    Hp0: C = M
    C --> CTI --> D
    O = M-D is O < T? IF so, we've found C, otherwise
   
    Hp1 C = C + O
    C --> CTI --> D
    O = M-D is O < T? IF so, we've found C, otherwise

    Hp2 C C + 0  
    C --> CTI --> D
    O = M-D
    And so on until we get a convergence and O<T  

    """
    import math
    # Initialise iterative method counter and losses
    counter = 0
    offsetKev = 10.
    en_iterationKev = en_measuredKev
    thresholdKev = threshold/1000.

    store_transfers = 719                           
    
    while (counter < iter_max) and (offsetKev > thresholdKev):
        
        capture_prob = 1. - math.exp(-pcti_alpha * en_iterationKev**(1.0-pcti_beta))
        losses = (p_transfers + store_transfers) * capture_prob * p_trap_density_pixel * 3.65/1000.
        damagedEnPKev =  (en_iterationKev - losses)
            
        serial_capture_prob = 1. - math.exp(-scti_alpha * damagedEnPKev**(1.0-scti_beta))
        serial_losses = s_transfers * serial_capture_prob * s_trap_density_pixel * 3.65/1000.
    
        damagedEnSKev =  damagedEnPKev - serial_losses

        offsetKev = (en_measuredKev - damagedEnSKev)
        counter = counter + 1
        en_iterationKev = en_iterationKev + offsetKev
        
    en_iteration = en_iterationKev * 1000.     
    return en_iteration

