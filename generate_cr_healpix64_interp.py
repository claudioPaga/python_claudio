#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:13:13 2019

Summary - ROSAT BG ESTIMATION PROCEDURES AND TESTING

Include functions to estimate the ROSAT X-ray diffuse background from pre-computed spectral fits of ROSAT data over an healpix map.
BG is calculated for the pixels of the map close to the input coordinate, these values are interpolated to estimate BG counts at input coords.

In generate_cr_healpix64_interp() the interpolation is done using the griddata() function.
In generate_cr_healpix64_interp_sphere() the interpolation is done using a function specific for points on an healpix map, but 
this requires the bg values over all the map to be computed before the interpolation.
The bg values over all pixels of the map are pre-computed using the function pre_compute_bg_counts_for_healpix_map()


@author: cp232
"""

def generate_cr_healpix64_interp(l_center, b_center, out_en_min, out_en_max, exposure_time, out_pixel_size_deg):
    """
    ;;; Summary - Determine Sky X-ray Background counts in SMILE (or XMM etc) pixels from ROSAT
    ;;;           Background model fits.
    ;;;           QDP model fits have been pre-computed for points on the
    ;;;           healpix grid (nside = 64) and archived in specific directory.  
    ;;;           The healpix grid points within a radius of the input
    ;;;           coordinates (l_center, b_center) are found using query_disc healpix
    ;;;           routine. Their correspondent QDP files are used to get the BG, and
    ;;;           the final estimated bg value is interpolated between
    ;;;           these points using the GRIDDATA routine.
    ;;;  
    ;;;
    ;;; Input - (l_center,b_center), Gal sky location 
    ;;;       - (out_en_min, out_en_max), energy range over which the bg is evaluated 
    ;;;       - exposure_time, the exposure time over which the bg counts are accumulated
    ;;;       - out_pixel_size_deg, the size of the "pixel" representing the FoV, in degrees.
    ;;;         For example, if out_pixel_size_deg = 1 degree, the FoV over which the BG is accumulated is a square of 1x1 degrees.
    ;;;
    ;;; Outputs - Text file with input pointings and estimated bg counts in output instrument FoV
    ;;; Returns - Estimated counts in FoV
    ;;;
    ;;; Procedure:  
    ;;; bg_cr_model (from ROSAT spectral fit)
    ;;; bg_cr_model/area (in sterad, ared of the ROSAT pointing)
    ;;; map new instrument (SMILE/XMM/etc) response to model spectrum
    ;;;
    ;;; NOTE - Spectral fits files from ROSAT script provide countrates
    ;;;        per arcmin (cts/s/arcimin) in ROSAT's 6 energy bins
    ;;;        So if I understood the units and what Steve needs correctly
    ;;;        to normalize the spectrum per steradiant the input size is
    ;;;        input_qdp_size_deg = 1 arcmin/60 = 0.016667
    ;;;        The 1 arcmin value is set by default if input keyword is not
    ;;;        provided explicetly
    ;;;
    ;;; EXAMPLE - Determine count rate at sky galactic coordinates (l, b)
    ;;;           = (280, -30), in the energy range [0.1, 2 keV], an
    ;;;           exposure of 300 seconds and an output pixel size of 2x2 degrees.
    ;;;           bg = generate_cr_healpix64_interp(280.0, -30.0, 0.1, 2.0, 300., 1.)
    ;;;           bg are the estimated background counts in 300 seconds in
    ;;;           the energy range [0.1-2] keV in the output instrument
    ;;;           FoV of 1x1 degrees
    ;;;
    ;;; HISTORY
    ;;; Adapted from equivalent IDL version
    ;;; generate_cr_point_healpix_interp_funct64.pro
    
    @author: cp232
"""

    
    import numpy as np
    import healpy as hp
    from astropy.io import ascii
    from scipy.interpolate import interp1d
    import sys
    from scipy.interpolate import griddata
    import math
    
    # Checks on the input parameters, should probably be done outside the function, before calling it.
    # Will figure this out later.
    
   
#    arguments = len(sys.argv) - 1
#    print ("the script is called with %i arguments" % (arguments))
#    print(sys.argv)
    
#    if arguments < 6:
#        print(len(sys.argv))
#        sys.exit(1)

    if isinstance(l_center, str):
        l_center = float(l_center)
        
    if isinstance(l_center, int) or isinstance(l_center, float):
        if l_center < 0 or l_center > 360:
            l_center = l_center % 360
    
    if not(isinstance(l_center, int) or isinstance(l_center, float)):   
        sys.exit(1)
    
    if isinstance(b_center, int) or isinstance(b_center, float):
        if b_center < -90 or b_center > 90:
            sys.exit(1)
    
    if not(isinstance(b_center, int) or isinstance(b_center, float)):   
        sys.exit(1)    
    
    if isinstance(out_en_min, int) or isinstance(out_en_min, float):
        if out_en_min < 0 or out_en_min > 10:
            out_en_min = 0
    
    if not(isinstance(out_en_min, int) or isinstance(out_en_min, float)):   
        sys.exit(1)    
    
    
    if isinstance(out_en_max, int) or isinstance(out_en_max, float):
        if out_en_max < 0:
            out_en_max = 10
    if out_en_max < out_en_min:
        sys.exit(1)   
                
    if not(isinstance(out_en_max, int) or isinstance(out_en_max, float)):   
        sys.exit(1)    
    
    if isinstance(exposure_time, int) or isinstance(exposure_time, float):
        if exposure_time < 0 :
            sys.exit(1)  
        
    if not(isinstance(exposure_time, int) or isinstance(exposure_time, float)):   
        sys.exit(1)           
        
       
    if isinstance(out_pixel_size_deg, int) or isinstance(out_pixel_size_deg, float):
        if out_pixel_size_deg < 0 :
            sys.exit(1)  
                
    if not(isinstance(out_pixel_size_deg, int) or isinstance(out_pixel_size_deg, float)):   
        sys.exit(1)           
        
       
    
    # The ROSAT qdp spectral files include counts in the ROSAT energy bins from an area of radius = 1 arcmin.
    input_qdp_size_deg = 0.01666
    # Determine area in sterad of ROSAT bg input models
    rosat_pointing_sqdeg = np.pi * (input_qdp_size_deg)**2
    rosat_pointing_sterad = rosat_pointing_sqdeg/3282.8 
    
    # Size of the "pixel" over which the background will be integrated.
    # For the estimate to make sense this should not be too different than the radius used to extract ROSAT bg spectra over the grid.
    # This should become an input value, with checks at the start, along with other parameters as the exposure time
    out_pixel_size_sterad = (out_pixel_size_deg)**2/3282.8   
      
    # Read SMILE response (EN, RESPONSE) from saved data file
    responseArea = ascii.read("/Users/cp232/python/mit_tutorial/sxi_resp.dat")
   
    ecommon = np.array(responseArea["col1"])
    response = np.array(responseArea["col2"])  

    # Use the healpix position to get the vector corresponding to the
    # input sky coordiantes
  
    vec = hp.ang2vec(l_center, b_center, lonlat=True)
    
    # Determine the pixels of the grid within a disk of radius radius_rad of the input coordinates (b_center, l_center) 
    NSIDE = 64
    nlist = 0
    radius_deg = 1.
    radius_rad = np.radians(radius_deg)
    while (nlist < 5):
        listpix = hp.query_disc(NSIDE, vec, radius_rad)
        nlist = len(listpix)
        radius_deg += 1
        radius_rad = np.radians(radius_deg)
    
    #Retrieve (l,b) cooridnates corresponding to listpix
    grid = ascii.read("/Users/cp232/python/mit_tutorial/healGrid64.dat")  
    lgrid = np.array(grid["l"])
    bgrid = np.array(grid["b"])
    l_disc = lgrid[listpix]
    b_disc = bgrid[listpix]
    sbgd_disc = np.zeros(nlist)

    # Get single healpix pixel corresponding to the input coordinates (l_center, b_center)
    ipring = hp.ang2pix(NSIDE, l_center, b_center, lonlat=True)

    # Get background estimates for pointings in the healpix grid in the proximity of the input location
    #print('l_disk  b_disk  bg_counts')
    for i, value in enumerate(listpix):
        qdp_file_name = '/Users/cp232/SMILE/ROSAT_bg/healpix/map64_r1/qdp_model_light/fit_bg_healpix' + str(value) + '_r1_light.qdp'
        qdp_table = ascii.read(qdp_file_name, data_start=3)
        qdp_energy = np.array(qdp_table["col1"])
        qdp_bin = np.array(qdp_table["col2"])
        qdp_model = np.array(qdp_table["col3"])
        qdp_model_sterad = qdp_model / rosat_pointing_sterad
        interpolatedResponseFunction = interp1d(ecommon, response)
        interpolatedResponse_at_qdp_energy = interpolatedResponseFunction(qdp_energy)
        sb_spec = qdp_model_sterad * interpolatedResponse_at_qdp_energy * exposure_time * out_pixel_size_sterad
        index = np.logical_and(qdp_energy >= out_en_min, qdp_energy <= out_en_max)
        sbgd_disc[i] = np.sum(sb_spec[index] * qdp_bin[index] * 2.) 
        if value == ipring:
            sbgd_pixel = sbgd_disc[i]
            
        #print('{:5.2f} {:7.2f} {:8.2f}'.format(l_disc[i], b_disc[i], sbgd_disc[i]))
     
    l_center_array = np.array(l_center)
    b_center_array = np.array(b_center)    
    coord_disc = (l_disc, b_disc)
    coord_disc_t = np.transpose(coord_disc)
    
    interpolatedBackground = griddata(coord_disc_t, sbgd_disc, (l_center_array, b_center_array), method='linear')
    # Try a different interpolation method is NaN is returned
    if math.isnan(interpolatedBackground):
        interpolatedBackground = griddata(coord_disc_t, sbgd_disc, (l_center_array, b_center_array), method='cubic')
        print(interpolatedBackground)
    if math.isnan(interpolatedBackground):
        interpolatedBackground = griddata(coord_disc_t, sbgd_disc, (l_center_array, b_center_array), method='nearest')
        print(interpolatedBackground)
    #If nothing works, simply return the background of the matching healpix pixel
    if math.isnan(interpolatedBackground):           
        interpolatedBackground = sbgd_pixel
        print(interpolatedBackground)
        
    # Output the estimated background to screen and file
#    print('')
#    print('Estimated background counts at (l, b) = ({:5.2f}, {:5.2f}): over a FoV area of {:1.6f} steradian.'.format(l_center, b_center, out_pixel_size_sterad))
#    print('Expected counts: {:.2f} '.format(interpolatedBackground))
#    print('Expected count rate in ' + str(exposure_time) + ' seconds: ' + str(interpolatedBackground/exposure_time))
#    print('')
#    
    # Output the estimated background counts to file
#    nameOutputFile = 'background_l' + str(l_center) + '_b'+ str(b_center) + '.dat'
#    outputFileHandle = open(nameOutputFile, 'w')
#    outputFileHandle.write('l   b   BackgroundCount    BackgroundCountRate' + '\n')
    #writestring = str(l_center) + ' ' + str(b_center)  + ' ' + str(interpolatedBackground) + ' ' + str(interpolatedBackground/exposure_time)
#    writestring = str(l_center) + ' ' + str(b_center)  + ' ' + str(interpolatedBackground) + ' ' + str(interpolatedBackground/exposure_time)
#    
#    outputFileHandle.write(writestring)
#    outputFileHandle.close()
#    table = {'l': (l_center), 'b': (b_center), 'BackgroundCount': (interpolatedBackground), 'BackgroundCountRate': (interpolatedBackground/exposure_time)}
#    ascii.write(table, nameOutputFile, formats={'l': '%.2f', 'b': '%.2f', 'BackgroundCount':'%.2f', 'BackgroundCountRate':'%.2f'})
#   
    estimatedBackground = [sbgd_pixel, interpolatedBackground]
    return estimatedBackground


def generate_cr_healpix64_interp_sphere(l_center, b_center, sbgd):
    
    """
    SUMMARY - Interpolates over a healpix map using the hp.interpolate_bilinear_lonlat() function
    
    RETURNS - Interplated bg value at the input coordinates.
    
    """
    
    from astropy_healpix import HEALPix
    from astropy.coordinates import Galactic
    import astropy.units as u
        
    # Use the healpix position to get the vector corresponding to the
    # input sky coordiantes  
#    vec = hp.ang2vec(l_center, b_center, lonlat=True)
    
    NSIDE = 64
    hp = HEALPix(nside=NSIDE, order='RING', frame=Galactic())
    #coord = SkyCoord(l_center, b_center, frame='galactic', unit = 'deg') 
    lon = l_center * u.deg
    lat = b_center * u.deg
    interpolatedBackground = hp.interpolate_bilinear_lonlat(lon, lat, sbgd)  
        
    return interpolatedBackground


def pre_compute_bg_counts_for_healpix_map(out_en_min, out_en_max, exposure_time, out_pixel_size_deg):
    
    """
    SUMMARY - Pre computes the bacgkround count rate on the healpix map for input energy range, exposure time and output FoV

    This will be used by the function:
        generate_cr_healpix64_interp_sphere(l_center, b_center, sbgd)
    That interplates over an healpix map
    
    RETURNS - the bg counts over the map
    
    """

    import numpy as np
    from astropy.io import ascii
    from scipy.interpolate import interp1d
    import healpy as hp
     # Read SMILE response (EN, RESPONSE) from saved data file
    responseArea = ascii.read("/Users/cp232/python/mit_tutorial/sxi_resp.dat")
   
    ecommon = np.array(responseArea["col1"])
    response = np.array(responseArea["col2"])
    interpolatedResponseFunction = interp1d(ecommon, response)

    # The ROSAT qdp spectral files include counts in the ROSAT energy bins from an area of radius = 1 arcmin.
    input_qdp_size_deg = 0.01666
    # Determine area in sterad of ROSAT bg input models
    rosat_pointing_sqdeg = np.pi * (input_qdp_size_deg)**2
    rosat_pointing_sterad = rosat_pointing_sqdeg/3282.8 
    
    # Size of the "pixel" over which the background will be integrated.
    # For the estimate to make sense this should not be too different than the radius used to extract ROSAT bg spectra over the grid.
    # This should become an input value, with checks at the start, along with other parameters as the exposure time
    out_pixel_size_sterad = (out_pixel_size_deg)**2/3282.8   

    # Pre-compute background counts for healpix map.
    NSIDE = 64
    NPIX = hp.nside2npix(NSIDE)
    print(NPIX)
    sbgd = np.zeros(NPIX)
    # Get background estimates for pointings in the healpix grid in the proximity of the input location
    #print('l_disk  b_disk  bg_counts')
    for i in range(NPIX):
        qdp_file_name = '/Users/cp232/SMILE/ROSAT_bg/healpix/map64_r1/qdp_model_light/fit_bg_healpix' + str(i) + '_r1_light.qdp'
        qdp_table = ascii.read(qdp_file_name, data_start=3)
        qdp_energy = np.array(qdp_table["col1"])
        qdp_bin = np.array(qdp_table["col2"])
        qdp_model = np.array(qdp_table["col3"])
        qdp_model_sterad = qdp_model / rosat_pointing_sterad      
        interpolatedResponse_at_qdp_energy = interpolatedResponseFunction(qdp_energy)
        sb_spec = qdp_model_sterad * interpolatedResponse_at_qdp_energy * exposure_time * out_pixel_size_sterad
        index = np.logical_and(qdp_energy >= out_en_min, qdp_energy <= out_en_max)
        sbgd[i] = np.sum(sb_spec[index] * qdp_bin[index] * 2.) 

    return sbgd


def test_bg_cr_input_table():
    """
    - SUMMARY -
    
    Test the procedures:
    generate_cr_healpix64_interp()
    generate_cr_healpix64_interp_sphere()
    
    for a set of random input coords from the input file:
    /Users/cp232/SMILE/ROSAT_bg/healpix/map64_r1/stats_healpix_distance_difference64_r1.txt
    """

    import numpy as np
    from astropy.io import ascii
        
    out_en_min = 0.1
    out_en_max = 2.0
    exposure_time = 300.
    out_pixel_size_deg  = 2.0

    # Read in input coordinates and IDL estimates
    # The values we want to compare against are the interpBg ones.
    bgTable = ascii.read("/Users/cp232/SMILE/ROSAT_bg/healpix/map64_r1/stats_healpix_distance_difference64_r1.txt")  
    linput = np.array(bgTable["l"])
    binput = np.array(bgTable["b"])
    pointBg = np.array(bgTable["bgPoint_r1"])
    interpBg = np.array(bgTable["bghealpixInterp"])
    
    # Pre-compute the background counts on the complete healpix map
    sbgd = pre_compute_bg_counts_for_healpix_map(out_en_min, out_en_max, exposure_time, out_pixel_size_deg)
    #l_center = 12
    #b_center = 44.45
    #bgInterSphere = generate_cr_healpix64_interp_sphere(l_center, b_center, sbgd)      
    
    # Setup the output file where the test results will be stored
    nameOutputTestFile = 'test_bg_idl_python.dat'
    outputFileHandle = open(nameOutputTestFile, 'w')
    outputFileHandle.write('l   b   BPGoint BGCountInterpIDL  BGCountInterpPython BGCountInterpSpherePython BGCountHealpixPixelPython ' + '\n')
    
    for index, value in enumerate(linput):
    # BG from interpolation on a plane + single healpix pixel
        bgPythonInterp = generate_cr_healpix64_interp(linput[index], binput[index], out_en_min, out_en_max, exposure_time, out_pixel_size_deg)
        # BG from interpolation on a sphere
        bgInterSphere = generate_cr_healpix64_interp_sphere(linput[index], binput[index], sbgd)
        print(index, interpBg[index], bgPythonInterp[1], bgPythonInterp[0], bgInterSphere)
        outputFileHandle.write('{:6.2f} {:6.2f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n'.format(float(linput[index]), float(binput[index]), float(pointBg[index]), float(interpBg[index]), float(bgPythonInterp[1]), float(str(bgInterSphere)), float(str(bgPythonInterp[0]))) )
        
    outputFileHandle.close()     

# Run the test
test_bg_cr_input_table()     