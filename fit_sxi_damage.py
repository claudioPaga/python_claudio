#! /usr/bin/python3

"""
fit_sci_damage

CTI parameters derivation fitting damaged spectrum

History:

V0- CP, 24 Nov 2024
    Replicates Andy Beardmore's fit_gauss.py, adapting it from a Gaussian to a
    CTI-damaging function 

"""


import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


from statistic import CStat

#Calculating the Gaussian PDF values given Gaussian parameters and random variable X
def gaus(X,C,X_mean,sigma):
    return C*np.exp(-(X-X_mean)**2/(2*sigma**2))  

def gauss(x, gc, gw, gn):
    """
    Gaussian : f(x) = gn * exp (-0.5 * z**2)

    where z = (x - gc) / gw

    with

       gc = gaussian centre
       gw = sigma
       gn = normalisation

    Integral is sqrt(2 * pi) * gn * gw
    """
    z = (x - gc) / gw
    model = gn * np.exp(-0.5 * z * z)
    return model


def eVDamage(en_trueKev, s_transfers, p_transfers, cti_alpha, cti_beta, p_trap_density_pixel, s_trap_density_pixel):
    """Perform the CTI damage.

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

    Apply the CTI damage based on CCD nY/pY and node output and CTI coeffs in
    in gain file.
    """
    import math
    import numpy as np

    store_transfers = 719     
    pcti_alpha = cti_alpha
    pcti_beta = cti_beta
    scti_alpha = cti_alpha
    scti_beta = cti_beta
    capture_prob = 1. - math.exp(-pcti_alpha * en_trueKev**(1.0-pcti_beta))
    losses = (p_transfers + store_transfers) * capture_prob * p_trap_density_pixel * 3.65/1000.
    damagedEnPKev =  (en_trueKev - losses)
            
    serial_capture_prob = 1. - math.exp(-scti_alpha * damagedEnPKev**(1.0-scti_beta))
    serial_losses = s_transfers * serial_capture_prob * s_trap_density_pixel * 3.65/1000.
    
    damagedEnSKev =  damagedEnPKev - serial_losses
    model = damagedEnSKev * 1000.
 
    return model



def fcn_andy(params, *args):
    """
    The function that is minimised - i.e. the CStat statistic.

    Args:
        params (list): the free parameters in the fit.
        args (tuple): the data (xdata and ydata).
    """
    # The data to fit passed in with args
    xd, yd = args
    # the free parameters to fit
    gc, gw, gn = params
    ym = gauss(xd, gc, gw, gn)
    cstat = CStat()
    stat = cstat(yd, ym)
    return stat


def fcn(params, *args):
    """
    The function that is minimised - i.e. the CStat statistic.

    Args:
        params (list): the free parameters in the fit.
        args (tuple): the data (xdata and ydata).
    """
    # The data to fit passed in with args
    xd, yd = args
    en_trueKev, s_transfers, p_transfers = xd
    cti_alpha, cti_beta, p_trap_density_pixel, s_trap_density_pixel = params
    # the free parameters to fit
    ym = np.zeros(len(en_trueKev))
    for i, en in enumerate(en_trueKev):
        ym[i] = eVDamage(en, s_transfers[i], p_transfers[i], cti_alpha, cti_beta, p_trap_density_pixel, s_trap_density_pixel)
    cstat = CStat()
    stat = cstat(yd, ym)
    print(stat)
    return stat



def fit_sxi_damage(xd, yd, *params):
    print('Input params:', params)
    result = minimize(fcn, params, args=(xd, yd), method='Nelder-Mead')
    return result


def fit_spec_andy(pha_bins, pha_spec, Eline=5895.45,
             lld=200, peak_width=300, ax=None):
    # bin centres.
    xc = (pha_bins[:-1] + pha_bins[1:]) * 0.5

    # Use the main peak maximum as the starting values,
    # but exclude the noise peak
    ind = xc > lld
    maxind = np.argmax(pha_spec[ind])
    gc = xc[ind][maxind]
    gn = pha_spec[ind][maxind]
    gw = 50.0

    # Select the pha range to fit over - use +/- 150 pha channels
    ind = (xc > gc - peak_width * 0.5) & (xc < gc + peak_width * 0.5)

    xd = xc[ind]
    yd = pha_spec[ind]

    # fit it
    res = fit_sxi_damage(xd, yd, gc, gw, gn)

    gain = Eline / res.x[0]
    fwhm_adu = res.x[1] * 2.355
    fwhm = fwhm_adu * gain

    print("Best fit parameters (gc, gw, gn):", res.x)
    print("Best fit cstat:", res.fun)
    print("Gain:", gain, "eV/DN")
    print("FWHM:", fwhm_adu, "ADU")
    print("FWHM:", fwhm, "eV")

    if ax is not None:
        # best fit model
        ym = gauss(xd, *tuple(res.x))
        ax.plot(xd, ym)

    return res


def fit_proc_call(evt_filename):
    
    """
    Procedure to set up the fitting, includes:
        Read in fits file
        Set starting CTI parameters
        Call minimization function
        
    Args:
        evt_filename (list): Damaged event file

    

    """
    import numpy as np
    from astropy.io import fits
    import math

    print("CTI calibration input evtfile =", evt_filename)

    
    # read necessary data columns of each event from the evtfile
    with fits.open(evt_filename, mode='update') as hdul:
        orig_cols = hdul['EVENTS'].data.columns
        ccdnr = hdul['EVENTS'].data.CCDNR
        amp = hdul['EVENTS'].data.AMP
        pha = hdul['EVENTS'].data.PHA
        pi = hdul['EVENTS'].data.PI
        rawx = hdul['EVENTS'].data.RAWX
        rawy = hdul['EVENTS'].data.RAWY
        mode = hdul['Primary'].header['DATAMODE']
        py_index = hdul['EVENTS'].data.CCDNR == 0
        ny_index = hdul['EVENTS'].data.CCDNR == 1
        obstime = hdul['Primary'].header['TSTART']

    # Scale CCD (X,Y) coordinates by the binning.
    binning = 1
    if mode == 'FT':
        binning = 6
    elif mode == 'FF':
        binning = 1
        
    # Spread events in the binned pixel
    rawxb = rawx*binning
    rawyb = rawy*binning
    

    # Initial guesses of CTI parameters, the same apart from their density for parallel and serial CTI
    p_trap_density_pixel = 0.03
    s_trap_density_pixel = 0.1
    cti_alpha_value = 0.15
    cti_beta_value = 0.2
    
    # Here select the Al, Mnka lines, and assign to xd the 'true' energy of 1.5 and 5.9 kev
    
    # Assumption: damaged spectra PI channel have a nominal gain of 10
    
    nominal_gain = 0.16
    pi = pi * nominal_gain
    
    # Select energies from the Al line at 1510 eV
    
    # 1 Initial gross selection
    al_index = np.logical_and(pi > 1000, pi < 1800)
    al_line = pi[al_index]
    
    # 2 Fit the energies using a Gaussian to determine centroid and width of the line
    mean = np.mean(al_line)
    sigma = np.std(al_line)
    counts, bins = np.histogram(al_line, bins = 100)
    hist=counts/sum(counts)
    n = len(hist)
    x_hist=np.zeros((n),dtype=float) 
    for ii in range(n):
        x_hist[ii]=(bins[ii+1]+bins[ii])/2   
    y_hist=hist
    #Gaussian least-square fitting process
    param_optimised,param_covariance_matrix = curve_fit(gauss,x_hist,y_hist,p0=[mean,sigma, max(y_hist)],maxfev=5000)
    print("Al line best fit pars:", param_optimised)
    
    # 3 Use the fit info to refine the Al line selection with energies within 3sigmas of Gaussian centroid value
    
    al_index = np.logical_and(pi > param_optimised[0]-3.0*param_optimised[1], pi < param_optimised[0]+3.0*param_optimised[1])
    al_line = pi[al_index]
    al_line_input = np.full_like(al_line, 1510)
    rawxb_al_line = rawxb[al_index]
    rawyb_al_line = rawyb[al_index]
    
    
    #########################
    # Select energies from the Mn line at 5900 eV
    
    # 1 Initial gross selection
    mn_index = np.logical_and(pi > 4000, pi < 6500)
    mn_line = pi[mn_index]
    
    # 2 Fit the energies using a Gaussian to determine centroid and width of the line
    mean = np.mean(mn_line)
    sigma = np.std(mn_line)
    counts, bins = np.histogram(mn_line, bins = 100)
    hist=counts/sum(counts)
    n = len(hist)
    x_hist=np.zeros((n),dtype=float) 
    for ii in range(n):
        x_hist[ii]=(bins[ii+1]+bins[ii])/2   
    y_hist=hist
    #Gaussian least-square fitting process
    param_optimised,param_covariance_matrix = curve_fit(gauss,x_hist,y_hist,p0=[mean,sigma, max(y_hist)],maxfev=5000)
    print("Mn line best fit pars:", param_optimised)

    # 3 Use the fit info to refine the Mn line selection with energies within 3sigmas of Gaussian centroid value
    
    mn_index = np.logical_and(pi > param_optimised[0]-3.0*param_optimised[1], pi < param_optimised[0]+param_optimised[1])
    mn_line = pi[mn_index]
    mn_line_input = np.full_like(mn_line, 5890)
    rawxb_mn_line = rawxb[mn_index]
    rawyb_mn_line = rawyb[mn_index]
    
    ###############
    # Concatenate the energies and coordinates of the two lines
    
    damaged_lines_en = np.concatenate((al_line, mn_line), axis = None)
    rawxb_lines = np.concatenate((rawxb_al_line, rawxb_mn_line), axis = None)
    rawyb_lines = np.concatenate((rawyb_al_line, rawyb_mn_line), axis = None)    
    reference_lines_en_kev = np.concatenate((al_line_input, mn_line_input), axis = None) / 1000.0

    #######################
    # Format the arrays for the fits
    # xd = damaged_energies (in KeV), rawxcoord, rawycoord
    # yd = reference energies (in eV)
    
    xd = reference_lines_en_kev, rawxb_lines, rawyb_lines
    yd = damaged_lines_en
    
    # fit it
    res = fit_sxi_damage(xd, yd, cti_alpha_value, cti_beta_value, p_trap_density_pixel, s_trap_density_pixel)

    print("Best fit parameters (cti_alpha_value, cti_beta_value, p_trap_density_pixel, s_trap_density_pixel):", res.x)
    print("Best fit cstat:", res.fun)

    
    return res

#fit_proc_call("/Users/cp232/SMILE/SXI_pipeline/TEST/FEV_Jan_ALL_A_OU_dmg.fits")
bfp = fit_proc_call("/Users/cp232/SMILE/SXI_pipeline/TEST/logo_Y1_test_dmg.evt")
print("Fitting procedure completed.")
