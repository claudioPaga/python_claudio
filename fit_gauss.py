import numpy as np
from scipy.optimize import minimize

from statistic import CStat


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


def fcn(params, *args):
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


def fit_gauss(xd, yd, *params):
    print('Input params:', params)
    result = minimize(fcn, params, args=(xd, yd), method='Nelder-Mead')
    return result


def fit_spec(pha_bins, pha_spec, Eline=5895.45,
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
    res = fit_gauss(xd, yd, gc, gw, gn)

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
