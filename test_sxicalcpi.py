#! /usr/bin/python3
"""Test sxicalcpi module - performs the gain and CTI correction"""

def test_sxicalcpi():
    """Calculates the PHA to PI gain/CTI correction

    INITIAL VERSION, CP, 22 Sept 2022.
    Test sxicalcpi using a prototype SXI event file,
    (created by Andy R, converted into fits by Andy B),
    and a simple Gain file.
    First check: assure that the corrected energies are
    higher than the measured ones.
    """
    import astropy
    from sxicalcpi import calcpi

    # Select event and gain files for testing.
    evt_filename = 'logo_Y1.evt'
    gain_filename = '/Users/cp232/SMILE/SXI_pipeline/CALDB/sxi_gain_v1.fits'

    # Read-in fits files, select relevant fields from 1st extension
    evt = astropy.io.fits.open(evt_filename)
    evt_data = evt[1].data
    pi = evt_data['PI']

    # Call sxicalcpi
    corrected_en = calcpi(evt_filename, gain_filename)

    # pyfits
    for i, measured_energies in enumerate(pi):
        bias = corrected_en[i]-measured_energies
        if abs(bias) < 0:
            print('WARNING, corrected energy lower than measured.')
            print('En_corr, En_measured')
            corrEnString = corrected_en[i]
            print('', corrected_en[i], ' ', measured_energies)
        else:
            print('En_corr, En_measured')
            print('', corrected_en[i], ' ', measured_energies)
            #print('' + corrEnString + ' ' + measured_energies)
                
    print('Test completed')
    # pyfits
    # pyfits
    # pyfits
    # pyfits
    # pyfits
    # pyfits
