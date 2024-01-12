#!/usr/bin/env python

import numpy as np
import astropy
import astropy.io.fits as pyfits
import datetime
import argparse
import sys
# import pylab as plt

description = "Create a FITs format SXI gain file.\n"


# file creation date (format is YYYY-MM-DDTHH:MM:SS)
creationdate = datetime.datetime.utcnow().replace(microsecond=0).isoformat()



#
# create a primary hdu
#
phdu = pyfits.PrimaryHDU()
phdu.header.set('TELESCOP', 'SMILE', 'Telescope (mission) name')
phdu.header.set('INSTRUME', 'SXI', 'Instrument name')

phdu_comments = [
    'FITS (Flexible Image Transport System) format is defined in \'Astronomy',
    'and Astrophysics\', volume 376, page 359; bibcode: 2001A&A...376..359H'
]

for comment in phdu_comments:
    phdu.header.add_comment(comment)

print("Primary HDU:")
for item in phdu.header.items():
    print(item)

#
# now the GAIN extension
#
# column format strings
form = str(1) + 'E'

# Initialise and fill the gain cti parameters, five time entries
met = np.zeros(5)
trap_density_image = np.zeros(5)
trap_density_serial = np.zeros(5)
cti_alpha_p = np.zeros(5)
cti_beta_p= np.zeros(5)
cti_apha_s = np.zeros(5)
cti_beta_s = np.zeros(5)
met = [800451200.0, 803451200.0, 806451200.0, 809451200.0, 812451200.0]
trap_density_image = [0.1, 0.1, 0.1, 0.1, 0.1]
trap_density_serial = [0.1, 0.1, 0.1, 0.1, 0.1]
cti_alpha_p = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_beta_p = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_alpha_s = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_beta_s = [0.2, 0.2, 0.2, 0.2, 0.2]

ttime = pyfits.Column(name='TIME', unit='s', format='D', array=met)
ttrap_density_image = pyfits.Column(name='IMA_DENSITY', format=form, array=trap_density_image)
trap_density_serial = pyfits.Column(name='SERIAL_DENSITY', format=form, array=trap_density_serial)
tcti_alpha_p = pyfits.Column(name='P_ALPHA', format=form, array=cti_alpha_p)
tcti_beta_p = pyfits.Column(name='P_BETA', format=form, array=cti_beta_p)
tcti_alpha_s = pyfits.Column(name='S_ALPHA', format=form, array=cti_alpha_s)
tcti_beta_s = pyfits.Column(name='S_BETA', format=form, array=cti_beta_s)


coldefs = pyfits.ColDefs([ttime,
                          ttrap_density_image,
                          trap_density_serial,
                          tcti_alpha_p,
                          tcti_beta_p,
                          tcti_alpha_s,
                          tcti_beta_s                          
                          ])

# pyfits
# ghdu = pyfits.new_table(coldefs)
# astropy.io.fits (or pyfits)
ghdu = pyfits.BinTableHDU.from_columns(coldefs)

ghdu.header.set('EXTNAME', 'GAIN',
                'Name of the binary table extension')

ghdu.header.set('TELESCOP', 'SMILE', 'Telescope (mission) name')
ghdu.header.set('INSTRUME', 'SXI', 'Instrument name')
ghdu.header.set('ORIGIN', 'LU', 'Source of FITS file')

# creator = str(pyfits.__package__) + ' ' + str(pyfits.__version__)
# astropy
creator = str(pyfits.__package__) + ' ' + str(astropy.__version__)
ghdu.header.set('CREATOR', creator, 'Creator')
ghdu.header.set('VERSION', 'v1', 'Extension version number')
ghdu.header.set('FILENAME', 'sxi_gain_v1.fits', 'File Name')

# mode
content = 'Photon Counting FF Gain'


ghdu.header.set('CONTENT', content, 'File content')
ghdu.header.set('TIMESYS', 'TT     ', 'Time system')
ghdu.header.set('MJDREFI', 51910, 'Reference MJD Integer part')
ghdu.header.set('MJDREFF', 7.4287037E-4, 'Reference MJD fractional')
ghdu.header.set('CLOCKAPP',
                False,
                'If clock corrections are applied (F/T)')

ghdu.header.set('CCLS0001', 'BCF',
                'Dataset is a Basic Calibration File')
ghdu.header.set('CDTP0001', 'DATA',
                'Calibration file contains data')
ghdu.header.set('CCNM0001', 'GAIN2', 'Type of Calibration data')

datamode = 'DATAMODE(' + 'PC' + ')'


ghdu.header.set('CDES0001',
                'SXI Gain CTI coefficients',
                'Description')

ghdu.header.set('NOM_GAIN', 10.,
                'Nominal Gain for PI channels')

gain_HDU_Comments = [
    "",
    "This extension provides CTI and charge trap coeff.s for the PI (pulse",
    "height invariant) energy scale calculation"
]

for comment in gain_HDU_Comments:
    ghdu.header.add_comment(comment)

ghdu.header.set('DATE', creationdate,
                'File creation date (YYYY-MM-DDTHH:MM:SS UT)')

for item in ghdu.header.items():
    print(item)

# Finally append the HDUs and create the FITs file.
fitsobj = pyfits.HDUList([phdu, ghdu])

fitsobj.info()

fitsobj.writeto('sxi_gain_v1.fits')
