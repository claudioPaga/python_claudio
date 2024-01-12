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
n_time_entries = 5

met = np.zeros(n_time_entries)
met = [800451200.0, 803451200.0, 806451200.0, 809451200.0, 812451200.0]

trap_density_image = np.zeros(n_time_entries)
trap_density_serial = np.zeros(n_time_entries)

trap_density_image = [0.1, 0.1, 0.1, 0.1, 0.1]
trap_density_serial = [0.05, 0.05, 0.05, 0.05, 0.05]

# Gain coeff, for CCDpy/ny, node E/F, node readout 1/2

gain_ccdpy_e_2nodes = np.zeros(n_time_entries)
gain_ccdpy_f_2nodes = np.zeros(n_time_entries)
gain_ccdpy_e_1node = np.zeros(n_time_entries)
gain_ccdpy_f_1node = np.zeros(n_time_entries)
gain_ccdny_e_2nodes = np.zeros(n_time_entries)
gain_ccdny_f_2nodes = np.zeros(n_time_entries)
gain_ccdny_e_1node = np.zeros(n_time_entries)
gain_ccdny_f_1node = np.zeros(n_time_entries)

gain_ccdpy_e_2nodes = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdpy_f_2nodes = [1.64, 1.64, 1.64, 1.64, 1.64]
gain_ccdpy_e_1node = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdpy_f_1node = [1.64, 1.64, 1.64, 1.64, 1.64]
gain_ccdny_e_2nodes = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdny_f_2nodes = [1.64, 1.64, 1.64, 1.64, 1.64]
gain_ccdny_e_1node = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdny_f_1node = [1.64, 1.64, 1.64, 1.64, 1.64]


# CTI coeff, for CCDpy/ny, node E/F, node readout 1/2
# Likey overkill, but justified by possible CTI non-uniformity.

cti_alpha_p_ccdpy_e_2nodes = np.zeros(n_time_entries)
cti_alpha_p_ccdpy_f_2nodes = np.zeros(n_time_entries)
cti_alpha_p_ccdpy_e_1node  = np.zeros(n_time_entries)
cti_alpha_p_ccdpy_f_1node  = np.zeros(n_time_entries)
cti_alpha_p_ccdny_e_2nodes = np.zeros(n_time_entries)
cti_alpha_p_ccdny_f_2nodes = np.zeros(n_time_entries)
cti_alpha_p_ccdny_e_1node  = np.zeros(n_time_entries)
cti_alpha_p_ccdny_f_1node  = np.zeros(n_time_entries)

cti_beta_p_ccdpy_e_2nodes = np.zeros(n_time_entries)
cti_beta_p_ccdpy_f_2nodes = np.zeros(n_time_entries)
cti_beta_p_ccdpy_e_1node  = np.zeros(n_time_entries)
cti_beta_p_ccdpy_f_1node  = np.zeros(n_time_entries)
cti_beta_p_ccdny_e_2nodes = np.zeros(n_time_entries)
cti_beta_p_ccdny_f_2nodes = np.zeros(n_time_entries)
cti_beta_p_ccdny_e_1node  = np.zeros(n_time_entries)
cti_beta_p_ccdny_f_1node  = np.zeros(n_time_entries)

cti_alpha_s_ccdpy_e_2nodes = np.zeros(n_time_entries)
cti_alpha_s_ccdpy_f_2nodes = np.zeros(n_time_entries)
cti_alpha_s_ccdpy_e_1node  = np.zeros(n_time_entries)
cti_alpha_s_ccdpy_f_1node  = np.zeros(n_time_entries)
cti_alpha_s_ccdny_e_2nodes = np.zeros(n_time_entries)
cti_alpha_s_ccdny_f_2nodes = np.zeros(n_time_entries)
cti_alpha_s_ccdny_e_1node  = np.zeros(n_time_entries)
cti_alpha_s_ccdny_f_1node  = np.zeros(n_time_entries)

cti_beta_s_ccdpy_e_2nodes = np.zeros(n_time_entries)
cti_beta_s_ccdpy_f_2nodes = np.zeros(n_time_entries)
cti_beta_s_ccdpy_e_1node  = np.zeros(n_time_entries)
cti_beta_s_ccdpy_f_1node  = np.zeros(n_time_entries)
cti_beta_s_ccdny_e_2nodes = np.zeros(n_time_entries)
cti_beta_s_ccdny_f_2nodes = np.zeros(n_time_entries)
cti_beta_s_ccdny_e_1node  = np.zeros(n_time_entries)
cti_beta_s_ccdny_f_1node  = np.zeros(n_time_entries)


cti_alpha_p_ccdpy_e_2nodes = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdpy_f_2nodes = [0.25, 0.25, 0.25, 0.25, 0.25] 
cti_alpha_p_ccdpy_e_1node  = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdpy_f_1node  = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdny_e_2nodes = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdny_f_2nodes = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdny_e_1node  = [0.25, 0.25, 0.25, 0.25, 0.25]
cti_alpha_p_ccdny_f_1node  = [0.25, 0.25, 0.25, 0.25, 0.25]

cti_beta_p_ccdpy_e_2nodes = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdpy_f_2nodes = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdpy_e_1node  = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdpy_f_1node  = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_e_2nodes = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_f_2nodes = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_e_1node  = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_f_1node  = [0.2, 0.2, 0.2, 0.2, 0.2]

cti_alpha_s_ccdpy_e_2nodes = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdpy_f_2nodes = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdpy_e_1node  = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdpy_f_1node  = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdny_e_2nodes = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdny_f_2nodes = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdny_e_1node  = [0.3, 0.3, 0.3, 0.3, 0.3]
cti_alpha_s_ccdny_f_1node  = [0.3, 0.3, 0.3, 0.3, 0.3]

cti_beta_s_ccdpy_e_2nodes = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdpy_f_2nodes = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdpy_e_1node  = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdpy_f_1node  = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdny_e_2nodes = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdny_f_2nodes = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdny_e_1node  = [0.18, 0.18, 0.18, 0.18, 0.18]
cti_beta_s_ccdny_f_1node  = [0.18, 0.18, 0.18, 0.18, 0.18]


ttime = pyfits.Column(name='TIME', unit='s', format='D', array=met)
ttrap_density_image = pyfits.Column(name='IMA_DENSITY', format=form, array=trap_density_image)
trap_density_serial = pyfits.Column(name='SERIAL_DENSITY', format=form, array=trap_density_serial)


tg_ccdpy_e2n = pyfits.Column(name='G_PY_E2N', format=form, array=gain_ccdpy_e_2nodes)
tg_ccdpy_f2n = pyfits.Column(name='G_PY_F2N', format=form, array=gain_ccdpy_f_2nodes) 
tg_ccdpy_e1n = pyfits.Column(name='G_PY_E1N', format=form, array=gain_ccdpy_e_1node) 
tg_ccdpy_f1n = pyfits.Column(name='G_PY_F1N', format=form, array=gain_ccdpy_f_1node) 
tg_ccdny_e2n = pyfits.Column(name='G_NY_E2N', format=form, array=gain_ccdny_e_2nodes) 
tg_ccdny_f2n = pyfits.Column(name='G_NY_F2N', format=form, array=gain_ccdny_f_2nodes) 
tg_ccdny_e1n = pyfits.Column(name='G_NY_E1N', format=form, array=gain_ccdny_e_1node) 
tg_ccdny_f1n = pyfits.Column(name='G_NY_F1N', format=form, array=gain_ccdny_f_1node) 

tcti_alpha_p_ccdpy_e2n = pyfits.Column(name='CTI_ALPHA_P_PY_E2N', format=form, array=cti_alpha_p_ccdpy_e_2nodes)
tcti_alpha_p_ccdpy_f2n = pyfits.Column(name='CTI_ALPHA_P_PY_F2N', format=form, array=cti_alpha_p_ccdpy_f_2nodes) 
tcti_alpha_p_ccdpy_e1n = pyfits.Column(name='CTI_ALPHA_P_PY_E1N', format=form, array=cti_alpha_p_ccdpy_e_1node) 
tcti_alpha_p_ccdpy_f1n = pyfits.Column(name='CTI_ALPHA_P_PY_F1N', format=form, array=cti_alpha_p_ccdpy_f_1node) 
tcti_alpha_p_ccdny_e2n = pyfits.Column(name='CTI_ALPHA_P_NY_E2N', format=form, array=cti_alpha_p_ccdny_e_2nodes) 
tcti_alpha_p_ccdny_f2n = pyfits.Column(name='CTI_ALPHA_P_NY_F2N', format=form, array=cti_alpha_p_ccdny_f_2nodes) 
tcti_alpha_p_ccdny_e1n = pyfits.Column(name='CTI_ALPHA_P_NY_E1N', format=form, array=cti_alpha_p_ccdny_e_1node) 
tcti_alpha_p_ccdny_f1n = pyfits.Column(name='CTI_ALPHA_P_NY_F1N', format=form, array=cti_alpha_p_ccdny_f_1node) 

tcti_beta_p_ccdpy_e2n = pyfits.Column(name='CTI_BETA_P_PY_E2N', format=form, array=cti_beta_p_ccdpy_e_2nodes)
tcti_beta_p_ccdpy_f2n = pyfits.Column(name='CTI_BETA_P_PY_F2N', format=form, array=cti_beta_p_ccdpy_f_2nodes) 
tcti_beta_p_ccdpy_e1n = pyfits.Column(name='CTI_BETA_P_PY_E1N', format=form, array=cti_beta_p_ccdpy_e_1node) 
tcti_beta_p_ccdpy_f1n = pyfits.Column(name='CTI_BETA_P_PY_F1N', format=form, array=cti_beta_p_ccdpy_f_1node) 
tcti_beta_p_ccdny_e2n = pyfits.Column(name='CTI_BETA_P_NY_E2N', format=form, array=cti_beta_p_ccdny_e_2nodes) 
tcti_beta_p_ccdny_f2n = pyfits.Column(name='CTI_BETA_P_NY_F2N', format=form, array=cti_beta_p_ccdny_f_2nodes) 
tcti_beta_p_ccdny_e1n = pyfits.Column(name='CTI_BETA_P_NY_E1N', format=form, array=cti_beta_p_ccdny_e_1node) 
tcti_beta_p_ccdny_f1n = pyfits.Column(name='CTI_BETA_P_NY_F1N', format=form, array=cti_alpha_p_ccdny_f_1node) 

tcti_alpha_s_ccdpy_e2n = pyfits.Column(name='CTI_ALPHA_S_PY_E2N', format=form, array=cti_alpha_s_ccdpy_e_2nodes)
tcti_alpha_s_ccdpy_f2n = pyfits.Column(name='CTI_ALPHA_S_PY_F2N', format=form, array=cti_alpha_s_ccdpy_f_2nodes) 
tcti_alpha_s_ccdpy_e1n = pyfits.Column(name='CTI_ALPHA_S_PY_E1N', format=form, array=cti_alpha_s_ccdpy_e_1node) 
tcti_alpha_s_ccdpy_f1n = pyfits.Column(name='CTI_ALPHA_S_PY_F1N', format=form, array=cti_alpha_s_ccdpy_f_1node) 
tcti_alpha_s_ccdny_e2n = pyfits.Column(name='CTI_ALPHA_S_NY_E2N', format=form, array=cti_alpha_s_ccdny_e_2nodes) 
tcti_alpha_s_ccdny_f2n = pyfits.Column(name='CTI_ALPHA_S_NY_F2N', format=form, array=cti_alpha_s_ccdny_f_2nodes) 
tcti_alpha_s_ccdny_e1n = pyfits.Column(name='CTI_ALPHA_S_NY_E1N', format=form, array=cti_alpha_s_ccdny_e_1node) 
tcti_alpha_s_ccdny_f1n = pyfits.Column(name='CTI_ALPHA_S_NY_F1N', format=form, array=cti_alpha_s_ccdny_f_1node) 

tcti_beta_s_ccdpy_e2n = pyfits.Column(name='CTI_BETA_S_PY_E2N', format=form, array=cti_beta_s_ccdpy_e_2nodes)
tcti_beta_s_ccdpy_f2n = pyfits.Column(name='CTI_BETA_S_PY_F2N', format=form, array=cti_beta_s_ccdpy_f_2nodes) 
tcti_beta_s_ccdpy_e1n = pyfits.Column(name='CTI_BETA_S_PY_E1N', format=form, array=cti_beta_s_ccdpy_e_1node) 
tcti_beta_s_ccdpy_f1n = pyfits.Column(name='CTI_BETA_S_PY_F1N', format=form, array=cti_beta_s_ccdpy_f_1node) 
tcti_beta_s_ccdny_e2n = pyfits.Column(name='CTI_BETA_S_NY_E2N', format=form, array=cti_beta_s_ccdny_e_2nodes) 
tcti_beta_s_ccdny_f2n = pyfits.Column(name='CTI_BETA_S_NY_F2N', format=form, array=cti_beta_s_ccdny_f_2nodes) 
tcti_beta_s_ccdny_e1n = pyfits.Column(name='CTI_BETA_S_NY_E1N', format=form, array=cti_beta_s_ccdny_e_1node) 
tcti_beta_s_ccdny_f1n = pyfits.Column(name='CTI_BETA_S_NY_F1N', format=form, array=cti_alpha_s_ccdny_f_1node) 

coldefs = pyfits.ColDefs([ttime, 
                          ttrap_density_image,
                          trap_density_serial,
                          tg_ccdpy_e2n,
                          tg_ccdpy_f2n,
                          tg_ccdpy_e1n,
                          tg_ccdpy_f1n,
                          tg_ccdny_e2n,
                          tg_ccdny_f2n,
                          tg_ccdny_e1n,
                          tg_ccdny_f1n,
                          tcti_alpha_p_ccdpy_e2n,
                          tcti_alpha_p_ccdpy_f2n,
                          tcti_alpha_p_ccdpy_e1n,
                          tcti_alpha_p_ccdpy_f1n,
                          tcti_alpha_p_ccdny_e2n,
                          tcti_alpha_p_ccdny_f2n,
                          tcti_alpha_p_ccdny_e1n,
                          tcti_alpha_p_ccdny_f1n,
                          tcti_beta_p_ccdpy_e2n,
                          tcti_beta_p_ccdpy_f2n,
                          tcti_beta_p_ccdpy_e1n,
                          tcti_beta_p_ccdpy_f1n,
                          tcti_beta_p_ccdny_e2n,
                          tcti_beta_p_ccdny_f2n,
                          tcti_beta_p_ccdny_e1n,
                          tcti_beta_p_ccdny_f1n,
                          tcti_alpha_s_ccdpy_e2n,
                          tcti_alpha_s_ccdpy_f2n,
                          tcti_alpha_s_ccdpy_e1n,
                          tcti_alpha_s_ccdpy_f1n,
                          tcti_alpha_s_ccdny_e2n,
                          tcti_alpha_s_ccdny_f2n,
                          tcti_alpha_s_ccdny_e1n,
                          tcti_alpha_s_ccdny_f1n,
                          tcti_beta_s_ccdpy_e2n,
                          tcti_beta_s_ccdpy_f2n,
                          tcti_beta_s_ccdpy_e1n,
                          tcti_beta_s_ccdpy_f1n,
                          tcti_beta_s_ccdny_e2n,
                          tcti_beta_s_ccdny_f2n,
                          tcti_beta_s_ccdny_e1n,
                          tcti_beta_s_ccdny_f1n
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

fitsobj.writeto('sxi_gain_v2.fits')
