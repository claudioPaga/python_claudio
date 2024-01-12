#!/usr/bin/env python
"""
Summary: Generate Gain file for SXI
Description:
17/4/2023 - Updated CTI pars based on OU calibration doc

2/12/2022 - Drop distinction between 1 and 2 nodes readout
as we will only change between the 2 rarely.
If this happens there will be entries in the gain file
at that time, preferably two very close in time
to ensure the time interpolation in calcpi
will be correct.
Date: 17 April 2023
Author: CP
"""

import numpy as np
import astropy
import astropy.io.fits as pyfits
import datetime
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

trap_density_image_py = np.zeros(n_time_entries)
trap_density_serial_py = np.zeros(n_time_entries)
trap_density_image_ny = np.zeros(n_time_entries)
trap_density_serial_ny = np.zeros(n_time_entries)

# Trap density values slightly different just to differentiate
# between the two halfes of the CCD.
# Parameters are initial guesses based on OU calib campaign
# /Users/cp232/Documents/SMILE/Open_SMILE_CAL_TR_01_01_00 FMCCD370_Calibration_Test_Report.pdf
# Using the script
# /Users/cp232/SMILE/TestRadHardness/cti_capture_p_simple_function.pro

trap_density_image_py = [0.003, 0.003, 0.003, 0.003, 0.003]
trap_density_serial_py = [0.01, 0.01, 0.01, 0.01, 0.01]
trap_density_image_ny = trap_density_image_py
trap_density_serial_ny = trap_density_serial_py
#trap_density_image_ny = [0.0032, 0.0032, 0.0032, 0.0032, 0.0032]
#trap_density_serial_ny = [0.011, 0.011, 0.011, 0.011, 0.011]

# Gain coeff, for CCDpy/ny, node E/F, node readout 1/2

gain_ccdpy_e = np.zeros(n_time_entries)
gain_ccdpy_f = np.zeros(n_time_entries)
gain_ccdny_e = np.zeros(n_time_entries)
gain_ccdny_f = np.zeros(n_time_entries)

gain_ccdpy_e = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdpy_f = gain_ccdpy_e
# gain_ccdpy_f = [1.64, 1.64, 1.64, 1.64, 1.64]
gain_ccdny_e = [1.6, 1.6, 1.6, 1.6, 1.6]
gain_ccdny_f = gain_ccdny_e
# gain_ccdny_f = [1.64, 1.64, 1.64, 1.64, 1.64]


# CTI coeff, for CCDpy/ny, node E/F, node readout 1/2
# Likey overkill, but justified by possible CTI non-uniformity.

cti_alpha_p_ccdpy_e = np.zeros(n_time_entries)
cti_alpha_p_ccdpy_f = np.zeros(n_time_entries)
cti_alpha_p_ccdny_e = np.zeros(n_time_entries)
cti_alpha_p_ccdny_f = np.zeros(n_time_entries)

cti_beta_p_ccdpy_e = np.zeros(n_time_entries)
cti_beta_p_ccdpy_f = np.zeros(n_time_entries)
cti_beta_p_ccdny_e = np.zeros(n_time_entries)
cti_beta_p_ccdny_f = np.zeros(n_time_entries)

cti_alpha_s_ccdpy_e= np.zeros(n_time_entries)
cti_alpha_s_ccdpy_f= np.zeros(n_time_entries)
cti_alpha_s_ccdny_e= np.zeros(n_time_entries)
cti_alpha_s_ccdny_f= np.zeros(n_time_entries)

cti_beta_s_ccdpy_e = np.zeros(n_time_entries)
cti_beta_s_ccdpy_f = np.zeros(n_time_entries)
cti_beta_s_ccdny_e = np.zeros(n_time_entries)
cti_beta_s_ccdny_f = np.zeros(n_time_entries)


cti_alpha_p_ccdpy_e = [0.15, 0.15, 0.15, 0.15, 0.15]
cti_alpha_p_ccdpy_f = [0.15, 0.15, 0.15, 0.15, 0.15] 
cti_alpha_p_ccdny_e = [0.15, 0.15, 0.15, 0.15, 0.15]
cti_alpha_p_ccdny_f = [0.15, 0.15, 0.15, 0.15, 0.15]

cti_beta_p_ccdpy_e = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdpy_f = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_e = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_p_ccdny_f = [0.2, 0.2, 0.2, 0.2, 0.2]

cti_alpha_s_ccdpy_e = [0.15, 0.15, 0.15, 0.15, 0.15]
cti_alpha_s_ccdpy_f = [0.15, 0.15, 0.15, 0.15, 0.15]
cti_alpha_s_ccdny_e = [0.15, 0.15, 0.15, 0.15, 0.15]
cti_alpha_s_ccdny_f = [0.15, 0.15, 0.15, 0.15, 0.15]

cti_beta_s_ccdpy_e = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_s_ccdpy_f = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_s_ccdny_e = [0.2, 0.2, 0.2, 0.2, 0.2]
cti_beta_s_ccdny_f = [0.2, 0.2, 0.2, 0.2, 0.2]


ttime = pyfits.Column(name='TIME', unit='s', format='D', array=met)
ttrap_density_image_py = pyfits.Column(name='IMA_DENSITY_PY', format=form, array=trap_density_image_py)
ttrap_density_serial_py = pyfits.Column(name='SERIAL_DENSITY_PY', format=form, array=trap_density_serial_py)
ttrap_density_image_ny = pyfits.Column(name='IMA_DENSITY_NY', format=form, array=trap_density_image_ny)
ttrap_density_serial_ny = pyfits.Column(name='SERIAL_DENSITY_NY', format=form, array=trap_density_serial_ny)


tg_ccdpy_e = pyfits.Column(name='G_PY_E', format=form, array=gain_ccdpy_e)
tg_ccdpy_f = pyfits.Column(name='G_PY_F', format=form, array=gain_ccdpy_f)
tg_ccdny_e = pyfits.Column(name='G_NY_E', format=form, array=gain_ccdny_e)
tg_ccdny_f = pyfits.Column(name='G_NY_F', format=form, array=gain_ccdny_f)

tcti_alpha_p_ccdpy_e = pyfits.Column(name='CTI_ALPHA_P_PY_E', format=form, array=cti_alpha_p_ccdpy_e)
tcti_alpha_p_ccdpy_f = pyfits.Column(name='CTI_ALPHA_P_PY_F', format=form, array=cti_alpha_p_ccdpy_f) 
tcti_alpha_p_ccdny_e = pyfits.Column(name='CTI_ALPHA_P_NY_E', format=form, array=cti_alpha_p_ccdny_e) 
tcti_alpha_p_ccdny_f = pyfits.Column(name='CTI_ALPHA_P_NY_F', format=form, array=cti_alpha_p_ccdny_f) 

tcti_beta_p_ccdpy_e = pyfits.Column(name='CTI_BETA_P_PY_E', format=form, array=cti_beta_p_ccdpy_e)
tcti_beta_p_ccdpy_f = pyfits.Column(name='CTI_BETA_P_PY_F', format=form, array=cti_beta_p_ccdpy_f) 
tcti_beta_p_ccdny_e = pyfits.Column(name='CTI_BETA_P_NY_E', format=form, array=cti_beta_p_ccdny_e) 
tcti_beta_p_ccdny_f = pyfits.Column(name='CTI_BETA_P_NY_F', format=form, array=cti_beta_p_ccdny_f) 

tcti_alpha_s_ccdpy_e = pyfits.Column(name='CTI_ALPHA_S_PY_E', format=form, array=cti_alpha_s_ccdpy_e)
tcti_alpha_s_ccdpy_f = pyfits.Column(name='CTI_ALPHA_S_PY_F', format=form, array=cti_alpha_s_ccdpy_f) 
tcti_alpha_s_ccdny_e = pyfits.Column(name='CTI_ALPHA_S_NY_E', format=form, array=cti_alpha_s_ccdny_e) 
tcti_alpha_s_ccdny_f = pyfits.Column(name='CTI_ALPHA_S_NY_F', format=form, array=cti_alpha_s_ccdny_f) 

tcti_beta_s_ccdpy_e = pyfits.Column(name='CTI_BETA_S_PY_E', format=form, array=cti_beta_s_ccdpy_e)
tcti_beta_s_ccdpy_f = pyfits.Column(name='CTI_BETA_S_PY_F', format=form, array=cti_beta_s_ccdpy_f) 
tcti_beta_s_ccdny_e = pyfits.Column(name='CTI_BETA_S_NY_E', format=form, array=cti_beta_s_ccdny_e) 
tcti_beta_s_ccdny_f = pyfits.Column(name='CTI_BETA_S_NY_F', format=form, array=cti_beta_s_ccdny_f) 

# GAIN EXTENSION FOR CCDnY

# GAIN EXTENSION FOR CCDnp

coldefs_py = pyfits.ColDefs([ttime, 
                          ttrap_density_image_py,
                          ttrap_density_serial_py,
                          tg_ccdpy_e,
                          tg_ccdpy_f,
                          tcti_alpha_p_ccdpy_e,
                          tcti_alpha_p_ccdpy_f,
                          tcti_beta_p_ccdpy_e,
                          tcti_beta_p_ccdpy_f,
                          tcti_alpha_s_ccdpy_e,
                          tcti_alpha_s_ccdpy_f,
                          tcti_beta_s_ccdpy_e,
                          tcti_beta_s_ccdpy_f,
                          ])

# pyfits
# ghdu = pyfits.new_table(coldefs)
# astropy.io.fits (or pyfits)
gpyhdu = pyfits.BinTableHDU.from_columns(coldefs_py)

gpyhdu.header.set('EXTNAME', 'GAIN_CCDpY',
                'Name of the binary table extension')

gpyhdu.header.set('TELESCOP', 'SMILE', 'Telescope (mission) name')
gpyhdu.header.set('INSTRUME', 'SXI', 'Instrument name')
gpyhdu.header.set('ORIGIN', 'LU', 'Source of FITS file')

# creator = str(pyfits.__package__) + ' ' + str(pyfits.__version__)
# astropy
creator = str(pyfits.__package__) + ' ' + str(astropy.__version__)
gpyhdu.header.set('CREATOR', creator, 'Creator')
gpyhdu.header.set('VERSION', 'v1', 'Extension version number')
gpyhdu.header.set('FILENAME', 'sxi_gain_v5.fits', 'File Name')

# mode
content = 'Photon Counting FF Gain CCDpY'


gpyhdu.header.set('CONTENT', content, 'File content')
gpyhdu.header.set('TIMESYS', 'TT     ', 'Time system')
gpyhdu.header.set('MJDREFI', 51910, 'Reference MJD Integer part')
gpyhdu.header.set('MJDREFF', 7.4287037E-4, 'Reference MJD fractional')
gpyhdu.header.set('CLOCKAPP',
                False,
                'If clock corrections are applied (F/T)')

gpyhdu.header.set('CCLS0001', 'BCF',
                'Dataset is a Basic Calibration File')
gpyhdu.header.set('CDTP0001', 'DATA',
                'Calibration file contains data')
gpyhdu.header.set('CCNM0001', 'GAIN2', 'Type of Calibration data')

datamode = 'DATAMODE(' + 'PC' + ')'


gpyhdu.header.set('CDES0001',
                'SXI CCDpY Gain CTI coefficients',
                'Description')

gpyhdu.header.set('NOM_GAIN', 10.,
                'Nominal CCDpY Gain for PI channels')

gain_py_HDU_Comments = [
    "",
    "This extension provides CCDpY Gain and CTI coefficients",
    "for the PI (pulse height invariant) energy scale calculation."
    ]

for comment in gain_py_HDU_Comments:
    gpyhdu.header.add_comment(comment)

gpyhdu.header.set('DATE', creationdate,
                'File creation date (YYYY-MM-DDTHH:MM:SS UT)')

for item in gpyhdu.header.items():
    print(item)


    

coldefs_ny = pyfits.ColDefs([ttime, 
                          ttrap_density_image_ny,
                          ttrap_density_serial_ny,
                          tg_ccdny_e,
                          tg_ccdny_f,
                          tcti_alpha_p_ccdny_e,
                          tcti_alpha_p_ccdny_f,
                          tcti_beta_p_ccdny_e,
                          tcti_beta_p_ccdny_f,
                          tcti_alpha_s_ccdny_e,
                          tcti_alpha_s_ccdny_f,
                          tcti_beta_s_ccdny_e,
                          tcti_beta_s_ccdny_f,
                          ])

# pyfits
# ghdu = pyfits.new_table(coldefs)
# astropy.io.fits (or pyfits)
gnyhdu = pyfits.BinTableHDU.from_columns(coldefs_ny)

gnyhdu.header.set('EXTNAME', 'GAIN_CCDnY',
                'Name of the binary table extension')

gnyhdu.header.set('TELESCOP', 'SMILE', 'Telescope (mission) name')
gnyhdu.header.set('INSTRUME', 'SXI', 'Instrument name')
gnyhdu.header.set('ORIGIN', 'LU', 'Source of FITS file')

# creator = str(pyfits.__package__) + ' ' + str(pyfits.__version__)
# astropy
creator = str(pyfits.__package__) + ' ' + str(astropy.__version__)
gnyhdu.header.set('CREATOR', creator, 'Creator')
gnyhdu.header.set('VERSION', 'v5', 'Extension version number')
gnyhdu.header.set('FILENAME', 'sxi_gain_v5.fits', 'File Name')

# mode
content = 'Photon Counting FF Gain'


gnyhdu.header.set('CONTENT', content, 'File content')
gnyhdu.header.set('TIMESYS', 'TT     ', 'Time system')
gnyhdu.header.set('MJDREFI', 51910, 'Reference MJD Integer part')
gnyhdu.header.set('MJDREFF', 7.4287037E-4, 'Reference MJD fractional')
gnyhdu.header.set('CLOCKAPP',
                False,
                'If clock corrections are applied (F/T)')

gnyhdu.header.set('CCLS0001', 'BCF',
                'Dataset is a Basic Calibration File')
gnyhdu.header.set('CDTP0001', 'DATA',
                'Calibration file contains data')
gnyhdu.header.set('CCNM0001', 'GAIN2', 'Type of Calibration data')

datamode = 'DATAMODE(' + 'PC' + ')'


gnyhdu.header.set('CDES0001',
                'SXI CCDnY Gain CTI coefficients',
                'Description')

gnyhdu.header.set('NOM_GAIN', 10.,
                'Nominal CCDnY Gain for PI channels')

gain_ny_HDU_Comments = [
    "",
    "This extension provides CCDnY Gain and CTI coefficients",
    "for the PI (pulse height invariant) energy scale calculation."
]

for comment in gain_ny_HDU_Comments:
    gnyhdu.header.add_comment(comment)

gnyhdu.header.set('DATE', creationdate,
                'File creation date (YYYY-MM-DDTHH:MM:SS UT)')

for item in gnyhdu.header.items():
    print(item)

# Finally append the HDUs and create the FITs file.
fitsobj = pyfits.HDUList([phdu, gpyhdu, gnyhdu])


fitsobj.info()

fitsobj.writeto('sxi_gain_v5.fits')
