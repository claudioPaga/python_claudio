#!/usr/bin/env python

import numpy as np
import astropy
import astropy.io.fits as pyfits
import datetime
import argparse
import sys
# import pylab as plt

description = "Create a FITs format XRT gain file.\n"

parser = argparse.ArgumentParser(description=description)

parser.add_argument("--gcnfile",
                    help="Input GCn coeffs file (required)",
                    default='none')

parser.add_argument("--trapfileinputs",
                    help="Input file of epoch, alpha1/2 ebreak traptables" +
                    " (required)",
                    default='none')

parser.add_argument("--trapfilepath",
                    help="Path to trap tables" +
                    " default: %(default)s)",
                    default='./traptables/')

parser.add_argument("--outfile",
                    help="Output FITs file" +
                    " (default: %(default)s)",
                    default='gain.fits')

parser.add_argument("--mode",
                    help="Readout mode (PHOTON or WINDOWED;" +
                    " default: %(default)s)",
                    default='PHOTON')

parser.add_argument("--vsub",
                    help="Substrate voltage (0 or 6;" +
                    " default: %(default)s)",
                    default='6')

parser.add_argument("--fromdate",
                    help="Date file to be used from" +
                    " (default: %(default)s)",
                    default='2001-01-01')

parser.add_argument("--fromtime",
                    help="Time file to be used from" +
                    " (default: %(default)s)",
                    default='00:00:00')

parser.add_argument("--version",
                    help="gain file version" +
                    " (default: %(default)s)",
                    default='0')

# NB - the beta1,2, etc are now read in from the gcnfile

# parser.add_argument("--beta1",
#                     help="beta1" +
#                     " (default: %(default)s)",
#                     default='0.25')

# parser.add_argument("--beta2",
#                     help="beta2" +
#                     " (default: %(default)s)",
#                     default='0.2')

# parser.add_argument("--e_cti",
#                     help="PI_NOM scale CTI break energy (eV)" +
#                     " (default: %(default)s)",
#                     default='5895.0')

# parser.add_argument("--beta1_trap",
#                     help="beta1_trap" +
#                     " (default: %(default)s)",
#                     default='0.25')

# parser.add_argument("--beta2_trap",
#                     help="beta2_trap" +
#                     " (default: %(default)s)",
#                     default='0.2')

# parser.add_argument("--e_cti_trap",
#                     help="PI scale CTI break energy (eV)" +
#                     " (default: %(default)s)",
#                     default='5895.0')

parser._optionals.title = "Arguments"

args = parser.parse_args()

print(args)

gcnfile = args.gcnfile

if gcnfile == 'none':
    sys.stderr.write('File of GCn coeffs required\n')
    sys.exit(1)

trapfileinputs = args.trapfileinputs

if trapfileinputs == 'none':
    sys.stderr.write('File of trap-table inputs required\n')
    sys.exit(1)

trapfilepath = args.trapfilepath

outfile = args.outfile

# readout mode
mode = args.mode.upper()
if not((mode == 'PHOTON') or (mode == 'WINDOWED')):
    sys.stderr.write('mode ' + mode + ' not recognised\n')
    sys.exit(1)

ngrade = 0
if mode == 'PHOTON':
    ngrade = 33
elif mode == 'WINDOWED':
    ngrade = 16
else:
    print('mode ' + mode + ' not recognised!')
    sys.exit(1)

# substrate voltage
vsub = str(args.vsub)
if not((vsub == '0') or (vsub == '6')):
    sys.stderr.write('vsub ' + vsub + ' not recognised\n')
    sys.exit(1)

# UTC date when file should first be used
fromdate = args.fromdate
# UTC time when file should first be used
fromtime = args.fromtime

# gain file version number
version = str(args.version)

# file creation date (format is YYYY-MM-DDTHH:MM:SS)
creationdate = datetime.datetime.utcnow().replace(microsecond=0).isoformat()

# Load new  coeffs. Should contain :
# MET, CCDTEMP_j, GCn_j, BETA1, BETA2, E_CTI, \
#                 GCn_TRAP_j, BETA1_TRAP, BETA2_TRAP, E_CTI_TRAP
# where  j = 0, 1, 2 (i.e. 3 CCDTEMPs ordered from high to low)
#        n = (0, 1, 2, 3) (i.e. coeffs
#                          GC0 GC1 GC2 GC3
#                          GC0_TRAP GC1_TRAP GC2_TRAP, GC3_TRAP
# each at 3 temperatures).
# NB this expects ccdtemp is stored in order high to low like the FITS file.
zz = np.loadtxt(gcnfile)

if zz.shape[1] != 34:
    msg = 'Incorrect number of columns in ' + gcnfile
    sys.stderr.write(msg + '\n')
    sys.exit(1)
    
i = 0
met = zz[:, i]
i = i + 1
ccdtemp = zz[:, i:i+3]
i = i + 3
gc0 = zz[:, i:i+3]
i = i + 3
gc1 = zz[:, i:i+3]
i = i + 3
gc2 = zz[:, i:i+3]
i = i + 3
gc3 = zz[:, i:i+3]
i = i + 3
beta1 = zz[:, i]
i = i + 1
beta2 = zz[:, i]
i = i + 1
# e_cti is in PI
e_cti = zz[:, i]
i = i + 1
gc0_trap = zz[:, i:i+3]
i = i + 3
gc1_trap = zz[:, i:i+3]
i = i + 3
gc2_trap = zz[:, i:i+3]
i = i + 3
gc3_trap = zz[:, i:i+3]
i = i + 3
beta1_trap = zz[:, i]
i = i + 1
beta2_trap = zz[:, i]
i = i + 1
# e_cti_trap is in PI
e_cti_trap = zz[:, i]
i = i + 1

# Don't forget GC4 and GC5
gc4 = np.zeros_like(gc0)
gc5 = np.zeros_like(gc0)
#
gc4_trap = np.zeros_like(gc0)
gc5_trap = np.zeros_like(gc0)

# Number of time intervals and CCDTemps
ntime, nccdtemp = ccdtemp.shape

# # Energy dependent CTI.
# # beta1, beta2, e_cti has never been changed
# # nominal
# beta1 = np.zeros((ntime), dtype=np.float)
# beta2 = np.zeros((ntime), dtype=np.float)
# e_cti = np.zeros((ntime), dtype=np.float)

# # e_cti is in PI
# e_cti[:] = 5895./10.
# beta1[:] = float(args.beta1)
# beta2[:] = float(args.beta2)
# # e_cti is in PI
# e_cti[:] = float(args.e_cti) / 10.

# # trapped
# beta1_trap = np.zeros((ntime), dtype=np.float)
# beta2_trap = np.zeros((ntime), dtype=np.float)
# e_cti_trap = np.zeros((ntime), dtype=np.float)

# beta1_trap[:] = float(args.beta1_trap)
# beta2_trap[:] = float(args.beta2_trap)
# # e_cti is in PI
# e_cti_trap[:] = float(args.e_cti_trap) / 10.

# Now for the trap tables.
# Read the file pointing to the actual epoch dependent trap table files
trapfiles_epoch = []
traps_alpha1 = []
traps_alpha2 = []
traps_ebreak = []
trapfiles = []
offset_scale_files = []

f = open(trapfileinputs)

for line in f.readlines():
    if line[0] == '#':
        continue

    if line == '':
        continue
    
    words = (line.strip()).split()

    print(words)

    trapfiles_epoch.append(float(words[0]))
    traps_alpha1.append(float(words[1]))
    traps_alpha2.append(float(words[2]))
    # convert ebreak from eV to PI
    traps_ebreak.append(float(words[3]) / 10.)
    trapfiles.append(words[4])

    # Is the offset_scale_file given ?
    # If not, set to 'none' for this trap epoch
    if len(words) == 6:
        offset_scale_files.append(words[5])
    else:
        offset_scale_files.append('none')
        
f.close()

# read the trap files - expected format DETX, DETY, YEXTENT, OFFSET (in eV)
# convert to RAWX, RAWY

# container for the trap data for different epochs
trapsdata = []

# maximum number of traps
maxntrap = 0

for trapfile in trapfiles:

    traps = np.loadtxt(trapfilepath + trapfile)

    traps_rawx = np.array(traps[:, 0], dtype=np.int) - 1
    traps_rawy = np.array(traps[:, 1], dtype=np.int) - 1
    traps_yextent = np.array(traps[:, 2], dtype=np.int)

    # convert from eV to PI
    # _0, 1, 2 are -48, -61.5, -75
    traps_offset_0 = traps[:, 3] / 10.
    traps_offset_1 = traps[:, 4] / 10.
    traps_offset_2 = traps[:, 5] / 10.

    trapsdata.append((traps_rawx, traps_rawy, traps_yextent,
                     traps_offset_0, traps_offset_1, traps_offset_2))

    # max number of traps
    ntrap = traps_offset_0.size

    if ntrap > maxntrap:
        maxntrap = ntrap

print("Maximum number of traps, maxntrap : ", maxntrap)

# Read in the per trap epoch grade dependent offset_scales, if they are
# provided.
# Otherwise offset_scales are set to 1.0 for that epoch.
offset_scale_data = []

for offset_scale_file in offset_scale_files:
    offset_scale_out = np.ones(ngrade, dtype=np.float)

    if offset_scale_file != 'none':
        
        # offset_scale_in contains array of (grade, scale)
        offset_scale_in = np.loadtxt(trapfilepath + offset_scale_file)

        # map input to full ngrade output array
        for ii in range(len(offset_scale_in)):
            tgrade, scale = offset_scale_in[ii]
            grade = int(tgrade)
            offset_scale_out[grade] = scale
        
    offset_scale_data.append(offset_scale_out)
    
# initialise arrays
rawx = np.zeros((ntime, maxntrap), dtype=np.int)
rawy = np.zeros((ntime, maxntrap), dtype=np.int)
yextent = np.zeros((ntime, maxntrap), dtype=np.int)
alpha1 = np.zeros((ntime), dtype=np.float)
alpha2 = np.zeros((ntime), dtype=np.float)
ebreak = np.ones((ntime), dtype=np.float)
offset = np.zeros((ntime, maxntrap*nccdtemp), dtype=np.float)
offset_scale = np.ones((ntime, ngrade), dtype=np.float)

# default values (no traps)
rawx[:, :] = -1
rawy[:, :] = -1
yextent[:, :] = -1

# Loop through all rows
for i in range(ntime):

    # Loop through trap epochs
    for ti in range(len(trapfiles)):

        trapsepoch = trapfiles_epoch[ti]

        if met[i] >= trapsepoch:

            (traps_rawx,
             traps_rawy,
             traps_yextent,
             traps_offset_0,
             traps_offset_1,
             traps_offset_2) = trapsdata[ti]

            alpha1[i] = traps_alpha1[ti]
            alpha2[i] = traps_alpha2[ti]
            ebreak[i] = traps_ebreak[ti]

            # reinitialise this row as it could have been set  in the
            # last iteration
            rawx[i, :] = -1
            rawy[i, :] = -1
            yextent[i, :] = -1
            offset[i, :] = -1

            rawx[i, 0:traps_rawx.size] = traps_rawx
            rawy[i, 0:traps_rawy.size] = traps_rawy
            yextent[i, 0:traps_yextent.size] = traps_yextent

            for j in range(traps_offset_0.size):
                offset[i, (j * nccdtemp + 0)] = traps_offset_0[j]
                offset[i, (j * nccdtemp + 1)] = traps_offset_1[j]
                offset[i, (j * nccdtemp + 2)] = traps_offset_2[j]
                
            # update the offset_scale values if we have them
            if len(offset_scale_data[ti]) > 0:
                offset_scale[i] = offset_scale_data[ti]

#
# create a primary hdu
#
phdu = pyfits.PrimaryHDU()
phdu.header.set('TELESCOP', 'SWIFT', 'Telescope (mission) name')
phdu.header.set('INSTRUME', 'XRT', 'Instrument name')

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
form = str(nccdtemp) + 'E'
form2 = str(maxntrap) + 'I'
form3 = str(maxntrap * nccdtemp) + 'E'
form4 = str(ngrade) + 'E'

ttime = pyfits.Column(name='TIME', unit='s', format='D', array=met)

tccdtemp = pyfits.Column(name='CCDTEMP', unit='C', format=form, array=ccdtemp)

tgc0 = pyfits.Column(name='GC0', format=form, array=gc0)
tgc1 = pyfits.Column(name='GC1', format=form, array=gc1)
tgc2 = pyfits.Column(name='GC2', format=form, array=gc2)
tgc3 = pyfits.Column(name='GC3', format=form, array=gc3)
tgc4 = pyfits.Column(name='GC4', format=form, array=gc4)
tgc5 = pyfits.Column(name='GC5', format=form, array=gc5)

tgc0_trap = pyfits.Column(name='GC0_TRAP', format=form, array=gc0_trap)
tgc1_trap = pyfits.Column(name='GC1_TRAP', format=form, array=gc1_trap)
tgc2_trap = pyfits.Column(name='GC2_TRAP', format=form, array=gc2_trap)
tgc3_trap = pyfits.Column(name='GC3_TRAP', format=form, array=gc3_trap)
tgc4_trap = pyfits.Column(name='GC4_TRAP', format=form, array=gc4_trap)
tgc5_trap = pyfits.Column(name='GC5_TRAP', format=form, array=gc5_trap)

tbeta1 = pyfits.Column(name='BETA1', format='1E', array=beta1)
tbeta2 = pyfits.Column(name='BETA2', format='1E', array=beta2)
te_cti = pyfits.Column(name='E_CTI', format='1E', array=e_cti)

tbeta1_trap = pyfits.Column(name='BETA1_TRAP', format='1E', array=beta1_trap)
tbeta2_trap = pyfits.Column(name='BETA2_TRAP', format='1E', array=beta2_trap)
te_cti_trap = pyfits.Column(name='E_CTI_TRAP', format='1E', array=e_cti_trap)

trawx = pyfits.Column(name='RAWX', format=form2, array=rawx)
trawy = pyfits.Column(name='RAWY', format=form2, array=rawy)
tyextent = pyfits.Column(name='YEXTENT', format=form2, array=yextent)

talpha1 = pyfits.Column(name='ALPHA1', format='1E', array=alpha1)
talpha2 = pyfits.Column(name='ALPHA2', format='1E', array=alpha2)
tebreak = pyfits.Column(name='EBREAK', format='1E', array=ebreak)

toffset = pyfits.Column(name='OFFSET', format=form3, array=offset)

toffset_scale = pyfits.Column(name='OFFSET_SCALE', format=form4,
                              array=offset_scale)

coldefs = pyfits.ColDefs([ttime,
                          tccdtemp,
                          tgc0,
                          tgc1,
                          tgc2,
                          tgc3,
                          tgc4,
                          tgc5,
                          tbeta1,
                          tbeta2,
                          te_cti,
                          tgc0_trap,
                          tgc1_trap,
                          tgc2_trap,
                          tgc3_trap,
                          tgc4_trap,
                          tgc5_trap,
                          tbeta1_trap,
                          tbeta2_trap,
                          te_cti_trap,
                          trawx,
                          trawy,
                          tyextent,
                          talpha1,
                          talpha2,
                          tebreak,
                          toffset,
                          toffset_scale
                          ])

# pyfits
# ghdu = pyfits.new_table(coldefs)
# astropy.io.fits (or pyfits)
ghdu = pyfits.BinTableHDU.from_columns(coldefs)

ghdu.header.set('EXTNAME', 'GAIN',
                'Name of the binary table extension')

ghdu.header.set('TELESCOP', 'SWIFT', 'Telescope (mission) name')
ghdu.header.set('INSTRUME', 'XRT', 'Instrument name')
ghdu.header.set('ORIGIN', 'LU', 'Source of FITS file')
# pyfits
# creator = str(pyfits.__package__) + ' ' + str(pyfits.__version__)
# astropy
creator = str(pyfits.__package__) + ' ' + str(astropy.__version__)
ghdu.header.set('CREATOR', creator, 'Creator')
ghdu.header.set('VERSION', version, 'Extension version number')
ghdu.header.set('FILENAME', outfile, 'File Name')

# mode
content = 'Photon Counting Gain'
if mode == 'WINDOWED':
    content = 'Windowed Timing Gain'

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

datamode = 'DATAMODE(' + mode + ')'

ghdu.header.set('CBD10001', datamode, 'Parameter Boundary')
xrtvsub = 'XRTVSUB(' + str(vsub) + ')'
ghdu.header.set('CBD20001', xrtvsub, 'Parameter Boundary')
ghdu.header.set('CVSD0001', fromdate,
                'UTC date when file should first be used')
ghdu.header.set('CVST0001', fromtime,
                'UTC time when file should first be used')

ghdu.header.set('CDES0001',
                'XRT Gain with Position & Time Dependent CTI coefficients',
                'Description')

ghdu.header.set('NOM_GAIN', 10.,
                'Nominal Gain for PI channels')

gain_HDU_Comments = [
    "",
    "This extension provides CTI and charge trap coeff.s for the PI (pulse",
    "height invariant) energy scale calculation as follows:",
    "",
    "1) PI calculation without charge trap correction:",
    "",
    "  PI(PHA) = (int)((PHA*(C0 + x*C1 + y*C2) + C3 + x*C4 + y*C5) / G)",
    "",
    "2) PI calculation with charge trap correction:",
    "",
    "  PI(PHA) = (int)((PHA*(GC0_TRAP+x*GC1_COR_TRAP+y*GC2_COR_TRAP)+",
    "                   GC3_TRAP+x*GC4_TRAP+y*GC5_TRAP)/G)+OFFSET_COR",
    "",
    "where",
    "",
    "  x == RAWX",
    "  y == RAWY",
    "  G == NOM_GAIN",
    "",
    "  GC1_COR = GC1*(E/E_CTI)**(-BETA1) E<=E_CTI",
    "  GC1_COR = GC1*(E/E_CTI)**(-BETA2) E> E_CTI",
    "  GC2_COR = GC2*(E/E_CTI)**(-BETA1) E<=E_CTI",
    "  GC2_COR = GC2*(E/E_CTI)**(-BETA2) E> E_CTI",
    "",
    "  GC1_COR_TRAP = GC1_TRAP*(E/E_CTI_TRAP)**(-BETA1_TRAP) for E<=E_CTI",
    "  GC1_COR_TRAP = GC1_TRAP*(E/E_CTI_TRAP)**(-BETA2_TRAP) for E> E_CTI",
    "  GC2_COR_TRAP = GC2_TRAP*(E/E_CTI_TRAP)**(-BETA1_TRAP) for E<=E_CTI",
    "  GC2_COR_TRAP = GC2_TRAP*(E/E_CTI_TRAP)**(-BETA2_TRAP) for E> E_CTI",
    "",
    "  OFFSET_COR = (OFFSET*OFFSET_SCALE)*(E/EBREAK)**(ALPHA1) for E<=EBREAK",
    "  OFFSET_COR = (OFFSET*OFFSET_SCALE)*(E/EBREAK)**(ALPHA2) for E> EBREAK",
    "",
    "with", 
    "  OFFSET_SCALE = grade dependent trap OFFSET scaling factor",
    ""
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

fitsobj.writeto(outfile)
