import numpy as np


# Parameters relative to all the relevant paths
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues 
PSFfile = 'PSF.npz'                                                   # processed PSF files
DIRout = 'test_results/'                                              # storage directory


# Parameters for the imagette and PSF decomposition
size_im_x = 6  # size of the imagette (x-direction)
size_im_y = 6  # size of the imagette (y-direction)
subres = 128   # resolution of the PSF
bsres = 20     # resolution of the b-spline decomposition of the PSF

# Parameters for the NSR
sb = (45. * 21)  # Background noise from zodiacal light in e-/px(poisson noise)times integration time (21 sec.)
sd = 50.2        # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2         # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
#dback = 85000  # transit depth in ppm
#td = 4         # transit duration in hours
ntr = 3        # number of transits in one hour

# Parameters for the magnitude intervals
n_tar = 200                         # number of targets per magnitude interval
Pmin = 10                                # minimum magnitude
Pmax = 13                               # maximum magnitude
binsize = 0.5                           # binsize around every magnitude value
