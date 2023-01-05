import spline2dbase
import h5py as h5py
import numpy as np
from fitting_psf import from_mm_2_pix, from_pix_2_mm
from imagette import list_psf, barycenter, gauss
import math
import scipy.signal

PSFDIR = './'
# PSFDIR = '/home/reza/plato/share/psf/Sep17_real_MC_T1413/'

# Let's define the PSF parameters
sizex = 8  # physical size of the PSF (x-direction)
sizey = 8  # physical size of the PSF (y-direction)
subres = 128  # resolution of the PSF
bsres = 20  # resolution of the b-spline decomposition of the PSF

# The second thing to do is to open the .hdf5 file that containg all the PSFs from biruni3
file_hdf5 = h5py.File(PSFDIR+'PSF.hdf5', 'r')

# The third thing to do is to define the parameters for the Diffusion Kernel to covolve the PSFs
DifKerSize = 3  # Size [pixel]
DifKerWidth = 0.2  # width [pixel] 0.1 -> 99.99% in the central pixel ; 0.2 -> 97.5% ; 0.3 -> 81.8%% ; 0.5 -> 46.8%

bsres = 20 # resolution adopted for the b-spline decomposition

# Now we build the diffusion kernel, a Gaussian function of size DifKerSize x DifKerSize centered on the middle of the central pixel
GaussKernel = gauss(math.floor(DifKerSize / 2.) + 0.5, math.floor(DifKerSize / 2.) + 0.5, DifKerWidth, DifKerSize,
                    subres=subres)
GaussKernel /= GaussKernel.sum()


# Second, we load the PSF coordinates from the list file
xpsf, ypsf = list_psf(PSFDIR+'list')

# Third, we convert the PSF coordinates on the Focal Planet from mm to pixel
xpsf_pix, ypsf_pix = from_mm_2_pix(xpsf, ypsf)

# We compute the folowing parameters
lx = bsres * sizex
ly = bsres * sizey
npsf = len(xpsf)
psfbs = np.zeros((npsf, lx, ly))
pxc = np.zeros(npsf)
pyc = np.zeros(npsf)
print('Processing the PSF')
for k in range(npsf):
    psf_id = str(k + 1)

    psf = np.array(file_hdf5[psf_id])
    # Convolving now the optical PSF by the Gaussian kernel defined previously
    psf = scipy.signal.fftconvolve(psf, GaussKernel, mode='same')
    # Now we normalize the psf
    psf /= psf.sum()
    # We compute the barycenter of the psf
    pxc[k], pyc[k] = barycenter(psf, subres=subres)

    # Then we convert the PSFs to b-spline
    psfbs[k] = spline2dbase.Pixel2Spline(psf, lx, ly)

np.savez('PSF.npz', psfbs=psfbs, pxc=pxc, pyc=pyc, xpsf_pix=xpsf_pix, ypsf_pix=ypsf_pix)


file_hdf5.close()

