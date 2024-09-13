import numpy as np  # type: ignore
import scipy.signal  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math
import h5py  # type: ignore
from imagette import gauss, barycenter
import spline2dbase  # type: ignore

PSFPATH = '/home/fercho/double-aperture-photometry/PSF_Focus_0mu.hdf5'
file_hdf5 = h5py.File(PSFPATH, 'r')

# Let's define the PSF parameters
sizex = 8  # physical size of the PSF (x-direction)
sizey = 8  # physical size of the PSF (y-direction)
subres = 64  # resolution of the PSF
bsres = 20  # resolution of the b-spline decomposition of the PSF
cmap = 'seismic'  # matplotlib color map for the plots
fsize= 15
# PSF indices you want to visualize
psf_indices = [2, 200, 4200]

# Define the parameters for the Diffusion Kernel to convolve the PSFs
DifKerSize = 3  # Size [pixel]
DifKerWidth = 0.2  # width [pixel]

# Build the diffusion kernel
GaussKernel = gauss(math.floor(DifKerSize / 2.) + 0.5, math.floor(DifKerSize / 2.) + 0.5, DifKerWidth, DifKerSize, subres=subres)
GaussKernel /= GaussKernel.sum()

# Create a 3x3 subplot (3 PSFs in 3 stages each)
fig, axs = plt.subplots(3, 3, sharex=False, sharey=False, figsize=(24, 14))

for i, k in enumerate(psf_indices):
    psf_id = str(k + 1)  # hdf5 file uses 1-based index
    print('PSF ID = %s ' % psf_id)

    psf_raw = np.array(file_hdf5[psf_id])

    # Convolve the optical PSF by the Gaussian kernel
    psf = scipy.signal.fftconvolve(psf_raw, GaussKernel, mode='same')
    psf /= psf.sum()  # Normalize the PSF

    # Compute the barycenter
    pxc, pyc = barycenter(psf, subres=subres)

    # Compute alpha angle
    xpsf_mm = float(file_hdf5[psf_id].attrs['centerCoordinates1'])
    ypsf_mm = float(file_hdf5[psf_id].attrs['centerCoordinates2'])
    Rfp = math.sqrt(xpsf_mm**2 + ypsf_mm**2)
    alpha = np.arctan(np.sqrt(xpsf_mm ** 2 + ypsf_mm ** 2) / 247.732)

    print('Xfp, Yfp, Rfp  [mm] = %g, %g, %g' % (xpsf_mm, ypsf_mm, Rfp))
    print('alpha [deg] = %g ' % (alpha * 180.0 / np.pi))

    # Compute the imagette
    lx = bsres * sizex
    ly = bsres * sizey
    psfbs = spline2dbase.Pixel2Spline(psf, lx, ly)
    imagette = spline2dbase.Spline2Imagette(np.ascontiguousarray(psfbs), bsres, 8, 8)

    # Plotting each PSF at different stages in the respective column
    axs[0, i].imshow(psf_raw[subres:-subres, subres:-subres], origin='lower', cmap=cmap)
    axs[0, i].set_title(r'$\mathbf{\alpha = %0.2f^\circ}$ (Raw)' % (alpha * 180.0 / np.pi), fontsize=fsize, fontweight='bold')

    axs[1, i].imshow(psf[subres:-subres, subres:-subres], origin='lower', cmap=cmap)
    axs[1, i].set_title(r'$\mathbf{\alpha = %0.2f^\circ}$ (Convolved)' % (alpha * 180.0 / np.pi), fontsize=fsize, fontweight='bold')

    axs[2, i].imshow(imagette[1:-1, 1:-1], origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
    axs[2, i].set_title(r'$\mathbf{\alpha = %0.2f^\circ}$ (Window)' % (alpha * 180.0 / np.pi), fontsize=fsize, fontweight='bold')

# Adjust layout and show the plots
plt.subplots_adjust(wspace=0.2, hspace=0.4)
plt.tight_layout()  # Automatically adjust to minimize white space
plt.savefig('PSFs_fig_2.pdf', dpi=300)
plt.show()