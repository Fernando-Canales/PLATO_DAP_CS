import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import h5py as h5py # type: ignore
import scipy.signal # type: ignore
import math # type: ignore
import struct
from astropy.io import fits # type: ignore
from astropy.io import fits as pyfits # type: ignore
import spline2dbase # type: ignore

#Directory with all the PSFs
dataDIR= '/home/fercho/double-aperture-photometry/psf_fits/'
dataDIR=''
fsize = 16
cmap='viridis'

# Let's define the PSF parameters
sizex = 8  # physical size of the PSF (x-direction)
sizey = 8  # physical size of the PSF (y-direction)
subres = 128  # resolution of the PSF
bsres = 20  # resolution of the b-spline decomposition of the PSF


# Imagette characteristics
sizex2 = 6  # physical size of the imagette
sizey2 = 6
xc,yc  = 3.,3. # star centroid within the imagette

# Now we load the list file
radius = np.loadtxt(dataDIR + 'list', unpack=True, usecols=[1])
#fits_filename = np.loadtxt(dataDIR + 'list', unpack=True, usecols=[0])

#file_hdf5 = h5py.File(dataDIR+'PSF.hdf5', 'r')

#for k in range(0, 439):
#    psf_id = str(k + 1)

#    psf = np.array(file_hdf5[psf_id])
#    plt.imshow(psf, cmap='gray')
#    plt.colorbar()
#    plt.show()
 
#fits_names = np.genfromtxt(dataDIR + 'list', dtype=str, usecols=(0))

#for i in fits_names:
#    fits_file = fits.open(dataDIR + i)
#    fits_info = fits_file[0].data
    
#    fits_file.close()
#    plt.imshow(fits_info, origin='lower', extent=(0, 1024,0,1024), cmap='plasma')
#    plt.title(i)
#    plt.colorbar()
#    plt.show()

fits_file_0 = fits.open(dataDIR + '6000-01174-045000.fits')

fits_file_10 = fits.open(dataDIR + '6000-10452-045000.fits')

fits_file_18 =  fits.open(dataDIR + '6000-18837-347329.fits') 

fits_info_0 = fits_file_0[0].data

fits_info_10 = fits_file_10[0].data

fits_info_18 = fits_file_18[0].data

   


# Defining the function that produces a Gaussian kernel
def gauss(xc, yc, width, size, subres=1):
    # 2D Gaussian function centered on (xc,yc)
    s = float(subres)
    p = subres * size
    (RX, RY) = np.meshgrid(np.arange(0, p) / s - xc + 0.5 / s,
                           np.arange(0, p) / s - yc + 0.5 / s)
    D2 = RX * RX + RY * RY

    return np.exp(-D2 / (2 * width * width))

# The third thing to do is to define the parameters for the Diffusion Kernel to covolve the PSFs
DifKerSize = 3  # Size [pixel]
DifKerWidth = 0.2  # width [pixel] 0.1 -> 99.99% in the central pixel ; 0.2 -> 97.5% ; 0.3 -> 81.8%% ; 0.5 -> 46.8%

# Now we build the diffusion kernel, a Gaussian function of size DifKerSize x DifKerSize centered in the middle of the central pixel
GaussKernel = gauss(math.floor(DifKerSize / 2.) + 0.5, math.floor(DifKerSize / 2.) + 0.5, DifKerWidth, DifKerSize,
                    subres=subres)
GaussKernel /= GaussKernel.sum()



# Defining the barycenter
def barycenter(array,mask=None,x=None,y=None,subres=1):
    if(isinstance(x, (np.ndarray, np.generic) )== False):
        x=(np.arange(0,array.shape[1])+0.5)/float(subres)
    if(isinstance(y, (np.ndarray, np.generic) )== False):
        y=(np.arange(0,array.shape[0])+0.5)/float(subres)

    if(x.ndim ==1 & y.ndim ==1):
        x,y = np.meshgrid(x,y)

    if(mask is not None):
        weight=np.sum(array*mask)
        bx= np.sum(array*x*mask)/weight
        by= np.sum(array*y*mask)/weight
    else:
        tmp=np.sum(array)
        bx= np.sum(array*x)/tmp
        by= np.sum(array*y)/tmp
    return bx,by

psf_0_cv = scipy.signal.fftconvolve(fits_info_0, GaussKernel, mode='same')
psf_0 = np.array(pyfits.open(dataDIR + '6000-01174-045000.fits')[0].data,dtype=np.float)
# Now we normalize the psf
psf_0 /= psf_0.sum()
pxc_0,pyc_0 = barycenter(psf_0,subres=subres)
lx_0 = bsres*sizex
ly_0 = bsres*sizey
psfbs_0 = spline2dbase.Pixel2Spline(psf_0, lx_0 ,ly_0)

offx_0 = xc-pxc_0
offy_0 = yc-pyc_0

# imagette obtained by integrating the b-spline decomposition of the PSF
imgagette_0 = spline2dbase.Spline2Imagette(psfbs_0,bsres,sizex2,sizey2,offx=offx_0,offy=offy_0)

psf_10_cv = scipy.signal.fftconvolve(fits_info_10, GaussKernel, mode='same')
psf_10 = np.array(pyfits.open(dataDIR + '6000-10452-045000.fits')[0].data,dtype=np.float)
# Now we normalize the psf
psf_10 /= psf_10.sum()
pxc_10,pyc_10 = barycenter(psf_10,subres=subres)
lx_10 = bsres*sizex
ly_10 = bsres*sizey
psfbs_10 = spline2dbase.Pixel2Spline(psf_10, lx_10 ,ly_10)

offx_10 = xc-pxc_10
offy_10 = yc-pyc_10

# imagette obtained by integrating the b-spline decomposition of the PSF
imgagette_10 = spline2dbase.Spline2Imagette(psfbs_10,bsres,sizex2,sizey2,offx=offx_10,offy=offy_10)


psf_18_cv = scipy.signal.fftconvolve(fits_info_18, GaussKernel, mode='same')
psf_18 = np.array(pyfits.open(dataDIR + '6000-18837-347329.fits')[0].data,dtype=np.float)
    # Now we normalize the psf
psf_18 /= psf_18.sum()
pxc_18,pyc_18 = barycenter(psf_18,subres=subres)
lx_18 = bsres*sizex
ly_18 = bsres*sizey
psfbs_18 = spline2dbase.Pixel2Spline(psf_18, lx_18 ,ly_18)

offx_18 = xc-pxc_18
offy_18 = yc-pyc_18

# imagette obtained by integrating the b-spline decomposition of the PSF
imgagette_18 = spline2dbase.Spline2Imagette(psfbs_18,bsres,sizex2,sizey2,offx=offx_18,offy=offy_18)

fits_file_0.close()
fits_file_10.close()
fits_file_18.close()


fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(14, 14))

# Plot images with grid
axs[0, 0].imshow(fits_info_0, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[0, 0].tick_params(axis='y', which='both', length=5)
axs[0, 0].set_title(r' $\alpha = 0.1^{\circ}$', fontsize=fsize+3)
axs[0, 0].set_xticklabels(axs[0, 0].get_xticks(), fontsize=fsize)
axs[0, 0].set_yticklabels(axs[0, 0].get_yticks(), fontsize=fsize)

axs[0, 1].imshow(fits_info_10, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[0, 1].set_title(r' $\alpha = 10.04^{\circ}$', fontsize=fsize+3)
axs[0, 1].set_xticklabels(axs[0, 1].get_xticks(), fontsize=fsize)
axs[0, 1].set_yticklabels(axs[0, 1].get_yticks(), fontsize=fsize)

axs[0, 2].imshow(fits_info_18, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[0, 2].set_title(r' $\alpha = 18.84^{\circ}$', fontsize=fsize+3)
axs[0, 2].set_xticklabels(axs[0, 2].get_xticks(), fontsize=fsize)
axs[0, 2].set_yticklabels(axs[0, 2].get_yticks(), fontsize=fsize)

axs[1, 0].imshow(psf_0_cv, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[1, 0].set_xticks(np.arange(0, 7, 1))
axs[1, 0].set_yticks(np.arange(0, 7, 1))
axs[1, 0].tick_params(axis='x', which='both', length=0)
axs[1, 0].set_xticklabels(axs[1, 0].get_xticks(), fontsize=fsize)
axs[1, 0].set_yticklabels(axs[1, 0].get_yticks(), fontsize=fsize)

axs[1, 1].imshow(psf_10_cv, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[1, 1].tick_params(axis='y', which='both', length=0)
axs[1, 1].set_xticklabels(axs[1, 1].get_xticks(), fontsize=fsize)
axs[1, 1].set_yticklabels(axs[1, 1].get_yticks(), fontsize=fsize)

axs[1, 2].imshow(psf_18_cv, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[1, 2].tick_params(axis='y', which='both', length=0)
axs[1, 2].set_xticklabels(axs[1, 2].get_xticks(), fontsize=fsize)
axs[1, 2].set_yticklabels(axs[1, 2].get_yticks(), fontsize=fsize)

axs[2, 0].imshow(imgagette_0, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[2, 0].set_xticks(np.arange(0, 7, 1))
axs[2, 0].set_yticks(np.arange(0, 7, 1))
axs[2, 0].tick_params(axis='both', which='both', length=5)
axs[2, 0].set_xticklabels(axs[2, 0].get_xticks(), fontsize=fsize)
axs[2, 0].set_yticklabels(axs[2, 0].get_yticks(), fontsize=fsize)

axs[2, 1].imshow(imgagette_0, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[2, 1].set_xticks(np.arange(0, 7, 1))
axs[2, 1].set_yticks(np.arange(0, 7, 1))
axs[2, 1].tick_params(axis='both', which='both', length=5)
axs[2, 1].set_xticklabels(axs[2, 1].get_xticks(), fontsize=fsize)
axs[2, 1].set_yticklabels(axs[2, 1].get_yticks(), fontsize=fsize)

axs[2, 2].imshow(imgagette_18, origin='lower', extent=(0, 6, 0, 6), cmap=cmap)
axs[2, 2].set_xticks(np.arange(0, 7, 1))
axs[2, 2].set_yticks(np.arange(0, 7, 1))
axs[2, 2].tick_params(axis='both', which='both', length=5)
axs[2, 2].set_xticklabels(axs[2, 2].get_xticks(), fontsize=fsize)
axs[2, 2].set_yticklabels(axs[2, 2].get_yticks(), fontsize=fsize)

# Add colorbar
# cbar_ax = fig.add_axes([0.948, 0.15, 0.02, 0.7])
# fig.colorbar(axs[1, 2].images[0], cax=cbar_ax, orientation='vertical')

plt.subplots_adjust(wspace=0.1, hspace=0.2)
plt.savefig('PSFs_corrected.pdf')
plt.show()