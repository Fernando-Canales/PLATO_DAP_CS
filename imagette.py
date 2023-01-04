# This module will generate the imagette of a target and a contaminant for a given PSF
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import spline2dbase
import math
import scipy.special

def psf_gauss_int(xc,yc,width_x,width_y,sizex,sizey):
    '''
    Gaussian PSF integrated over the pixels of an imagette
    '''

    zx = 1./(math.sqrt(2.)*width_x)
    zy = 1./(math.sqrt(2.)*width_y)
    # Error function
    erfx = scipy.special.erf((np.arange(0,sizex+1) - xc)*zx)
    erfy = scipy.special.erf((np.arange(0,sizey+1) - yc)*zy)
    # normalization required to have the flux density equals to one  at the center
    ## u = math.sqrt(math.pi/2.)*width
    u = 0.5 # normalization required to have the mask integral equals to one
    maskx = (erfx[1:sizex+1] - erfx[0:sizex])*u
    masky = (erfy[1:sizey+1] - erfy[0:sizey])*u
    maskx = np.reshape(maskx,(1,sizex))
    masky = np.reshape(masky,(sizey,1))
    maskx = np.repeat(maskx,sizey,axis=0)
    masky = np.repeat(masky,sizex,axis=1)
    return masky*maskx

# Let's load the data that we need from the catalogue
def catalogue(path):
    data_gaia = np.load(path)
    return data_gaia


# Let's load the data that we need from the list of the PSFS
def list_psf(path):
    Xpsf, Ypsf = np.loadtxt(path, unpack=True, usecols=[4,5])
    return Xpsf, Ypsf


# Defining the barycenter
def barycenter(array, mask=None, x=None, y=None, subres=1):
    if isinstance(x, (np.ndarray, np.generic)) == False:
        x = (np.arange(0, array.shape[1]) + 0.5) / float(subres)
    if isinstance(y, (np.ndarray, np.generic)) == False:
        y = (np.arange(0, array.shape[0]) + 0.5) / float(subres)

    if x.ndim == 1 & y.ndim == 1:
        x, y = np.meshgrid(x, y)

    if mask is not None:
        weight = np.sum(array * mask)
        bx = np.sum(array * x * mask) / weight
        by = np.sum(array * y * mask) / weight
    else:
        tmp = np.sum(array)
        bx = np.sum(array * x) / tmp
        by = np.sum(array * y) / tmp
    return bx, by


# Defiinig the function that produces a Gaussian kernel
def gauss(xc, yc, width, size, subres=1):
    # 2D Gaussian function centered on (xc,yc)
    s = float(subres)
    p = subres * size
    (RX, RY) = np.meshgrid(np.arange(0, p) / s - xc + 0.5 / s,
                           np.arange(0, p) / s - yc + 0.5 / s)
    D2 = RX * RX + RY * RY

    return np.exp(-D2 / (2 * width * width))


# Let's define the imagette window
def window(xt, yt, sx, sy):
    i0 = np.round(xt - sx / 2)
    j0 = np.round(yt - sy / 2)
    x = xt - i0
    y = yt - j0
    return x, y, i0, j0

# Let's create a function that draws a randomly  targets from my interval
def ran_unique_int(n, interval):
    """
    Generate n unique random integer numbers in the given interval
    """
    r = np.random.randint(interval[0], interval[1], size=n)
    if interval[1] - interval[0] + 1 < n:
        print('more requested random integers than possible given the interval')
        return r
    p = 0
    while p < n:
        u = np.unique(r)
        p = len(u)
        r[0:p] = u
        if p < n:
            r[p:] = np.random.randint(interval[0], interval[1], size=n-p)
    return r


# This function plots the imagette and the PSF
def ploting_initial(rows, cols, psf, imagette, i, j):
    fig, axs = plt.subplots(rows, cols)
    axs[0].imshow(psf, origin='lower', interpolation=None)
    axs[0].set_title(i)
    axs[1].imshow(imagette, origin='lower')
    axs[1].set_title(j)
    axs[1].set_xticks([0, 1, 2, 3, 4, 5])
    axs[1].set_xticklabels([0, 1, 2, 3, 4, 5])
    fig.tight_layout()
    plt.show()


# This function plots the imagettes for both target and contaminant
def ploting_imagettes(rows, cols, ft, fc):
    fig, axs = plt.subplots(rows, cols)
    axs[0].imshow(ft, origin='lower', cmap='viridis')
    axs[0].set_title(f'For the Target')
    axs[0].set_xticks([0, 1, 2, 3, 4, 5])
    axs[0].set_xticklabels([0, 1, 2, 3, 4, 5])
    axs[1].imshow(fc, origin='lower')
    axs[1].set_title(f'For the Contaminant')
    axs[1].set_xticks([0, 1, 2, 3, 4, 5])
    axs[1].set_xticklabels([0, 1, 2, 3, 4, 5])
    fig.tight_layout()
    plt.show()


# This function plots the unsorted nsr imagette
def ploting_nsr(n, i):
    plt.imshow(n, origin='lower', cmap='viridis')
    plt.title(i)
    plt.show()


# This function plots the sorted nsr imagette
def ploting_nsr_s(n, i):
    plt.imshow(n, origin='lower', cmap='viridis')
    plt.title(i)
    plt.show()
