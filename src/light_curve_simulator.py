"""
Script that creates a lightcurve for
the nominal and extended mask
Fernando 16th Nov. 2024
"""

import numpy as np #type: ignore
import spline2dbase # type: ignore
import h5py # type: ignore


psf_file = 'PSF_Focus_0mu_0.2pxdif.npz'
psfidx = 100 # index of the PSF to be used
subres = 64 # psf sub-pixel resolution
bsres = 20 # b-spline resolution
PSFsizex = 8 # PSF size in pixels
PSFsizey = 8

nexp = 1500 #  number of exposures
driftamplitude = 1. # amplitude of the drif in pix/90 days

readout_noise = 52.  # [e-]
zero_point = 20.62 # camera zero point
integration_time = 21. # integration time [s]
focal_lenght = 247.5 # focal length [mm]
gain = 25. # full gain [e-/ADU]
camera_nuber = 24 # number of cameras
addnoise = True # add random noise or not (photon noise + detector noise)
background_level = 2500. # background level [e-/pix]
background_residual_error = 100. # background residual error [e-/pix]

imagette_size = 6 # size of the imagette
target_mag = 11 # PLATO magnitude , target
contaminant_mag = target_mag + 2 # magnitude contaminant
distance_to_contaminant_from_target = 1. # distance of the contaminant from the target [pixels]
x0,y0 = 3.1, 2.9 # target position in the imagette

# stellar variability
flux_variation_source = 'fluxtransitcorot7b_90d.hdf5' # Corot 7b depth=3.958e-4, ffluxvar
flux_contaminant_transit = 0.5/(3.958e-4) #  scaling factor for the false transit (transit in the contaminant), fluxvarsf
flux_target_transit = 0. # scaling factor for the true transit (transit in the target),  fluxvarsftf

# Name of the output file
fname = 'simulated_light_curves'

def add_poisson_noise(image):
    '''
    image_new =  add_poisson_noise(image)
    add poisson noise to an input image

    :param image: the input image
    :return: an image with poisson noise added
    '''
    image_new = np.random.poisson(lam=image.flatten(), size=(image.size))
    image_new = np.array(image_new.reshape(image.shape),dtype=image.dtype)
    return image_new

def add_gauss_noise(image,sigma):
    '''
    image_new = add_gauss_noise(image,sigma)
    add Gaussian noise to an input image

    :param image: the input image
    :param sigma: standard deviation associated with the Gaussian distribution
    :return: an image with Gaussian noise added
    '''
    noise = np.random.normal(size=image.size,scale=sigma)
    return image + noise.reshape(image.shape)


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

def NSR(M,sb, sd, sq, It, Ic):
    nsr = np.sqrt(np.sum((It + Ic + sb + sd ** 2 + sq ** 2)*M))/ np.sum(It*M)
    return nsr


def aperture_computation(ft, fc, sb, sd, sq):
    """ We are following the procedure described
    in subsection 4.6.3. of Marchiori paper
    (Binary mask)

    Args:
        ft (_float_): target flux
        fc (_float_): contaminant flux
        sb (_float_): background noise
        sd (_float_): detector noise
        sq (_float_): quantization noise
    """
    # First we compute the NSR of the system. Eq. (36) of Marchiori's paper
    nsr = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft

    # Then we flatten the nsr and also the fluxes
    nsr_1d = nsr.flatten()
    ft_1d = ft.flatten()
    fc_1d = fc.flatten()

    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting them
    nsr_1d_index = np.argsort(nsr_1d)

    # Then we obtain the target and contaminant flux for such indexes
    ft_1d = ft_1d[nsr_1d_index]
    fc_1d = fc_1d[nsr_1d_index]

    # Then we compute the aggregat noise-to-signal ratio. Eq. (37) in Marchiori's paper
    nsr_agg = np.zeros(len(ft_1d))
    for i in range(1, len(ft_1d) + 1):
        nsr_agg[i - 1] = np.sqrt(np.sum(ft_1d[:i] + fc_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i])

    # We create now a vector full with zeros
    aperture = np.zeros(36)

    # Then we create a boolean mask
    boolean_mask_for_the_aperture = nsr_1d_index[:np.argmin(nsr_agg) + 1]

    # Then we obtain the nominal mask
    aperture[boolean_mask_for_the_aperture] = 1

    # Then we reshape the aperture
    aperture = aperture.reshape((6, 6))

    return aperture


# We define here a function for obtaining an extended mask given a nominal-binary mask
def extended_binary_mask(mask, W):
    ny, nx = mask.shape
    maske = np.zeros((ny, nx))
    maske[:, :] = mask[:, :]
    for j in range(ny):
        for i in range(nx):
            if mask[j, i] > 1.0 - 1e-5:
                for k in range(-W, W + 1):
                    if (j + k >= 0) & (j + k < ny):
                        for m in range(-W, W + 1):
                            if (i + m >= 0) & (i + m < nx):
                                maske[j + k, i + m] = 1
    return maske

# Let's load the psf files
psfdata = np.load(psf_file)

psfbs = np.array(psfdata['psfbs'][psfidx])
pxc = float(psfdata['pxc'][psfidx])
pyc = float(psfdata['pyc'][psfidx])
xpsf_pix = float(psfdata['xpsf_pix'][psfidx]) # x-coordinate PSF center in pixels
ypsf_pix = float(psfdata['ypsf_pix'][psfidx]) # y-coordinate PSF center in pixels
