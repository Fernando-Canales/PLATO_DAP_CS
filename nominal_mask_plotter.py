import numpy as np # type: ignore
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant
import matplotlib.pyplot as plt # type: ignore
from astropy.io import fits # type: ignore
from astropy.io import fits as pyfits # type: ignore
import spline2dbase # type: ignore


# CONFIGURATION PARAMETERS
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues 
dataDIR= '/home/fercho/double-aperture-photometry/psf_fits/'                                                # processed PSF files
xc, yc = 3., 3.

# Target characteristics
data = np.load(cataDIR + 'SFP_DR3_20230101.npy') # star catalogue from GAIA
#fits_file_10 = fits.open(dataDIR + '6000-10452-045000.fits')
fits_file_10 = fits.open(dataDIR + '6000-18837-257329.fits')
# PSF characteristics
sizex_psf = 8
sizey_psf = 8

fits_info_10 = fits_file_10[0].data
subres = 128
bsres = 20

# Parameters for the NSR
sb = (45 * 21) # Background noise from zodiacal light in e-/px(poisson noise)times integration time (21 sec.)
sd = 50.2      # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2       # Quantization noise in units of e-rms/px

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



psf_10 = np.array(pyfits.open(dataDIR + '6000-10452-045000.fits')[0].data,dtype=np.float)
psf_10 /= psf_10.sum()
pxc_10,pyc_10 = barycenter(psf_10,subres=subres)
lx_10 = bsres*sizex_psf
ly_10 = bsres*sizey_psf
psfbs_10 = spline2dbase.Pixel2Spline(psf_10, lx_10 ,ly_10)
offx_10 = xc-pxc_10
offy_10 = yc-pyc_10
# Imagette characteristics
sizex_imagette = 6
sizey_imagette = 6



def nominal_mask_computation(ft, fc, sb, sd, sq):
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
    #fc_1d = fc.flatten()
    
    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting them
    nsr_1d_index = np.argsort(nsr_1d)
    
    # Then we obtain the target and contaminant flux for such indexes
    ft_1d = ft_1d[nsr_1d_index]
    #fc_1d = fc_1d[nsr_1d_index]
    
    # Then we compute the aggregat noise-to-signal ratio. Eq. (37) in Marchiori's paper
    nsr_agg = np.zeros(len(ft_1d))
    for i in range(1, len(ft_1d) + 1):
        nsr_agg[i - 1] = np.sqrt(np.sum(ft_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i])
    
    # We create now a vector full with zeros
    aperture = np.zeros(36)
    
    # Then we create a boolean mask
    boolean_mask_for_the_aperture = nsr_1d_index[:np.argmin(nsr_agg) + 1]
    
    # Then we obtain the nominal mask
    aperture[boolean_mask_for_the_aperture] = 1
    
    # Then we reshape the aperture
    aperture = aperture.reshape((6,6))
    
    return aperture, nsr_agg

Pmin = 10                               # minimum magnitude
Pmax = 11                               # maximum magnitude
binsize = 0.5                           # binsize around every magnitude value

mask = (data[:, 2] >= Pmin - binsize / 2.) & (data[:, 2] <= Pmax + binsize / 2.)

targets_P5 = data[mask, :]   # magnitude of every target in the current bin

x_tar = targets_P5[:, 3]     # x-coordinate in the focal plane for every TARGET
y_tar = targets_P5[:, 4]     # y-coordinate in the focal plane for every TARGET
# We convert the coordinates of the randomly chosen targets to mm to obtain the vignetting
x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)
# Let's pick one target
m_t = targets_P5[:, 2][1]                                                 # magnitude of the given TARGET
alpha = np.arctan(np.sqrt(x_tar_mm[1] ** 2 + y_tar_mm[1] ** 2) / 247.732)
# Flux of the target
# Now we define the window (imagette) and the coordinates of the TARGET inside
# Then we compute the imagette for the TARGET by integrating the b-spline decomposition of the PSF
imagette = spline2dbase.Spline2Imagette(psfbs_10, bsres, sizex_imagette, sizey_imagette, offx=offx_10, offy=offy_10)

f_ref_t = reference_flux_target(m_t) * (np.cos(alpha) ** 2) # reference flux with vignetting after integration time for the TARGET
It = f_ref_t * imagette    


nominal_mask, nsr_agg_nominal = nominal_mask_computation(ft=It, fc=0, sb=sb, sd=sd, sq=sq)
min_index = np.argmin(nsr_agg_nominal)
min_value = nsr_agg_nominal[min_index]
print('Target magnitude =', m_t)
print('Size of the mask =', nominal_mask.sum())
plt.imshow(nominal_mask, origin='lower', extent=(0,6,0,6))
plt.grid(True, linewidth=2)
plt.figure(figsize=(8, 6), dpi=300)  # 8x6 inches, 300 DPI
plt.plot(((10 ** 6) / (12 * np.sqrt(24))) * nsr_agg_nominal, 'o-')
# Add arrow and text indicating the minimum value
plt.ylabel(r'$NSR_{agg}$ over 1h and 24 cameras $[ppm hr^{\frac{1}{2}}]$')
plt.xlabel('Number of pixels composing the binary mask')
# Use tight layout to avoid cropping
#plt.tight_layout()
# Manually adjust bottom margin to prevent x-label from being cropped
plt.subplots_adjust(bottom=0.2)
plt.savefig('NSR_agg_evolution.png', dpi=300, bbox_inches='tight', pad_inches=0.2)  # Save as PNG
plt.show()
