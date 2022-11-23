import matplotlib.pyplot as plt
import spline2dbase
import scipy.signal
import h5py as h5py
import numpy as np
from fitting_psf import from_mm_2_pix, from_pix_2_mm, closest_psf, contaminants, reference_flux_target, reference_flux_contaminant
from imagette import catalogue, list_psf, barycenter, gauss, window, ploting_imagettes, ploting_nsr, ploting_nsr_s, ploting_initial
from NSR import spr_crit, NSRn, nsr1h, nsr_AGG
from pylab import *

# The first thing to do is to load the GAIA catalogue with all the stars
data = catalogue('SFP_DR3_20220831.npy')

# The second thing to do is to open the .hdf5 file that containg all the PSFs from biruni3
file_hdf5 = h5py.File('PSF.hdf5', 'r+')

# The third thing to do is to define the parameters of the Diffusion Kernel to covolve the PSFs
DifKerSize = 3  # Size [pixel]
DifKerWidth = 0.2  # width [pixel] 0.1 -> 99.99% in the central pixel ; 0.2 -> 97.5% ; 0.3 -> 81.8%% ; 0.5 -> 46.8%

# Now we build the diffusion kernel, a Gaussian function of size DifKerSize x DifKerSize centered on the middle of the central pixel
GaussKernel = gauss(math.floor(DifKerSize / 2.) + 0.5, math.floor(DifKerSize / 2.) + 0.5, DifKerWidth, DifKerSize,
                    subres=128)
GaussKernel /= GaussKernel.sum()

# Parameters for the NSR
sb = 45.  # Background noise form zodiacal light
sd = 50.2  # Overall detector noise (including readout, smearing and dark current)
sq = 7.2  # Quantization noise

# Now we save the x and y coordinates on the focal plane of all the stars in the catalogue
x_star = data[:, 3]
y_star = data[:, 4]

# Second, we load the PSF coordinates from the list file
xpsf, ypsf = list_psf('list')

# Third, we convert the PSF coordinates on the Focal Planet from mm to pixel
xpsf_pix, ypsf_pix = from_mm_2_pix(xpsf, ypsf)

# Now we can create the mask for getting just the stars from magnitude 9 to 16
mask = (data[:, 2] >= 9) & (data[:, 2] < 10)
targets = data[mask, :]

# Now we choose the random targets
np.random.seed(300)
j = np.random.randint(0, targets.shape[0] - 1, size=300)
targets = targets[j]
# Now we obtain the x and y coordinates of the targets on the focal plane
x_tar = targets[:, 3]
y_tar = targets[:, 4]

# We convert the coordinates of the randomly chosen targets to mm for obtaining the vignetting afterwards
x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)

# Now we start the loop over all the randomly chosen targets
for i in range(0, len(x_tar)):
    # First we compute the angle for obtaining the vignetting
    alpha = np.arctan(np.sqrt(x_tar_mm[i] ** 2 + y_tar_mm[i] ** 2) / 247.732)
    # Now we found the closest psf to every target (we will use this psf for the contaminants as well)
    sd = (xpsf_pix - x_tar[i]) ** 2 + (ypsf_pix - y_tar[i]) ** 2
    # We store as a string the psf index of the closest psf
    k = str(np.argmin(sd) + 1)
    #print("The closest PSF, with difussion, is the one with number:", k)
    # We select the closest psf from the .hdf5 file
    psf = np.array(file_hdf5[k])
    # Convolving now the optical PSF by the Gaussian kernel defined previously
    psf = scipy.signal.fftconvolve(psf, GaussKernel, mode='same')
    # We plot the PSF (optional)
    #plt.imshow(psf, origin='lower')
    #plt.title('PSF after adding the difussion')
    #plt.show()
    # Now we normalize the psf
    psf /= psf.sum()
    # We compute the barycenter of the psf
    pxc, pyc = barycenter(psf, subres=128)
    # Now we we define the window (imagette) and find the coordinates of the target inside of it
    x_im, y_im, i0, j0 = window(x_tar[i], y_tar[i], 6, 6)
    # Then we obtain the offset between the center of the imagette and the center of the PSF
    offx = x_im - pxc
    offy = y_im - pyc
    # Then we convert the PSFs to b-spline
    psfbs = spline2dbase.Pixel2Spline(psf, lx=20 * 8, ly=20 * 8)
    # Then we finally compute the imagette for the target by integrating the b-spline decomposition of the PSF
    imagette = spline2dbase.Spline2Imagette(psfbs, 20, 6, 6, offx=offx, offy=offy)
    #print('Target imagette:', imagette)
    # Then we can print the coordinates of the C.O.B.
    COBx, COBy = barycenter(imagette, subres=1)
    # Let's obtain the value of the reference flux after the integration time for the target star including the vignetting
    f_ref_t = reference_flux_target(targets[:, 2][i]) * (np.cos(alpha) ** 2)
    # Let's obtain the flux per pixel of the target
    I_pix_t = f_ref_t * imagette
    # Now it is time to find the contaminants sorrounding each target. We write the distance condition (10 pixels)
    dist = np.sqrt((x_star - x_tar[i]) ** 2 + (y_star - y_tar[i]) ** 2)
    # We define a useful mask now
    m = (dist > 0) & (dist <= 10)
    # We get the the index of all the contaminants now with the following line
    j = np.where(m)[0]
    #print('index of the contaminants', j)
    # Now we find the magnitude of each contaminant
    m_c = data[:, 2][j]
    #print('Magnitudes', m_c)
    #m_c = m_c[~np.isnan(m_c)]
    # Getting rid of the nans
    # Now we get the coordinates of all the contaminants for the given target as well as their total number
    x_c = x_star[j]
    y_c = y_star[j]
    n_c = len(x_c)
    # Now we find the coordinates of each contaminant inside the window
    x_c_im = x_c - i0
    y_c_im = y_c - j0
    # Now we compute the offset between the center of each contaminant and the center of the PSF
    offx_c = x_c_im - pxc
    offy_c = y_c_im - pyc
    # We define an array that will contain the 'imagettes' of every contaminant
    Ic = np.zeros((n_c, 6, 6))
    #print(Ic)
    # We define an array that will contain the 'reference fluxes' of every contaminant
    #f_ref_c = np.zeros(n_c)
    for l in range(0, n_c):
        #print('magnitude', m_c[l])
        # Then we finally compute the imagette for each contaminant by integrating the b-spline decomposition of the PSF
        Ic[l, :, :] = spline2dbase.Spline2Imagette(psfbs, 20, 6, 6, offx=offx_c[l], offy=offy_c[l])
        # Now we make sure to deal only with positive magnitudes
        if m_c[l] > 0:
            COBx_c, COBy_c = barycenter(Ic[l], subres=1)
            # Let's obtain the value of the reference flux for every contaminant star
            f_ref_c = reference_flux_contaminant(f_ref_t, m_c[l], targets[:, 2][i])
            # Let's calculate the Intensity per pixel of the imagette of every contaminant star
            Ic[l, :, :] = f_ref_c * Ic[l]

    # Now we define an array with the contribution from all the stars to each pixel
    IC = np.nansum(Ic, axis=0)
    #print('Your summed array is: {}'.format(IC))
    #print(fpixt_1d)
    #print(Ic_1d)
    # print(np.sum(Ic_1d))
    # Let's calculate the NSR of the imagette where the target is the target
    #nsr_t = NSRn(sb=sb, sd=sd, sq=sq, ft=f_pix_t, fc=IC)
    nsr_t = NSRn(45, 50.2, 7.2, I_pix_t, IC)
    # Let's plot now the NSR of the imagette of the system
    ploting_nsr(nsr_t, i='NSR of the Target')

    # Now the most intelligent thing to do right now is to convert nsr_t into a 1-D array for the rest of the script and
    # then, at the end, we convert it again into a 2-D array. Of course, it will be useful to do the same thing for the
    # fluxes per pixel of both contaminant and target
    It_1d = I_pix_t.flatten()
    Ic_1d = IC.flatten()
    nsrt_1d = nsr_t.flatten()

    # Now we sort the nsr in increasing order
    nsrt_1d_sorted = np.sort(nsrt_1d)

    # Now, in order to compute NSR_agg, I need to know the values of the flux of every pixel in the sorted NSR. So I need
    # to have an array with the sorted indexes of the original flux array, as follows:
    It_index = np.argsort(nsrt_1d)

    # We can do the same for the flux of the contaminant
    Ic_index = np.argsort(nsrt_1d)

    # Now I re-arrange the original flux array with the indexes of the sorted array
    It_1d = It_1d[It_index]
    Ic_1d = Ic_1d[Ic_index]
    # Now we compute the aggregate noise to signal ratio
    NSR_agg = nsr_AGG(It_1d, Ic_1d, 45, 50.2, 7.2)
    NSR_agg = np.array(NSR_agg)
    # Let's obtain NSR_1h
    NSR1h = nsr1h((10 ** 6) / (12 * np.sqrt(24)), NSR_agg)

    plt.plot(NSR1h, 'o-')
    plt.xlabel('Number of pixels $m$  composing the binary mask', fontsize=14)
    plt.ylabel('$ NSR_{1h}[ppm \, \sqrt{hr}] $', fontsize=14)
    plt.title('NSR evolution curve for the target mask', fontsize=14)
    # plt.savefig('binary_mask.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()

    # Now it's time to compute the mask vector
    print(f'The size of the mask is:\n', np.argmin(NSR1h) + 1, 'pixels')

    # Now it's time to compute the flux over the mask, i.e. The Flux
    f_t = It_1d[:np.argmin(NSR1h) + 1]
    f_c = Ic_1d[:np.argmin(NSR1h) + 1]

    # First we create a vector wiht only zeroes
    w_t = np.zeros(36)

    # Then we create a vector with just the amount of indexes of the mask
    #mask = fpixt_index[0:np.argmin(NSR1h) + 1]
    mask = It_index[:np.argmin(NSR1h) + 1]

    # Then we create our mask, we show the index where the mask vector has a value of 1
    w_t[mask] = 1

    print(f'This is the mask for the target:\n', mask)

    # We compute the critical SPR
    SPR_crit = spr_crit(dback=85000, nsr=NSR1h, td=4, ntr=3)

    # We create an array to store the imagettes of the contaminants which sprk > spr_crit
    n_bad = np.zeros((len(targets), 6, 6))
    # Now we can compute the SPRk for every contaminant
    for n in range(0, len(Ic)):
        # We compute the sprk of every contaminant
        sprk = np.sum(Ic[n].flatten() * w_t) / (np.sum(f_t) + np.sum(f_c))
        print('The sprk of the contaminant number', n, ' is:', sprk)
        if sprk > SPR_crit:
            print("The previous star can produce a false positive")
            # We store the imagette of the contaminant with sprk > spr_crit
            n_bad[n] = Ic[n]




