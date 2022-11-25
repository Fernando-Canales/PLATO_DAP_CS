import matplotlib.pyplot as plt
import spline2dbase
import scipy.signal
import h5py as h5py
import numpy as np
from fitting_psf import from_mm_2_pix, from_pix_2_mm, closest_psf, contaminants, reference_flux_target, \
    reference_flux_contaminant
from imagette import catalogue, list_psf, barycenter, gauss, window, ploting_imagettes, ploting_nsr, ploting_nsr_s, \
    ploting_initial
from NSR import spr_crit, aperture, NSRn, nsr_AGG
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
sb = (45. * 21)  # Background noise form zodiacal light in units of e-/px after multiplying by the integration time (poisson noise)
sd = 50.2  # Overall detector noise (including readout at beginning of life, smearing and dark current) in units of e- rms/px
sq = 7.2  # Quantization noise in units of e-rms/px

# Parameters for the eclispigin binaries
dback = 85000  # transit depth in ppm
td = 4  # transit duration in hours
ntr = 3  # number of transits

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
    s_d = (xpsf_pix - x_tar[i]) ** 2 + (ypsf_pix - y_tar[i]) ** 2
    # We store as a string the psf index of the closest psf
    k = str(np.argmin(s_d) + 1)
    # print("The closest PSF, with difussion, is the one with number:", k)
    # We select the closest psf from the .hdf5 file
    psf = np.array(file_hdf5[k])
    # Convolving now the optical PSF by the Gaussian kernel defined previously
    psf = scipy.signal.fftconvolve(psf, GaussKernel, mode='same')
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
    ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
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
    # print('index of the contaminants', j)
    # Now we find the magnitude of each contaminant
    m_c = data[:, 2][j]
    # m_c = m_c[~np.isnan(m_c)]
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
    for l in range(0, n_c):
        # print('magnitude', m_c[l])
        # Then we finally compute the imagette for each contaminant by integrating the b-spline decomposition of the PSF
        Ic[l, :, :] = spline2dbase.Spline2Imagette(psfbs, 20, 6, 6, offx=offx_c[l], offy=offy_c[l])
        # Now we make sure to deal only with stars with positive magnitudes
        if m_c[l] > 0:
            COBx_c, COBy_c = barycenter(Ic[l], subres=1)
            # Let's obtain the value of the reference flux for every contaminant star
            f_ref_c = reference_flux_contaminant(f_ref_t, m_c[l], targets[:, 2][i])
            # Let's calculate the Intensity per pixel of the imagette of every contaminant star
            Ic[l, :, :] = f_ref_c * Ic[l]

    # Now we define an array with the contribution from all the stars to each pixel
    Ic_acc = np.sum(Ic, axis=0)

    # Let's compute the aperture of the target
    NSR1h, w_t = aperture(ft=I_pix_t, fc=Ic_acc, sb=sb, sd=sd, sq=sq)

    # We compute now the flux over the mask pixels for both target and all contaminants
    f_t = I_pix_t.flatten() * w_t
    f_c = Ic_acc.flatten() * w_t

    # Let's plot the aperture
    plt.plot(NSR1h, 'o-')
    plt.xlabel('Number of pixels $m$  composing the binary mask', fontsize=14)
    plt.ylabel('$ NSR_{1h}[ppm \, \sqrt{hr}] $', fontsize=14)
    plt.title('Aperture size for the target', fontsize=14)
    # plt.savefig('binary_mask.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()

    # We compute the critical SPR now
    SPR_crit = spr_crit(dback=dback, nsr=NSR1h, td=td, ntr=ntr)

    # We compute the denominator of the sprk outside the loop
    f_tot_t = (np.sum(f_t) + np.sum(f_c))

    # We create a vector for saving the spr of every contaminant
    sprk = np.zeros(n_c)
    # Now we can compute the SPRk for every contaminant
    for n in range(0, n_c):
        # We compute the sprk of every contaminant
        sprk[n] = np.sum(Ic[n].flatten() * w_t) / f_tot_t

    # We define now a mask for getting the index of all the contaminants that fulfills spr > SPR_crit
    p = np.where(sprk > SPR_crit)[0]
    # Now we define the number of contaminant stars for which sprk > SPR_crit
    n_bad = len(p)
    # And then we compute the aperture for every of this contaminants (n_bad)
    for o in p:
        # We define the term that englobes the sigma of the target and the accumulated flux of the contaminants without the contaminant of interest
        Itc_acc = I_pix_t + Ic_acc - Ic[o]

        # Then we procedd to compute the secondary aperture
        NSR1h_c, w_c = aperture(ft=Ic[o], fc=Itc_acc, sb=sb, sd=sd, sq=sq)

        # Let's plot the aperture
        plt.plot(NSR1h_c, 'o-')
        plt.xlabel('Number of pixels $m$  composing the binary mask', fontsize=14)
        plt.ylabel('$ NSR_{1h}[ppm \, \sqrt{hr}] $', fontsize=14)
        plt.title('Aperture size for a contaminant that $SPR_{k} > SPR_{crit}$', fontsize=14)
        # plt.savefig('binary_mask.pdf', format='pdf', bbox_inches='tight', dpi=300)
        plt.show()

        # We compute the flux over the secondary mask
        f_beb = Ic_acc.flatten() * w_c
        f_t_c = I_pix_t.flatten() * w_c

        # We define the denominator of the spr_c calculation
        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))

        # We compute spr_c
        spr_c = np.sum(Itc_acc.flatten() * w_c) / f_tot_c

        # We compute now the delta_obs for the two apertures
        delta_obs_t = sprk[o] * dback
        delta_obs_c = (1 - spr_c) * dback
        print(min(NSR1h))
        print(min(NSR1h_c))
        print(m_c[o])
        print(sprk[o])
        print(spr_c)
        print(w_t)
        print(w_c)
        print('delta obs of the target is:', delta_obs_t)
        print('delta obs of the contaminant is:', delta_obs_c)

        # We compute now the statistical significances for a given transit event
        eta_t = sprk[o] * np.sqrt(td * ntr) * dback / min(NSR1h)
        eta_c = (1 - spr_c) * np.sqrt(td * ntr) * dback / min(NSR1h_c)

        print('eta of the target is:', eta_t)
        print('eta of the contaminant is:', eta_c)