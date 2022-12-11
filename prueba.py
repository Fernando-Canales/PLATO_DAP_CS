import matplotlib.pyplot as plt
import spline2dbase
import scipy.signal
import h5py as h5py
import numpy as np
from fitting_psf import from_mm_2_pix, from_pix_2_mm, closest_psf, contaminants, reference_flux_target, \
    reference_flux_contaminant
from imagette import catalogue, list_psf, barycenter, gauss, window, ran_unique_int, ploting_imagettes, ploting_nsr, \
    ploting_nsr_s, \
    ploting_initial
from NSR import spr_crit, aperture, NSRn, nsr_AGG, SPR, mask_to_bitmask, bitmask_to_mask, extended_binary_mask
from pylab import *

# The first thing to do is to load the GAIA catalogue with all the stars
data = catalogue('SFP_DR3_20220831.npy')

# The second thing to do is to open the .hdf5 file that containg all the PSFs from biruni3
file_hdf5 = h5py.File('PSF.hdf5', 'r+')

# The third thing to do is to define the parameters for the Diffusion Kernel to covolve the PSFs
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

# Parameters for the eclipsing binaries
dback = 85000  # transit depth in ppm
td = 4  # transit duration in hours
ntr = 3  # number of transits

# Define an ID for every target
ID = np.arange(0, data.shape[0])

# Now we save the x and y coordinates on the focal plane of all the stars in the catalogue
x_star = data[:, 3]
y_star = data[:, 4]

# Second, we load the PSF coordinates from the list file
xpsf, ypsf = list_psf('list')

# Third, we convert the PSF coordinates on the Focal Planet from mm to pixel
xpsf_pix, ypsf_pix = from_mm_2_pix(xpsf, ypsf)

# Define a numpy array for saving the metrics of interest (Target ID, magnitude, N_bad, etc.) (Is hard-coded now)
save_info = np.zeros((300 * 8, 17))
# The same for the extended mask
save_info_ext = np.zeros((300 * 8, 8))
# The same for bray's et al. assumption of using 2 x 2 masks
save_info_bray = np.zeros((300 * 8, 6))

# Now we choose the random targets using Réza's function
np.random.seed(300)

# We define a counter to store our data
counter = 0
# Now we can create the mask for getting only stars from P5 sample magnitude range
for i in range(7, 14):
    mask = (data[:, 2] >= i - 0.5) & (data[:, 2] < i)
    targets_P5 = data[mask, :]
    ID_target = ID[mask]

    j = ran_unique_int(n=300, interval=[0, targets_P5.shape[0] - 1])
    targets_P5 = targets_P5[j]
    ID_target = ID_target[j]
    # Now we obtain the x and y coordinates of the targets on the focal plane
    x_tar = targets_P5[:, 3]
    y_tar = targets_P5[:, 4]

    # We convert the coordinates of the randomly chosen targets to mm for obtaining the vignetting afterwards
    x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)

    # Now we start the loop over all the randomly chosen targets
    for k in range(0, len(x_tar)):
        # Define target ID
        ID_target[k]
        # First we compute the angle for obtaining the vignetting
        alpha = np.arctan(np.sqrt(x_tar_mm[k] ** 2 + y_tar_mm[k] ** 2) / 247.732)
        # Now we found the closest psf to every target (we will use this psf for the contaminants as well)
        s_d = (xpsf_pix - x_tar[k]) ** 2 + (ypsf_pix - y_tar[k]) ** 2
        # We store as a string the psf index of the closest psf
        l = str(np.argmin(s_d) + 1)
        # We select the closest psf from the .hdf5 file
        psf = np.array(file_hdf5[l])
        # Convolving now the optical PSF by the Gaussian kernel defined previously
        psf = scipy.signal.fftconvolve(psf, GaussKernel, mode='same')
        # Now we normalize the psf
        psf /= psf.sum()
        # We compute the barycenter of the psf
        pxc, pyc = barycenter(psf, subres=128)
        # Now we we define the window (imagette) and find the coordinates of the target inside of it
        x_t_im, y_t_im, i0, j0 = window(x_tar[k], y_tar[k], 6, 6)
        # Then we obtain the offset between the center of the imagette and the center of the PSF
        offx = x_t_im - pxc
        offy = y_t_im - pyc
        # Then we convert the PSFs to b-spline
        psfbs = spline2dbase.Pixel2Spline(psf, lx=20 * 8, ly=20 * 8)
        # Then we finally compute the imagette for the target by integrating the b-spline decomposition of the PSF
        imagette = spline2dbase.Spline2Imagette(psfbs, 20, 6, 6, offx=offx, offy=offy)
        # ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
        # Then we can print the coordinates of the C.O.B.
        COBx, COBy = barycenter(imagette, subres=1)
        # Let's obtain the value of the reference flux after the integration time for the target star including the vignetting
        f_ref_t = reference_flux_target(targets_P5[:, 2][k]) * (np.cos(alpha) ** 2)
        # Let's obtain the flux per pixel of the target
        It = f_ref_t * imagette
        # Now it is time to find the contaminants sorrounding each target. We write the distance condition (10 pixels)
        dist = np.sqrt((x_star - x_tar[k]) ** 2 + (y_star - y_tar[k]) ** 2)
        # We define a useful mask now
        m = (dist > 0) & (dist < 10)
        # We get the the index of all the contaminants now with the following line
        n = np.where(m)[0]
        # Now we find the magnitude of each contaminant
        m_c = data[:, 2][n]
        # m_c = m_c[~np.isnan(m_c)]
        # Getting rid of the nans
        # Now we get the coordinates of all the contaminants for the given target as well as their total number
        x_c = x_star[n]
        y_c = y_star[n]
        # We define now the number of contaminants
        n_c = len(x_c)
        # Now we find the coordinates of each contaminant inside the window
        x_c_im = x_c - i0
        y_c_im = y_c - j0
        # Now we compute the offset between the center of each contaminant and the center of the PSF
        offx_c = x_c_im - pxc
        offy_c = y_c_im - pyc
        # We define an array that will contain the 'imagettes' of every contaminant
        Ic = np.zeros((n_c, 6, 6))
        for o in range(0, n_c):
            # Then we finally compute the imagette for each contaminant by integrating the b-spline decomposition of the PSF
            Ic[o, :, :] = spline2dbase.Spline2Imagette(psfbs, 20, 6, 6, offx=offx_c[o], offy=offy_c[o])
            # Now we make sure to deal only with stars with positive magnitudes
            if m_c[o] > 0:
                COBx_c, COBy_c = barycenter(Ic[o], subres=1)
                # Let's obtain the value of the reference flux for every contaminant star
                f_ref_c = reference_flux_contaminant(f_ref_t, m_c[o], targets_P5[:, 2][k])
                # Let's calculate the Intensity per pixel of the imagette of every contaminant star
                Ic[o, :, :] = f_ref_c * Ic[o]

        # Now we define an array with the contribution from all the stars to each pixel
        Ic_acc = np.sum(Ic, axis=0)

        # Let's compute the aperture of the target
        NSR1h, w_t = aperture(ft=It, fc=Ic_acc, sb=sb, sd=sd, sq=sq)

        # Now we store the nominal mask into a mask_key
        w_t_key = mask_to_bitmask(w_t)

########################################################################################################################
#                                   NOW THE EXTENDED MASK METHOD                                                       #
########################################################################################################################

        # First we create the extended mask given the nominal mask
        w_ext = extended_binary_mask(w_t, W=1)

        # Now we compute all the metrics associated with this mask. Let's begin with the NSR
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext)) / np.sum(It * w_ext)
        NSR_ext_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext

        # Then compute the critical SPR
        SPR_crit_ext = spr_crit(dback=dback, nsr=NSR_ext_1h, td=td, ntr=ntr)

        # Then we compute the sprk over the extended mask for all the contaminants for a this target
        sprk_ext, sprk_max_ext, SPR_tot_ext, n_bad_ext = SPR(SPR_crit=SPR_crit_ext, n_c=n_c, f_contaminant=Ic,
                                                             f_tot=(It + Ic_acc), w=w_ext)

        # Now we compute the observed transit depth given the maximum value of sprk
        delta_obs_ext = sprk_max_ext * dback

        # And now we compute the statistical significance of an Earth-like planet over the extended mask
        eta_ext = sprk_max_ext * np.sqrt(td * ntr) * dback / NSR_ext_1h

########################################################################################################################
#                                     END OF THE EXTENDED MASK METHOD                                                  #
########################################################################################################################

########################################################################################################################
#                                   TESTING  J.C. Bray et al's ASSUMPTION OF A 2 x 2 MASK                              #
########################################################################################################################

        # The mask has to contain the 4 pixels around the center,
        w_bray = np.zeros((6, 6))
        w_bray[2:4, 2:4] = 1

        NSR_bray = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_bray)) / np.sum(It * w_bray)
        NSR_bray_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_bray

        # We compute the critical SPR now
        SPR_crit_bray = spr_crit(dback=dback, nsr=NSR_bray_1h, td=td, ntr=ntr)

        sprk_bray, sprk_max_bray, SPR_tot_bray, n_bad_bray = SPR(SPR_crit=SPR_crit_bray, n_c=n_c, f_contaminant=Ic,
                                                                 f_tot=(It + Ic_acc), w=w_bray)

########################################################################################################################
#                                          END OF TESTING Bray et al's ASSUMPTION                                      #
########################################################################################################################

        # Now in this part of the code we present the calculations for the sprk of every contaminant as well as the
        # calculation of the SPR_crit.

        # We compute the critical SPR now
        SPR_crit = spr_crit(dback=dback, nsr=NSR1h, td=td, ntr=ntr)

        sprk, sprk_max, SPR_tot, n_bad = SPR(SPR_crit=SPR_crit, n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)

        # Now we have to know for which contaminant corresponds the highest sprk value we just found
        ind_sprk = np.argmax(sprk)
        Ic_max = Ic[ind_sprk]

        # Now we compute the secondary aperture for this contamninant with the highets sprk

        # We define the term that englobes the sigma of the target and the accumulated flux of the contaminants without
        # the contaminant of interest
        Itc_acc = It + Ic_acc - Ic_max

        # Then we procedd to compute the secondary aperture
        NSR1h_c, w_c = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq)

        # Now we store this secondary mask in a mask key
        w_c_key = mask_to_bitmask(w_c)

        # We compute the flux over the secondary mask
        f_beb = Ic_acc * w_c
        f_t_c = It * w_c

        # We define the denominator of the spr_c calculation
        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))

        # We compute spr_c
        spr_c = np.sum(Itc_acc * w_c) / f_tot_c

        # We compute now the delta_obs for the two apertures
        delta_obs_t = sprk_max * dback
        delta_obs_c = (1 - spr_c) * dback

        # We compute now the statistical significances for a given transit event
        eta_t = sprk_max * np.sqrt(td * ntr) * dback / NSR1h
        eta_c = (1 - spr_c) * np.sqrt(td * ntr) * dback / NSR1h_c

        print('Delta P is:', m_c[ind_sprk] - targets_P5[:, 2][k])
        print('Distance between the target and contaminant',
              (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        # if eta_t > eta_c:
        print('NSR_T is:', NSR1h)
        print('NSR_c is:', NSR1h_c)
        print('eta_t is:', eta_t)
        print('eta_c is:', eta_c)
        print('spr_t', sprk_max)
        print('spr_c', spr_c)

        # Now we need to compute the efficiency of the extended mask, this is donde by computing the ratio of the number
        # of false positives given by the extended mask such that eta_ext > eta_t over the number of false positives
        # given by the nominal mask.

        # The number of false positives given by the extended mask such that etat_ext > eta_t is given by
        n_eff_ext = len(np.where((sprk_ext > SPR_crit_ext) & (sprk > SPR_crit) & (sprk_ext > sprk))[0])


        save_info[counter, :] = [ID_target[k], targets_P5[:, 2][k], l, n_c, w_t_key, NSR1h, n_bad, SPR_crit,
                                 m_c[ind_sprk], sprk_max, NSR1h_c, w_c_key, SPR_tot, eta_t, eta_c, eta_ext, n_eff_ext]
        save_info_ext[counter, :] = [ID_target[k], targets_P5[:, 2][k], n_c, NSR_ext_1h, n_bad_ext, SPR_crit_ext,
                                     SPR_tot_ext, eta_ext]
        save_info_bray[counter, :] = [ID_target[k], targets_P5[:, 2][k], n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray]
        counter = counter + 1

save_info = save_info[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]
# Now it is time to save the metrics of interest into a.npy file
np.save('targets_P5.npy', save_info)
np.save('targets_P5_extended.npy', save_info_ext)
np.save('targets_P5_bray.npy', save_info_bray)
