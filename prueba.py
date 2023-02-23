# noinspection PyUnresolvedReferences
import numpy as np
import spline2dbase
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant
from imagette import catalogue, barycenter, window, ran_unique_int, centroid_shift
from NSR import spr_crit, aperture, SPR, mask_to_bitmask, extended_binary_mask
from pylab import *

# First we establish all relevant paths: star catalogue, PSF data file and storage directory
cataloguesDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/catalogues_stars/'
PSFfile = 'PSF.npz'
DIRout = 'test_results/'

# Now we call the catalogue, the PSF file and the transit depths and transit durations catalogue for the eclipsing
# binaries
data = np.load(cataloguesDIR + 'SFP_DR3_20220831.npy')
psfdata = np.load(PSFfile)
# delta_back, transit_dur = np.loadtxt(cataloguesDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt', unpack=True,
# usecols=[0, 1])

# Parameters for the imagette and PSF
size_im_x = 6  # size of the imagette (x-direction)
size_im_y = 6  # size of the imagette (y-direction)
subres = 128  # resolution of the PSF
bsres = 20  # resolution of the b-spline decomposition of the PSF

# Parameters for the NSR
sb = (45. * 21)  # Background noise from zodiacal light in units of e-/px(poisson noise)times integration time (21 sec.)
sd = 50.2  # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2  # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries
dback = 85000  # transit depth in ppm
td = 4  # transit duration in hours
ntr = 3  # number of transits in one hour

# Define an ID for every target
ID = np.arange(0, data.shape[0])
# Define the number of targets
n_tar = 300
# Now we save the x and y coordinates on the focal plane for all the stars in the catalogue
x_star = data[:, 3]
y_star = data[:, 4]

# We load the PSF data after running Réza's script (process_psf.py)
psfbs = psfdata['psfbs']
pxc = psfdata['pxc']
pyc = psfdata['pyc']
xpsf_pix = psfdata['xpsf_pix']
ypsf_pix = psfdata['ypsf_pix']

# Now we use a seed for obtaining the same number of targets and dback and td values
np.random.seed(n_tar)

# We set the minimum value of magnitude
Pmin = 8
# We set the maximum value of magnitude
Pmax = 13
# We se the binsize value
binsize = 0.5
# We set the number of intervals
nP = int((Pmax - Pmin) / binsize + 1)

# Define a numpy array for saving the metrics of interest (Target ID, magnitude, N_bad, etc.)
save_info = np.zeros((n_tar * nP, 17))
# The same for the secondary/contaminant mask
save_info_contaminant = np.zeros((n_tar * nP, 20))
# The same for the extended mask
save_info_ext = np.zeros((n_tar * nP, 28))
# The same for bray's et al. assumption of using 2 x 2 masks
save_info_bray = np.zeros((n_tar * nP, 8))

# We create a numpy array for saving the number of stars per magnitude bin
n_star_p_bin = np.zeros(nP)

# We define this counter in order to store our data
counter = 0

# Now we can create the mask for getting only stars from P5 sample magnitude range
for i in range(nP):
    Pi = Pmin + i * binsize
    mask = (data[:, 2] >= Pi - binsize / 2.) & (data[:, 2] <= Pi + binsize / 2.)
    n_star_p_bin[i] = mask.sum()
    targets_P5 = data[mask, :]
    ID_target = ID[mask]

    j = ran_unique_int(n=n_tar, interval=[0, targets_P5.shape[0] - 1])
    targets_P5 = targets_P5[j]
    ID_target = ID_target[j]
    # Now we obtain the x and y coordinates of the targets on the focal plane
    x_tar = targets_P5[:, 3]
    y_tar = targets_P5[:, 4]

    # We convert the coordinates of the randomly chosen targets to mm for obtaining the vignetting afterwards
    x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)
    print(
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    print('Beginning the calculations for the targets of magnitude', Pi, '\n')
    print(
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    # Now we start the loop over all the randomly chosen targets
    for k in range(0, len(x_tar)):
        # Define target ID
        ID_target[k]
        # First we compute the angle for obtaining the vignetting
        alpha = np.arctan(np.sqrt(x_tar_mm[k] ** 2 + y_tar_mm[k] ** 2) / 247.732)
        # Now we found the closest psf to every target (we will use this psf for the contaminants as well)
        s_d = (xpsf_pix - x_tar[k]) ** 2 + (ypsf_pix - y_tar[k]) ** 2
        # Now we obtain the index of the closest PSF
        psf_idx = np.argmin(s_d)
        # Now we compute the barycenter coordinates of the PSF
        pxc[psf_idx]
        pyc[psf_idx]
        # Now we get the b-spline decomposition of the closest PSF
        psfbs[psf_idx]
        # Now we define the window (imagette) and find the coordinates of the target inside of it
        x_t_im, y_t_im, i0, j0 = window(x_tar[k], y_tar[k], 6, 6)
        # Then we obtain the offset between the center of the imagette and the center of the PSF
        offx = x_t_im - pxc[psf_idx]
        offy = y_t_im - pyc[psf_idx]
        # Then we finally compute the imagette for the target by integrating the b-spline decomposition of the PSF
        imagette = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y, offx=offx, offy=offy)
        # ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
        # Then we can print the coordinates of the C.O.B.
        COBx, COBy = barycenter(imagette, subres=1)
        # Now we define a variable for the magnitude of every target we are going to analyze
        m_t = targets_P5[:, 2][k]
        # Let's obtain the value of the reference flux after the integration time for the target star with vignetting
        f_ref_t = reference_flux_target(m_t) * (np.cos(alpha) ** 2)
        # Let's obtain the flux per pixel of the target:;
        It = f_ref_t * imagette
        # Now it is time to find all the contaminants surrounding each target. We put the distance condition (10 pixels)
        dist = np.sqrt((x_star - x_tar[k]) ** 2 + (y_star - y_tar[k]) ** 2)
        # Now we define a useful mask
        m = (dist > 0) & (dist < 10)
        # Now we get the index of all the contaminants that fulfills the mask requirement
        n = np.where(m)[0]
        # Now we find the magnitude of each contaminant of the given target
        m_c = data[:, 2][n]
        # Now we get the coordinates of all the contaminants for the given target as well as their total number
        x_c = x_star[n]
        y_c = y_star[n]
        # We define now the number of contaminants for a given target
        n_c = len(x_c)
        # Now we find the coordinates of each contaminant inside the window
        x_c_im = x_c - i0
        y_c_im = y_c - j0
        # Now we compute the offset between the center of each contaminant and the center of the PSF
        offx_c = x_c_im - pxc[psf_idx]
        offy_c = y_c_im - pyc[psf_idx]
        # We define an array that will contain the 'imagettes' of every contaminant
        Ic = np.zeros((n_c, 6, 6))
        for o in range(1, n_c + 1):
            # Now we make sure to deal only with stars with positive magnitudes
            if m_c[o - 1] > 0:
                # Then we compute the imagette for each contaminant by integrating the b-spline decomposition of the PSF
                Ic[o - 1, :, :] = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y,
                                                               offx=offx_c[o - 1], offy=offy_c[o - 1])
                COBx_c, COBy_c = barycenter(Ic[o - 1], subres=1)
                # Let's obtain the value of the reference flux for every contaminant star
                f_ref_c = reference_flux_contaminant(f_ref_t, m_c[o - 1], m_t)
                # Let's calculate the Intensity per pixel of the imagette of every contaminant star
                Ic[o - 1, :, :] = f_ref_c * Ic[o - 1]

        # Now we define an array with the contribution from all the stars to each pixel
        Ic_acc = np.sum(Ic, axis=0)

        # Now we define the total flux (target and all contaminants)
        f_tot = It + Ic_acc
        # Let's compute the aperture of the target
        NSR1h, w_t = aperture(ft=It, fc=Ic_acc, sb=sb, sd=sd, sq=sq)

        # Now we store the nominal mask into a mask_key
        w_t_key = mask_to_bitmask(w_t)

        # Now we store the size of each nominal mask
        w_t_size = np.count_nonzero(w_t)

        # Now we present the calculation for the sprk of every contaminant as well as the one from SPR_crit.

        # sprk, sprk_max, SPR_tot, n_bad = SPR(SPR_crit=SPR_crit, n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)
        sprk, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_t)
        # Now we choose randomly a value for dback and for td
        # dback_index = np.random.randint(0, len(delta_back))
        # td_index = np.random.randint(0, len(transit_dur))
        # We obtain the values of dback and td from the catalogue
        # dback = delta_back[dback_index]
        # td = transit_dur[td_index]
        # We compute the critical SPR now
        SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr)

        # Now we compute the number of contaminant stars above SPR_crit (i.e. N_bad)
        n_bad = np.sum(sprk > SPR_crit)
        # Now we get the index of the contaminant star with the highest sprk value with respect to the nominal mask
        ind_sprk = np.argmax(sprk)
        # Now we select the (intensity per pixel array) imagette of the contaminant star with the highest sprk value
        # with respect to the nominal mask
        Ic_max = Ic[ind_sprk]
        # Now we obtain the magnitude for the contaminant with the highest sprk value with respect to the nominal mask
        m_c_bad = m_c[ind_sprk]
        # Now we find the distance between the given target and the contaminant with the highest value of sprk
        dist_bad = (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2
        # Now we define a term that contains the flux of the targets and their respective contaminants except the
        # one with the highest value of SPRk (i.e. the contaminant of interest)
        Itc_acc = It + Ic_acc - Ic_max

        # Then we proceed to compute the secondary aperture
        NSR1h_c, w_c = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq)

        # Now we store each secondary mask in a mask key
        w_c_key = mask_to_bitmask(w_c)

        # Now we store the size of each secondary mask
        w_c_size = np.count_nonzero(w_c)

        # Now, as a test, we create an 'extended' secondary mask; just to see how it goes
        w_c_2 = extended_binary_mask(w_c, W=1)

        # Now we store each secondary mask in a mask key
        w_c_key_2 = mask_to_bitmask(w_c)

        # Now we store the size of each secondary mask
        w_c_size_2 = np.count_nonzero(w_c)

        # Now we compute the NSR of this extended secondary mask
        NSR_c_2 = np.sqrt(np.sum((Ic_max + Itc_acc + sb + sd ** 2 + sq ** 2) * w_c_2)) / np.sum(Ic_max * w_c_2)
        # Now we compute the NSR over 1hr
        NSR1h_c_2 = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_c_2
        print(type(NSR_c_2))
        print(type(NSR1h_c_2))

        # We compute the sprk with respect to the secondary mask
        sprk_sec, SPR_tot_sec = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_c)
        sprk_sec_2, SPR_tot_sec_2 = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_c_2)

        # We compute the flux over the secondary mask and the extended secondary mask
        f_beb = Ic_acc * w_c
        f_t_c = It * w_c
        f_beb_2 = Ic_acc * w_c_2
        f_t_c_2 = It * w_c_2

        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))
        f_tot_c_2 = (np.sum(f_t_c_2) + np.sum(f_beb_2))

        # We compute spr_tot_c, that is, the expression given in Marchiori presentation for PLATO week #8
        spr_tot_c = np.sum(Itc_acc * w_c) / f_tot_c
        spr_tot_c_2 = np.sum(Itc_acc * w_c_2) / f_tot_c_2

        # We compute now the delta_obs for the two apertures
        delta_obs_t = sprk[ind_sprk] * dback
        delta_obs_c = (1 - spr_tot_c) * dback
        delta_obs_c_2 = (1 - spr_tot_c_2) * dback

        # We compute now the statistical significances for a given transit event
        eta_t = sprk[ind_sprk] * np.sqrt(td * ntr) * dback / (NSR1h * (1 - SPR_tot))
        eta_c = np.sqrt(td * ntr) * dback / NSR1h_c
        eta_c_2 = np.sqrt(td * ntr) * dback / NSR1h_c_2

        # -------------------------------------------NOMINAL COB-------------------------------------------------------#
        eta_cob, sigma_1_24, abs_cob = centroid_shift(w=w_t, Ik=Ic_max, I_t=It, I_contaminants=Ic_acc,
                                                      sprk=sprk[ind_sprk], dback=dback, sb=sb, sd=sd, sq=sq, td=td,
                                                      ntr=ntr)
        # -------------------------------------------NOMINAL COB-------------------------------------------------------#

        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        eta_cob_c, sigma_1_24_c, abs_cob_c = centroid_shift(w=w_c, Ik=Ic_max, I_t=It, I_contaminants=Ic_acc,
                                                            sprk=sprk_sec[ind_sprk], dback=dback, sb=sb, sd=sd, sq=sq,
                                                            td=td, ntr=ntr)

        eta_cob_c_2, sigma_1_24_c_2, abs_cob_c_2 = centroid_shift(w=w_c_2, Ik=Ic_max, I_t=It, I_contaminants=Ic_acc,
                                                                  sprk=sprk_sec_2[ind_sprk], dback=dback, sb=sb, sd=sd,
                                                                  sq=sq, td=td, ntr=ntr)
        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        ################################################################################################################
        #                                   NOW THE EXTENDED MASK METHOD                                               #
        ################################################################################################################

        """
        Now we create different extended masks given the nominal mask
        """
        # We create the first extended mask by adding 1 pixel to the nominal mask
        w_ext = extended_binary_mask(w_t, W=1)
        # We create the second extended mask by adding 2 pixels to the nominal mask
        w_ext_2 = extended_binary_mask(w_t, W=2)
        # We create the second extended mask by adding 2 pixels to the nominal mask
        w_ext_3 = extended_binary_mask(w_t, W=3)

        # Now we store each extended mask in a mask key
        w_ext_key = mask_to_bitmask(w_ext)
        w_ext_key_2 = mask_to_bitmask(w_ext_2)
        w_ext_key_3 = mask_to_bitmask(w_ext_3)

        # Now we store the size of each extended mask
        w_ext_size = np.count_nonzero(w_ext)
        w_ext_size_2 = np.count_nonzero(w_ext_2)
        w_ext_size_3 = np.count_nonzero(w_ext_3)

        # Now we compute all the metrics associated with every extended mask. Let's begin with the NSR
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext)) / np.sum(It * w_ext)
        NSR_ext_2 = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext_2)) / np.sum(It * w_ext_2)
        NSR_ext_3 = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext_3)) / np.sum(It * w_ext_3)
        NSR_ext_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext
        NSR_ext_1h_2 = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext_2
        NSR_ext_1h_3 = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext_3

        # Then we compute the sprk over each extended mask for all the contaminants for a given target
        sprk_ext, SPR_tot_ext = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_ext)
        sprk_ext_2, SPR_tot_ext_2 = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_ext_2)
        sprk_ext_3, SPR_tot_ext_3 = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_ext_3)

        # Then compute the critical SPR of the extended mask (SPR_crit_ext)
        SPR_crit_ext = spr_crit(dback=dback, SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h, td=td, ntr=ntr)
        SPR_crit_ext_2 = spr_crit(dback=dback, SPR_tot=SPR_tot_ext_2, nsr=NSR_ext_1h_2, td=td, ntr=ntr)
        SPR_crit_ext_3 = spr_crit(dback=dback, SPR_tot=SPR_tot_ext_3, nsr=NSR_ext_1h_3, td=td, ntr=ntr)

        # Now we compute the observed transit depth given the contaminant with the highest sprk on the nominal mask
        delta_obs_ext = sprk_ext[ind_sprk] * dback
        delta_obs_ext_2 = sprk_ext_2[ind_sprk] * dback
        delta_obs_ext_3 = sprk_ext_3[ind_sprk] * dback

        # Now we compute the statistical significance of the signal over the extended mask
        eta_ext = sprk_ext[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h
        eta_ext_2 = sprk_ext_2[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h_2
        eta_ext_3 = sprk_ext_3[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h_3

        # ------------------------------------------EXTENDED COB--------------------------------------------------------#
        eta_cob_ext, sigma_1_24_ext, abs_cob_ext = centroid_shift(w=w_ext, Ik=Ic_max, I_t=It, I_contaminants=Ic_acc,
                                                                  sprk=sprk_ext[ind_sprk], dback=dback, sb=sb, sd=sd,
                                                                  sq=sq, td=td, ntr=ntr)

        eta_cob_ext_2, sigma_1_24_ext_2, abs_cob_ext_2 = centroid_shift(w=w_ext_2, Ik=Ic_max, I_t=It,
                                                                        I_contaminants=Ic_acc,
                                                                        sprk=sprk_ext_2[ind_sprk],
                                                                        dback=dback, sb=sb, sd=sd, sq=sq, td=td,
                                                                        ntr=ntr)

        eta_cob_ext_3, sigma_1_24_ext_3, abs_cob_ext_3 = centroid_shift(w=w_ext_3, Ik=Ic_max, I_t=It,
                                                                        I_contaminants=Ic_acc,
                                                                        sprk=sprk_ext_3[ind_sprk],
                                                                        dback=dback, sb=sb, sd=sd, sq=sq, td=td,
                                                                        ntr=ntr)
        # ------------------------------------------EXTENDED COB--------------------------------------------------------#

        ################################################################################################################
        #                                     END OF THE EXTENDED MASK METHOD                                          #
        ################################################################################################################
        ################################################################################################################
        #                                   TESTING  J.C. Bray et al.'s ASSUMPTION OF A 2 x 2 MASK                     #
        ################################################################################################################

        # The mask has to contain the 4 pixels around the center,
        w_bray = np.zeros((6, 6))
        w_bray[2:4, 2:4] = 1

        NSR_bray = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_bray)) / np.sum(It * w_bray)
        NSR_bray_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_bray

        sprk_bray, SPR_tot_bray = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_bray)

        SPR_crit_bray = spr_crit(dback=dback, SPR_tot=SPR_tot_bray, nsr=NSR_bray_1h, td=td, ntr=ntr)

        # Now we compute the number of contaminant stars above SPR_crit (i.e. N_bad)
        n_bad_bray = np.sum(sprk_bray > SPR_crit_bray)

        ################################################################################################################
        #                                          END OF TESTING Bray et al.'s ASSUMPTION                             #
        ################################################################################################################
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              '++++++')
        print('Target ID:', ID_target[k])
        print('Total number of contaminant stars:', n_c)
        print('Number of contaminant stars with sprk values above SPR_crit:', n_bad)
        print('Magnitude difference between the target and the contaminant with the highest value of sprK:',
              m_c[ind_sprk] - m_t)
        print('The distance between the target and the contaminant with the highest value of sprk is:',
              (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        print('Transit depth of the eclipsing binary with the highest value of sprk:', dback, 'ppm')
        print('Transit duration of the eclipsing binary with the highest value of sprk:', td, 'hours')
        print('Nominal COB shift significance:', abs_cob)
        print('Secondary COB shift significance:', abs_cob_c)
        print('Extended COB shift significance:', abs_cob_ext)
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              '++++++\n')
        # The number of false positives given by the extended mask such that eta_ext > eta_t is given by
        # n_eff_ext = len(np.where((sprk_ext > SPR_crit_ext) & (sprk[ind_sprk] > SPR_crit) & (sprk_ext > sprk[ind_sprk])
        # )[0])

        save_info[counter, :] = [ID_target[k], m_t, n_c, m_c_bad, dist_bad, w_t_key, w_t_size, NSR1h, n_bad, SPR_crit,
                                 sprk[ind_sprk], SPR_tot, eta_t, delta_obs_t, abs_cob, eta_cob, sigma_1_24]

        save_info_contaminant[counter, :] = [ID_target[k], m_t, w_c_key, w_c_size, NSR1h_c, spr_tot_c, eta_c,
                                             delta_obs_c, abs_cob_c, eta_cob_c, sigma_1_24_c, w_c_key_2, w_c_size_2,
                                             NSR1h_c_2, spr_tot_c_2, eta_c_2, delta_obs_c_2, abs_cob_c_2, eta_cob_c_2,
                                             sigma_1_24_c_2]

        save_info_ext[counter, :] = [ID_target[k], m_t, w_ext_key, w_ext_size, NSR_ext_1h, sprk_ext[ind_sprk],
                                     SPR_crit_ext, eta_ext, delta_obs_ext, abs_cob_ext, eta_cob_ext, sigma_1_24_ext,
                                     w_ext_key_2, w_ext_size_2, NSR_ext_1h_2, eta_ext_2, delta_obs_ext_2, abs_cob_ext_2,
                                     eta_cob_ext_2, sigma_1_24_ext_2, w_ext_key_3, w_ext_size_3, NSR_ext_1h_3,
                                     eta_ext_3, delta_obs_ext_3, abs_cob_ext_3, eta_cob_ext_3, sigma_1_24_ext_3]

        save_info_bray[counter, :] = [ID_target[k], m_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray,
                                      sprk_bray[ind_sprk], SPR_tot_bray]
        counter = counter + 1

save_info = save_info[0:counter]
save_info_contaminant = save_info_contaminant[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]

# Now it is time to save the metrics of interest into several .npy files
np.save(DIRout + 'targets_P5.npy', save_info)
np.save(DIRout + 'targets_P5_contaminant.npy', save_info_contaminant)
np.save(DIRout + 'targets_P5_extended.npy', save_info_ext)
np.save(DIRout + 'targets_P5_bray.npy', save_info_bray)

fout = open(DIRout + 'star_count.txt', 'w')
for i in range(nP):
    Pi = Pmin + i * binsize
    fout.write('%.2f %i\n' % (Pi, n_star_p_bin[i]))
fout.close()
