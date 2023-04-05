import numpy as np
import spline2dbase
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant
from imagette import barycenter, window, ran_unique_int, centroid_shift
from NSR import spr_crit, aperture, SPR, mask_to_bitmask, extended_binary_mask
from pylab import *

# -----------------------------------------------
# CONFIGURATION PARAMETERS

# Parameters relative to all the relevant paths
cataDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/catalogues_stars/' # directory with all star catalogues 
PSFfile = 'PSF.npz'                                                        # processed PSF files
DIRout = 'test_results/'                                                   # storage directory

# Parameters for the imagette and PSF decomposition
size_im_x = 6  # size of the imagette (x-direction)
size_im_y = 6  # size of the imagette (y-direction)
subres = 128   # resolution of the PSF
bsres = 20     # resolution of the b-spline decomposition of the PSF

# Parameters for the NSR
sb = (45. * 21)  # Background noise from zodiacal light in e-/px(poisson noise)times integration time (21 sec.)
sd = 50.2        # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2         # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
dback = 85000  # transit depth in ppm
td = 4         # transit duration in hours
ntr = 3        # number of transits in one hour

# Parameters for the magnitude intervals
Pmin = 8                                # minimum magnitude
Pmax = 13                               # maximum magnitude
binsize = 0.5                           # binsize around every magnitude value
nP = int((Pmax - Pmin) / binsize + 1)   # number of bins

# -----------------------------------------------

data = np.load(cataDIR + 'SFP_DR3_20230101.npy') # star catalogue from GAIA
psfdata = np.load(PSFfile)                       # processed PSFs 
# del_back, tr_dur = np.loadtxt(cataDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt', unpack=True, usecols=[0, 1]) # transit depth and duration from Kepler Eclipsing Binary Catalogue

ID = np.arange(0, data.shape[0]) # ID for every star in the catalogue
n_tar = 300                      # number of targets per magnitude interval
x_star = data[:, 3]              # x-coordinate on the focal plane for every star in the catalogue
y_star = data[:, 4]              # y-coordinate on the focal plane for every star in the catalogue
psfbs = psfdata['psfbs']         # PSF b-spline decomposition
pxc = psfdata['pxc']             # x-coordinate of the PSF b-spline decomposition centroid
pyc = psfdata['pyc']             # y-cooridnate of the PSF b-spline decomposition centroid
xpsf_pix = psfdata['xpsf_pix']   # x-coordinate of the PSF in pixel
ypsf_pix = psfdata['ypsf_pix']   # y-coordinate of the PSF in pixel


# Now we use a seed for obtaining the same number of targets and dback and td values
np.random.seed(n_tar)

save_info = np.zeros((n_tar * nP, 19))             # numpy array to store the metrics of interest for the nominal mask
save_info_contaminant = np.zeros((n_tar * nP, 13)) # numpy array to store the metrics of interest for the secondary mask
save_info_ext = np.zeros((n_tar * nP, 15))         # numpy array to store the metrics of interest for the extended mask
save_info_bray = np.zeros((n_tar * nP, 8))         # numpy array to store the metrics of interest for Bray's 2 x 2 mask
n_star_p_bin = np.zeros(nP)                        # numpy array to store the number of stars per bin

# We define a counter before the main loop
counter = 0

# We start the main for cycle, to loop over each bin
for i in range(nP):
    Pi = Pmin + i * binsize
    # Now we define a mask to get only the targets whithin the i-th magnitude bin
    mask = (data[:, 2] >= Pi - binsize / 2.) & (data[:, 2] <= Pi + binsize / 2.)
    n_star_p_bin[i] = mask.sum() # number of targets in the current bin
    targets_P5 = data[mask, :]   # magnitude of every target in the current bin
    ID_target = ID[mask]         # ID of every target in the current bin
    # Now we choose (randomly) only 'n_tar' targets from every bin. These are our new "TARGETS".
    j = ran_unique_int(n=n_tar, interval=[0, targets_P5.shape[0] - 1])
    targets_P5 = targets_P5[j]   # magnitude of every TARGET
    ID_target = ID_target[j]     # ID of every new TARGET
    n_t = len(j)                 # number of new TARGET
    x_tar = targets_P5[:, 3]     # x-coordinate in the focal plane for every TARGET
    y_tar = targets_P5[:, 4]     # y-coordinate in the focal plane for every TARGET
    # We convert the coordinates of the randomly chosen targets to mm to obtain the vignetting
    x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)
    
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    print('Beginning Calculations for Targets of Magnitude:', Pi, '\n')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    
    def process_target(k): # main processing routine (works on a given TARGET)
        
        ID_t = ID_target[k]                                                       # TARGET ID
        m_t = targets_P5[:, 2][k]                                                 # magnitude of the given TARGET
        alpha = np.arctan(np.sqrt(x_tar_mm[k] ** 2 + y_tar_mm[k] ** 2) / 247.732) # angle to compute the vignetting
        s_d = (xpsf_pix - x_tar[k]) ** 2 + (ypsf_pix - y_tar[k]) ** 2             # distance btwn all the PSFs and the given TARGET
        psf_idx = np.argmin(s_d)                                                  # index of the closest PSF
        pxc[psf_idx]                                                              # x-coordinate of the barycenter of the closest PSF
        pyc[psf_idx]                                                              # y-coordinate of the barycenter of the closest PSF
        psfbs[psf_idx]                                                            # b-spline decomposition of the closest PSF

        # Now we define the window (imagette) and the coordinates of the TARGET inside
        x_t_im, y_t_im, i0, j0 = window(x_tar[k], y_tar[k], 6, 6)
        offx = x_t_im - pxc[psf_idx] # x-coordinate of the offset between the center of the imagette and the center of the PSF
        offy = y_t_im - pyc[psf_idx] # y-coordinate of the offset between the center of the imagette and the center of the PSF
        # Then we compute the imagette for the TARGET by integrating the b-spline decomposition of the PSF
        imagette = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y, offx=offx, offy=offy)
        # ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
        #COBx, COBy = barycenter(imagette, subres=1)
        f_ref_t = reference_flux_target(m_t) * (np.cos(alpha) ** 2) # reference flux with vignetting after integration time for the TARGET
        It = f_ref_t * imagette                                     # Intensity of the TARGET (flux per pixel)
        
        # Now we find all the contaminants surrounding each TARGET. We put a distance condition (10 pixels)
        dist = np.sqrt((x_star - x_tar[k]) ** 2 + (y_star - y_tar[k]) ** 2)

        m = (dist > 0) & (dist < 10)   # mask containing the distance condition for contaminants
        n = np.where(m)[0]             # index of the contaminants
        m_c = data[:, 2][n]            # magnitude of each contaminant star for a given TARGET
        x_c = x_star[n]                # x-coordinate of a given contaminant in the focal plane
        y_c = y_star[n]                # y-coordinate of a given contaminant in the focal plane
        n_c = len(x_c)                 # number of contaminants for a given TARGET
        x_c_im = x_c - i0              # x-coordinate of a given contaminant inside the imagette (window)
        y_c_im = y_c - j0              # y-coordinate of a given contaminant inside the imagette (window)
        offx_c = x_c_im - pxc[psf_idx] # x-coordinate of the offset between the center of each contaminant and each PSF
        offy_c = y_c_im - pyc[psf_idx] # y-coordinate of the offset between the center of each contaminant and each PSF

        # Now we have to build an 'imagette' for every contaminant
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

        Ic_acc = np.sum(Ic, axis=0)                                                  # numpy array containing the contribution from all contaminants
        f_tot = It + Ic_acc                                                          # total flux of an imagette (TARGET + all contaminants)
        NSR1h, w_t = aperture(ft=It, fc=Ic_acc, sb=sb, sd=sd, sq=sq)                 # nominal mask for a given TARGET
        w_t_key = mask_to_bitmask(w_t)                                               # mask key for the nominal mask
        w_t_size = np.count_nonzero(w_t)                                             # mask size for the nominal mask
        # sprk, sprk_max, SPR_tot, n_bad = SPR(SPR_crit=SPR_crit, n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)
        sprk, sprk_max, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_t) # sprk for every contaminant (as well as sprkmax and spr-tot)
        # Now we choose randomly a value for dback and for td
        # dback_index = np.random.randint(0, len(del_back))
        # td_index = np.random.randint(0, len(tr_dur))
        # We obtain the values of dback and td from the catalogue
        # dback = del_back[dback_index]
        # td = tr_dur[td_index]
        # We compute the critical SPR now
        SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr) # spr_crit w.r.t. the nominal mask
        n_bad = np.sum(sprk > SPR_crit) # N_bad
        ind_sprk = np.argmax(sprk)      # index of the contaminant with the highest sprk w.r.t. nominal mask
        m_c_bad = m_c[ind_sprk]         # magnitude for the contaminant with the highest sprk value w.r.t the nominal mask
        # Now we select the (intensity per pixel array) 'imagette' of the contaminant star with the highest sprk value
        # w.r.t. the nominal mask
        Ic_max = Ic[ind_sprk]

        # Now we find the distance between the target and the contaminant star with the highest value of sprk w.r.t the
        # nominal mask
        dist_bad = (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2
        # Now we define a term that contains the flux of the targets and their respective contaminants except the
        # one with the highest value of SPRk (i.e. the contaminant of interest)
        Itc_acc = It + Ic_acc - Ic_max

        NSR1h_c, w_c = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq) # secondary mask for the most problematic contaminant
        w_c_key = mask_to_bitmask(w_c)                                      # mask key for the secondary mask
        w_c_size = np.count_nonzero(w_c)                                    # mask size for the secondary mask
        sprk_sec, sprk_max_sec, SPR_tot_sec = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_c) # sprk for every contaminant (as well as sprkmax and spr-tot) w.r.t. secondary mask

        # We compute the flux over the secondary mask and the extended secondary mask
        f_beb = Ic_acc * w_c
        f_t_c = It * w_c
        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))

        # We compute spr_tot_c, that is, the expression given in Marchiori presentation for PLATO week #8
        spr_tot_c = np.sum(Itc_acc * w_c) / f_tot_c

        # We compute now the delta_obs and etas for the two apertures
        delta_obs_t = sprk[ind_sprk] * dback                                         # observed transit depth in the nominal mask
        delta_obs_c = (1 - spr_tot_c) * dback                                        # observed transit depth in the secondary mask
        eta_t = sprk[ind_sprk] * np.sqrt(td * ntr) * dback / (NSR1h * (1 - SPR_tot)) # signal statistical significance in the nominal mask
        eta_c = np.sqrt(td * ntr) * dback / NSR1h_c                                  # signal statistical significance in the secondary mask

        # -------------------------------------------NOMINAL COB-------------------------------------------------------#
        eta_cob, sigma_1_24, abs_cob, eta_cob_wrong, sigma_1_24_wrong = centroid_shift(w=w_t, Ik=Ic_max, I_t=It,
        I_contaminants=Ic_acc, sprk=sprk[ind_sprk], dback=dback, sb=sb, sd=sd, sq=sq, td=td, ntr=ntr)
        # -------------------------------------------NOMINAL COB-------------------------------------------------------#

        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        eta_cob_c, sigma_1_24_c, abs_cob_c, eta_cob_c_wrong, sigma_1_24_c_wrong = centroid_shift(w=w_c, Ik=Ic_max,
        I_t=It, I_contaminants=Ic_acc, sprk=sprk_sec[ind_sprk], dback=dback, sb=sb, sd=sd, sq=sq, td=td, ntr=ntr)
        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        ################################################################################################################
        #                                               EXTENDED MASK                                                  #
        ################################################################################################################
        w_ext = extended_binary_mask(w_t, W=1) # extended mask (1 pixel ring around the nominal mask)
        w_ext_key = mask_to_bitmask(w_ext)     # mask key for the extended mask
        w_ext_size = np.count_nonzero(w_ext)   # mask size for the extended mask

        # Now we compute all the metrics associated with every extended mask.
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext)) / np.sum(It * w_ext)
        NSR_ext_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext
        sprk_ext, sprk_max_ext, SPR_tot_ext = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_ext)
        SPR_crit_ext = spr_crit(dback=dback, SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h, td=td, ntr=ntr)
        n_bad_ext = np.sum(sprk_ext > SPR_crit_ext)                           # N_bad for the extended mask
        delta_obs_ext = sprk_ext[ind_sprk] * dback                            # w.r.t. the nominal mask
        eta_ext = sprk_ext[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h # w.r.t the nominal mask 


        # ------------------------------------------EXTENDED COB-------------------------------------------------------#
        eta_cob_ext, sigma_1_24_ext, abs_cob_ext, eta_cob_ext_wrong, sigma_1_24_ext_wrong = centroid_shift(w=w_ext, Ik=Ic_max, 
        I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[ind_sprk], dback=dback, sb=sb, sd=sd, sq=sq, td=td, ntr=ntr)
        # ------------------------------------------EXTENDED COB-------------------------------------------------------#
        ################################################################################################################
        #                                               EXTENDED MASK                                                  #
        ################################################################################################################

        ################################################################################################################
        #                                   TESTING  J.C. Bray et al.'s ASSUMPTION OF A 2 x 2 MASK                     #
        ################################################################################################################

        # The mask has to contain the 4 centered pixels
        w_bray = np.zeros((6, 6))
        w_bray[2:4, 2:4] = 1

        # Now we compute the NSR w.r.t bray's mask
        NSR_bray = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_bray)) / np.sum(It * w_bray)
        NSR_bray_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_bray
        sprk_bray, sprk_max_bray, SPR_tot_bray = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_bray)
        SPR_crit_bray = spr_crit(dback=dback, SPR_tot=SPR_tot_bray, nsr=NSR_bray_1h, td=td, ntr=ntr)
        n_bad_bray = np.sum(sprk_bray > SPR_crit_bray) # N_bad for Bray's mask

        ################################################################################################################
        #                                          END OF TESTING Bray et al.'s ASSUMPTION                             #
        ################################################################################################################

        # The number of false positives given by the extended mask such that eta_ext > eta_t is given by
        # n_eff_ext = len(np.where((sprk_ext > SPR_crit_ext) & (sprk[ind_sprk] > SPR_crit) & (sprk_ext > sprk[ind_sprk])
        # )[0])
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              '++++++')
        print('Target ID = ', ID_target[k])
        print('Total number of contaminant stars for this target = ', n_c)
        print('Nbad for this target = ', n_bad)
        print('Magnitude difference between the target and the contaminant with the highest sprk = ',
              m_c[ind_sprk] - m_t)
        print('Distance between the target and the contaminant with the highest sprk = ',
              (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        print('Transit depth of the eclipsing binary with the highest sprk = ', dback, 'ppm')
        print('Transit duration of the eclipsing binary with the highest sprk = ', td, 'hours')
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              '++++++\n')

        # Now we save the important metrics for every target w.r.t the nominal mask
        save_info = np.array([ID_t, m_t, n_c, m_c_bad, dist_bad, w_t_key, w_t_size, NSR1h, n_bad, SPR_crit, sprk[ind_sprk], SPR_tot, eta_t, delta_obs_t, abs_cob, eta_cob, 
                              sigma_1_24, eta_cob_wrong, sigma_1_24_wrong])
        
        # Now we save the important metrics w.r.t the secondary mask
        save_info_contaminant = np.array([ID_t, m_t, w_c_key, w_c_size, NSR1h_c, spr_tot_c, eta_c, delta_obs_c, abs_cob_c, eta_cob_c, sigma_1_24_c, eta_cob_c_wrong, sigma_1_24_c_wrong])

        # Now we save the important metrics w.r.t the extended mask
        save_info_ext = np.array([ID_t, m_t, w_ext_key, w_ext_size, NSR_ext_1h, sprk_ext[ind_sprk], SPR_crit_ext, eta_ext, delta_obs_ext, abs_cob_ext, eta_cob_ext, 
                                  sigma_1_24_ext, n_bad_ext, eta_cob_ext_wrong, sigma_1_24_ext_wrong])

        # Now we save the import metrics w.r.t the 4 pixel mask described by Bray et al. 2023
        save_info_bray = np.array([ID_t, m_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray, sprk_bray[ind_sprk], SPR_tot_bray])
        
        return save_info, save_info_contaminant, save_info_ext, save_info_bray     

    for k in range(n_t):
        result = process_target(k)
        save_info[counter] = result[0]
        save_info_contaminant[counter] = result[1]
        save_info_ext[counter] = result[2]
        save_info_bray[counter] = result[3]
        counter = counter + 1
        

    print("NUMBER OF PROCESSED TARGETS: %i" % counter)

save_info = save_info[0:counter]
save_info_contaminant = save_info_contaminant[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]

# Now it is time to save the metrics into several .npy files, each file for each mask
np.save(DIRout + 'targets_P5.npy', save_info)
np.save(DIRout + 'targets_P5_contaminant.npy', save_info_contaminant)
np.save(DIRout + 'targets_P5_extended.npy', save_info_ext)
np.save(DIRout + 'targets_P5_bray.npy', save_info_bray)

fout = open(DIRout + 'star_count.txt', 'w')
for i in range(nP):
    Pi = Pmin + i * binsize
    fout.write('%.2f %i\n' % (Pi, n_star_p_bin[i]))
fout.close()