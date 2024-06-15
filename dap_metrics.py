# ------------------------------------------------
# First we import the main libraries and modules
import numpy as np # type: ignore
import spline2dbase # type: ignore
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant
from imagette import barycenter, window, ran_unique_int, centroid_shift
from NSR import spr_crit, aperture_computation, SPR, mask_to_bitmask, extended_binary_mask
from pylab import * # type: ignore
import sys
import os
# ------------------------------------------------

# ------------------------------------------------
# CONFIGURATION PARAMETERS

# Parameters relative to all the relevant paths
##Rezadata = np.load('/home/fercho/double-aperture-photometry/test_results_multiprocessing/' + 'targets_P5.npy')
#RezaID = int(Rezadata[:, 0])
#print(RezaID)
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues 
PSFfile = 'PSF.npz'                                                   # processed PSF files
DIRout = 'test_results/' # storage directory
#output_run_dir = DIRout + sys.argv[1] + "/"
#os.makedirs(output_run_dir, exist_ok = True)

# Parameters for the imagette and PSF decomposition
size_im_x = 6  # size of the imagette (x-direction)
size_im_y = 6  # size of the imagette (y-direction)
subres = 128   # resolution of the PSF
bsres = 20     # resolution of the b-spline decomposition of the PSF

# Parameters for the NSR
sb = (45 * 21) # Background noise from zodiacal light in e-/px(poisson noise)times integration time (21 sec.)
sd = 50.2      # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2       # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
#dback = 85000  # transit depth in ppm
#td = 4         # transit duration in hours
ntr = 3        # number of transits in one hour

# Parameters for the magnitude intervals
n_tar = 10                            # number of targets per magnitude interval
Pmin = 10                               # minimum magnitude
Pmax = 13                               # maximum magnitude
binsize = 0.5                           # binsize around every magnitude value
nP = int((Pmax - Pmin) / binsize + 1)   # number of bins
# ------------------------------------------------

# ------------------------------------------------
# MORE CONFIGURATION PARAMETERS

data = np.load(cataDIR + 'SFP_DR3_20230101.npy') # star catalogue from GAIA
psfdata = np.load(PSFfile)                       # processed PSFs 
del_back, tr_dur = np.loadtxt(cataDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt', unpack=True, usecols=[0, 1]) # transit depth and duration from Kepler Eclipsing Binary Catalogue

ID = np.arange(0, data.shape[0]) # ID for every star in the catalogue

# down select the full catalog to a single target
#RI = 0
#data = data[ID==RezaID[RI],:]
#ID = ID[RezaID[RI]]

x_star = data[:, 3]              # x-coordinate on the focal plane for every star in the catalogue
y_star = data[:, 4]              # y-coordinate on the focal plane for every star in the catalogue
psfbs = psfdata['psfbs']         # PSF b-spline decomposition
pxc = psfdata['pxc']             # x-coordinate of the PSF b-spline decomposition centroid
pyc = psfdata['pyc']             # y-cooridnate of the PSF b-spline decomposition centroid
xpsf_pix = psfdata['xpsf_pix']   # x-coordinate of the PSF in pixel
ypsf_pix = psfdata['ypsf_pix']   # y-coordinate of the PSF in pixel

# Now we use a seed for obtaining the same number of targets and dback and td values
np.random.seed(300)

file_out = open(DIRout + 'metrics_fer.txt', 'w')
save_info = np.zeros((n_tar * nP, 106))     # numpy arr. to store the metrics for the nominal mask
save_info_sec = np.zeros((n_tar * nP, 16)) # numpy arr. to store the metrics for the secondary mask
save_info_ext = np.zeros((n_tar * nP, 105)) # numpy arr. to store the metrics for the extended mask
save_info_bray = np.zeros((n_tar * nP, 8)) # numpy arr. to store the metrics for Bray's 2 x 2 mask
n_star_p_bin = np.zeros(nP)                # numpy arr. to store the number of stars per bin

# ------------------------------------------------

# ------------------------------------------------
# BEGIN OF THE METRICS COMPUTATION

# We define a counter before the main for loop
counter = 0

# We start the main for loop. This goes for every magnitude bin
for i in range(nP):
    Pi = Pmin + i * binsize
    # Now we define a mask to get only the targets whithin the i-th magnitude bin
    mask = (data[:, 2] >= Pi - binsize / 2.) & (data[:, 2] <= Pi + binsize / 2.)
    n_star_p_bin[i] = mask.sum() # number of targets in the current bin
    targets_P5 = data[mask, :]   # magnitude of every target in the current bin
    ID_target = ID[mask]         # ID of every target in the current bin
    # Now we choose (randomly) only 'n_tar' targets from every bin. These are our new "TARGETS".
    j = ran_unique_int(n=n_tar, interval=[0, targets_P5.shape[0] - 1])
    print(j)
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
    
    # We start the principal processing routine. This goes for every target.
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
                Ic[o - 1, :, :] = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y, offx=offx_c[o - 1], offy=offy_c[o - 1])
                COBx_c, COBy_c = barycenter(Ic[o - 1], subres=1)
                # Let's obtain the value of the reference flux for every contaminant star
                f_ref_c = reference_flux_contaminant(f_ref_t, m_c[o - 1], m_t)
                # Let's calculate the Intensity per pixel of the imagette of every contaminant star
                Ic[o - 1, :, :] = f_ref_c * Ic[o - 1]

        Ic_acc = np.sum(Ic, axis=0)                                                  # numpy array containing the contribution from all contaminants
        f_tot = It + Ic_acc                                                          # total flux of an imagette (TARGET + all contaminants)
        nominal_mask = aperture_computation(ft=It, fc=Ic_acc, sb=sb, sd=sd, sq=sq)
        nominal_mask_key = mask_to_bitmask(nominal_mask)
        nominal_mask_size = np.count_nonzero(nominal_mask)
        nsr_nominal_mask = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * (nominal_mask**2))) / np.sum(It * nominal_mask) # E1. 11 in Marchiori paper
        nsr_1h_24_cameras = ((10 ** 6) / (12 * np.sqrt(24))) * nsr_nominal_mask
        nsr_1h_6_cameras = ((10 ** 6) / (12 * np.sqrt(6))) * nsr_nominal_mask

        # sprk, sprk_max, SPR_tot, n_bad = SPR(SPR_crit=SPR_crit, n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)
        sprk, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=nominal_mask) # sprk for every contaminant (as well as sprkmax and spr-tot)
        sprk_6_cameras,  SPR_tot_6_cameras = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=nominal_mask)
        # 10 values of transit depth
        dback = np.zeros(n_c)
        td = np.zeros(n_c)
        # Now we choose randomly a single value for dback and for td
        #dback_index = np.random.randint(0, len(del_back))
        #td_index = np.random.randint(0, len(tr_dur))
        # We obtain then single values of dback and td from the catalogue
        #dback = del_back[dback_index]
        #td = tr_dur[td_index]
        #print(dback)
        # We obtain 10 random values for transit depth and transit duration
        for n in range(n_c):
            dback_index = np.random.randint(0, len(del_back))
            td_index = np.random.randint(0, len(tr_dur))
            
            dback[n] = del_back[dback_index]
            td[n] = tr_dur[td_index]
        
        # We compute the critical SPR now
        #SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr) # spr_crit w.r.t. the nominal mask
        SPR_crit_24_cameras = spr_crit(dback=dback[0], SPR_tot=SPR_tot, nsr=nsr_1h_24_cameras, td=td[0], ntr=ntr) # spr_crit w.r.t. the nominal mask
        SPR_crit_6_cameras = spr_crit(dback=dback[0], SPR_tot=SPR_tot, nsr=nsr_1h_6_cameras, td=td[0], ntr=ntr)
        
        n_bad = np.sum(sprk > SPR_crit_24_cameras) # N_bad
        n_bad_24_cameras = np.sum(sprk > SPR_crit_24_cameras)
        n_bad_6_cameras = np.sum(sprk > SPR_crit_6_cameras)
        ind_sprk = np.argmax(sprk)      # index of the contaminant with the highest sprk w.r.t. nominal mask
        m_c_bad = m_c[ind_sprk]         # magnitude for the contaminant with the highest sprk value w.r.t the nominal mask
        # Now we select the (intensity per pixel array) 'imagette' of the contaminant star with the highest sprk value
        # w.r.t. the nominal mask
        Ic_max = Ic[ind_sprk]
        dback_contaminant_highest_spr = dback[ind_sprk]
        td_contaminant_highest_spr = td[ind_sprk]

        # Now we find the distance between the target and the contaminant star with the highest value of sprk w.r.t the
        # nominal mask
        dist_bad = (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2
        # Now we define a term that contains the flux of the targets and their respective contaminants except the
        # one with the highest value of SPRk (i.e. the contaminant of interest)
        Itc_acc = It + Ic_acc - Ic_max

        #NSR1h_c, NSR1h_c_6, w_c, w_c_6 = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq) # secondary mask for the most problematic contaminant
        secondary_mask = aperture_computation(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq)
        secondary_mask_key = mask_to_bitmask(secondary_mask)
        secondary_mask_size = np.count_nonzero(secondary_mask)
       
        nsr_secondary_mask = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * (secondary_mask**2))) / np.sum(It * secondary_mask) # E1. 11 in Marchiori paper
        nsr_1h_24_cameras_secondary_mask = ((10 ** 6) / (12 * np.sqrt(24))) * nsr_secondary_mask 
        nsr_1h_6_cameras_secondary_mask = ((10 ** 6) / (12 * np.sqrt(6))) * nsr_secondary_mask 

        sprk_sec, SPR_tot_sec = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=secondary_mask) # sprk for every contaminant (as well as sprkmax and spr-tot) w.r.t. secondary mask
        # We compute the flux over the secondary mask and the extended secondary mask
        f_beb = Ic_acc * secondary_mask
        f_t_c = It * secondary_mask
        #f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))
        f_tot_c = np.sum(f_t_c + f_beb)
        #if np.sum(f_t_c + f_beb) == f_tot_c:
        #    print('A goevo, putos')
        #else:
        #    breakpoint

        # We compute spr_tot_c, that is, the expression given in Marchiori presentation for PLATO week #8
        spr_tot_secondary_mask = np.sum(Itc_acc * secondary_mask) / f_tot_c

        # We compute now the delta_obs and etas for the two apertures
        delta_obs_nominal_mask = sprk[ind_sprk] * dback_contaminant_highest_spr                                         # observed transit depth in the nominal mask
        delta_obs_secondary_mask = (1 - spr_tot_secondary_mask) * dback_contaminant_highest_spr                                        # observed transit depth in the secondary mask
        eta_t = sprk[ind_sprk] * np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / (nsr_1h_24_cameras * (1 - SPR_tot)) # signal statistical significance in the nominal mask
        eta_t_6_cameras = sprk[ind_sprk] * np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr/ (nsr_1h_6_cameras * (1 - SPR_tot)) # signal statistical significance in the nominal mask
        eta_c = np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / nsr_1h_24_cameras_secondary_mask # signal statistical significance in the secondary mask
        eta_c_6_cameras = np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / nsr_1h_6_cameras_secondary_mask
    
        
        eta_true_positive_24_cameras_earth_like = 84 * np.sqrt(13 * 3) / nsr_1h_24_cameras
        eta_true_positive_6_cameras_jovian_planet =  10100 * np.sqrt(29.6 * 3) / nsr_1h_6_cameras
        eta_true_positive_24_cameras_super_earth = 522 * np.sqrt(42 * 3) / nsr_1h_24_cameras
        eta_true_positive_6_cameras_super_earth = 522 * np.sqrt(42 * 3) / nsr_1h_6_cameras
        
        nsprmax = min(10, n_c)
        SPRK_10first = np.zeros(10)
        eta_10first = np.zeros(10)
        eta_cob_10first = np.zeros(10)
        sigma_cob_10first = np.zeros(10)
        abs_cob_shift_10first = np.zeros(10)
        eta_cob_10first_6_cameras = np.zeros(10)
        sigma_cob_10first_6_cameras = np.zeros(10)
        abs_cob_shift_10first_6_cameras = np.zeros(10)
        dback_10first = np.zeros(10)
        td_10first = np.zeros(10)
        #dback_10first_index = np.random.randint(0, len(del_back), size=10)
        # sorting the SPRk by decreasing order and taking the 10 first values
        sprk_sorted_index = (np.argsort(sprk)[::-1])
        SPRK_10first[0:nsprmax] = sprk[sprk_sorted_index[0:nsprmax]]
        dback_10first[0:nsprmax] = dback[sprk_sorted_index[0:nsprmax]]
        td_10first[0:nsprmax] = td[sprk_sorted_index[0:nsprmax]]
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            eta_10first[l] = sprk[m] * np.sqrt(td_10first[l] * ntr) * dback_10first[l] / nsr_1h_24_cameras / (1. - SPR_tot)
            eta_cob_10first[l], sigma_cob_10first[l], abs_cob_shift_10first[l] = centroid_shift(w=nominal_mask, Ik=Ic[m],
            n_cam=24, I_t=It, I_contaminants=Ic_acc, sprk=sprk[m], dback=dback_10first[l], sb=sb, sd=sd, sq=sq, td=td_10first[l], ntr=ntr)
            eta_cob_10first_6_cameras[l], sigma_cob_10first_6_cameras[l], abs_cob_shift_10first_6_cameras[l] = centroid_shift(w=nominal_mask, Ik=Ic[m],
            n_cam=6, I_t=It, I_contaminants=Ic_acc, sprk=sprk[m], dback=dback_10first[l], sb=sb, sd=sd, sq=sq, td=td_10first[l], ntr=ntr)

        # -------------------------------------------NOMINAL COB-------------------------------------------------------#
        ## 24 cameras
        eta_cob, sigma_1_24, abs_cob = centroid_shift(w=nominal_mask, Ik=Ic_max, n_cam=24, I_t=It, I_contaminants=Ic_acc, 
                                            sprk=sprk[ind_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
        
        ## 6 cameras
        eta_cob_6_cameras, sigma_1_6_cameras, abs_cob_6_cameras = centroid_shift(w=nominal_mask, Ik=Ic_max, n_cam=6, I_t=It,
        I_contaminants=Ic_acc, sprk=sprk[ind_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr) 
        # -------------------------------------------NOMINAL COB-------------------------------------------------------#

        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        ## 24 cameras
        eta_cob_c, sigma_1_24_c, abs_cob_c = centroid_shift(w=secondary_mask, Ik=Ic_max, n_cam=24, I_t=It, I_contaminants=Ic_acc, 
                                            sprk=sprk_sec[ind_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
        
        ## 6 cameras
        eta_cob_c_6_cameras, sigma_1_6_cameras_c, abs_cob_c_6_cameras = centroid_shift(w=secondary_mask, Ik=Ic_max, n_cam=6, I_t=It, I_contaminants=Ic_acc, 
                                            sprk=sprk_sec[ind_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
        # ------------------------------------------SECONDARY COB------------------------------------------------------#
        
        ################################################################################################################
        #                                               EXTENDED MASK                                                  #
        ################################################################################################################
        extended_mask = extended_binary_mask(nominal_mask, W=1) # extended mask (1 pixel ring around the nominal mask)
        extended_mask_key = mask_to_bitmask(extended_mask)     # mask key for the extended mask
        extended_mask_size = np.count_nonzero(extended_mask)   # mask size for the extended mask

        w_ext_6_cameras = extended_binary_mask(nominal_mask, W=1)
        
        # Now we compute all the metrics associated with every extended mask.
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * extended_mask)) / np.sum(It * extended_mask) # E1. 11 in Marchiori paper
        NSR_ext_1h_24_cameras = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext
        NSR_ext_1h_6_cameras = ((10 ** 6) / (12 * np.sqrt(6))) * NSR_ext
        sprk_ext, SPR_tot_ext = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=extended_mask)
        SPR_crit_ext = spr_crit(dback=dback_10first[0], SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h_24_cameras, td=td_10first[0], ntr=ntr)
        SPR_crit_ext_6_cameras = spr_crit(dback=dback_10first[0], SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h_6_cameras, td=td_10first[0], ntr=ntr)
        n_bad_ext = np.sum(sprk_ext > SPR_crit_ext)                           # N_bad for the extended mask
        delta_obs_ext = sprk_ext[ind_sprk] * dback_10first[0]                            # w.r.t. the nominal mask
        eta_ext = sprk_ext[ind_sprk] * np.sqrt(td_10first[0] * ntr) * dback_10first[0] / (NSR_ext_1h_24_cameras *(1 - SPR_tot_ext)) # w.r.t the nominal mask
        #eta_ext_6_cameras = sprk_ext[ind_sprk] * np.sqrt(td * ntr) * dback / (NSR_ext_1h_6_cameras *(1 - SPR_tot_ext))

        # ------------------------------------------EXTENDED COB-------------------------------------------------------#
        ## 24 cameras w.r.t. the most significant contaminant
        eta_cob_ext, sigma_1_24_ext, abs_cob_ext = centroid_shift(w=extended_mask, Ik=Ic_max, n_cam=24, 
        I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[ind_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
        # ------------------------------------------EXTENDED COB-------------------------------------------------------#
        
        # Now we will take the first 10 contaminants in terms of SPR
        
        # We take the 10 first contaminants in terms of SPRk to use them for the extended mask
        nsprmax = min(10, n_c)
        sprk_sorted_index = (np.argsort(sprk)[::-1])
        
        SPRK_ext_10first = np.zeros(10)
        eta_ext_10first = np.zeros(10)
        eta_ext_10first_6_cameras = np.zeros(10)
        eta_cob_ext_10first = np.zeros(10)
        eta_cob_ext_10first_6_cameras = np.zeros(10)
        sigma_cob_ext_10first = np.zeros(10)
        sigma_cob_ext_10first_6_cameras = np.zeros(10)
        abs_cob_shift_ext_10first = np.zeros(10)
        abs_cob_shift_ext_10first_6_cameras = np.zeros(10)
        dback_ext_10first = np.zeros(10)
        td_ext_10first = np.zeros(10)

        
        SPRK_ext_10first[0:nsprmax] = sprk_ext[sprk_sorted_index[0:nsprmax]]
        dback_ext_10first[0:nsprmax] = dback[sprk_sorted_index[0:nsprmax]]
        td_ext_10first[0:nsprmax] = td[sprk_sorted_index[0:nsprmax]]
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            eta_ext_10first[l] = sprk_ext[m] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_1h_24_cameras * (1 - SPR_tot_ext))       
            eta_ext_10first_6_cameras[l] = sprk_ext[l] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_1h_6_cameras * (1 - SPR_tot_ext))
            eta_cob_ext_10first[l], sigma_cob_ext_10first[l], abs_cob_shift_ext_10first[l] = centroid_shift(w=extended_mask, Ik=Ic[m],  
            n_cam=24, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)
            eta_cob_ext_10first_6_cameras[l], sigma_cob_ext_10first_6_cameras[l], abs_cob_shift_ext_10first_6_cameras[l] = centroid_shift(w=extended_mask, Ik=Ic[m],
            n_cam=6, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)
       
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
        sprk_bray,  SPR_tot_bray = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_bray)
        SPR_crit_bray = spr_crit(dback=dback_10first[0], SPR_tot=SPR_tot_bray, nsr=NSR_bray_1h, td=td_10first[0], ntr=ntr)
        n_bad_bray = np.sum(sprk_bray > SPR_crit_bray) # N_bad for Bray's mask
        ################################################################################################################
        #                                          END OF TESTING Bray et al.'s ASSUMPTION                             #
        ################################################################################################################
        # The number of false positives given by the extended mask such that eta_ext > eta_t is given by
        # n_eff_ext = len(np.where((sprk_ext > SPR_crit_ext) & (sprk[ind_sprk] > SPR_crit) & (sprk_ext > sprk[ind_sprk])
        # )[0])
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('Target ID = ', ID_target[k])
        print('Total number of contaminant stars for this target = ', n_c)
        print('Nbad for this target = ', n_bad)
        print('Nbad for this target (24 cameras) =', n_bad_24_cameras)
        print('Nbad for this target (6 cameras) =', n_bad_6_cameras)
        print('Transit depth of the contaminant with the highest spr =', dback_contaminant_highest_spr, 'ppm')
        print('Magnitude difference between the target and the contaminant with the highest sprk = ',
              m_c[ind_sprk] - m_t)
        print('Distance between the target and the contaminant with the highest sprk = ',
              (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        print('Secondary mask size =', secondary_mask_size)
        print('sprk_tot_24_cameras:', SPR_tot)
        print('sprk_tot_6_cameras', SPR_tot_6_cameras)
        print('Spr_crit =', SPR_crit_24_cameras)
        print('NSR_T is:', nsr_1h_24_cameras)
        print('min(nsr_agg_1h_24_cameras)', nsr_1h_24_cameras)
        print('nsr', nsr_1h_24_cameras)
        print('nsr_6_cameras', nsr_1h_6_cameras)
        print('NSR_ext', NSR_ext)
        print('NSR_c is:',  nsr_1h_24_cameras_secondary_mask)
        print('eta_t is:', eta_t)
        print('eta_c_6_cameras is:', eta_c_6_cameras)
        print('NSR1H_6_C',  nsr_1h_6_cameras_secondary_mask)
        print('spr_tot_c', spr_tot_secondary_mask)
        print('spr_tot_c_6_cameras', SPR_tot_sec)
        print('spr_t is:', sprk[ind_sprk])
        print('SPRk_10_first_ext are:', SPRK_ext_10first)
        print('eta_10first_ext are:', eta_ext_10first)
        print('eta_10first are:', eta_10first)
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

        # Now we save the important metrics for every target w.r.t the nominal mask
        save_info = np.array([ID_t, m_t, n_c, m_c_bad, dist_bad, nominal_mask_key, nominal_mask_size, nsr_1h_24_cameras, n_bad, SPR_crit_24_cameras, sprk[ind_sprk], SPR_tot, eta_t, delta_obs_nominal_mask, abs_cob, eta_cob, 
                              sigma_1_24])
        
        save_info = np.append(save_info, SPRK_10first)
        save_info = np.append(save_info, eta_10first)
        save_info = np.append(save_info, eta_true_positive_24_cameras_earth_like)
        save_info = np.append(save_info, eta_true_positive_6_cameras_jovian_planet)
        save_info = np.append(save_info, eta_true_positive_24_cameras_super_earth)
        save_info = np.append(save_info, eta_true_positive_6_cameras_super_earth)
        save_info = np.append(save_info, eta_cob_6_cameras)
        save_info = np.append(save_info, abs_cob_6_cameras)
        save_info = np.append(save_info, sigma_1_6_cameras)
        save_info = np.append(save_info, SPR_crit_6_cameras)
        save_info = np.append(save_info, eta_t_6_cameras)
        save_info = np.append(save_info, eta_cob_10first)
        save_info = np.append(save_info, sigma_cob_10first)
        save_info = np.append(save_info, abs_cob_shift_10first)
        save_info = np.append(save_info, eta_cob_10first_6_cameras)
        save_info = np.append(save_info, sigma_cob_10first_6_cameras)
        save_info = np.append(save_info, abs_cob_shift_10first_6_cameras)

        # Now we save the important metrics w.r.t the secondary mask
        save_info_contaminant = np.array([ID_t, m_t, secondary_mask_key, secondary_mask_size, nsr_1h_24_cameras_secondary_mask, spr_tot_secondary_mask, eta_c, delta_obs_secondary_mask, abs_cob_c, eta_cob_c, sigma_1_24_c])
        save_info_contaminant = np.append(save_info_contaminant, eta_c_6_cameras)
        save_info_contaminant = np.append(save_info_contaminant, eta_cob_c_6_cameras)
        save_info_contaminant = np.append(save_info_contaminant, abs_cob_c_6_cameras)
        save_info_contaminant = np.append(save_info_contaminant, sigma_1_6_cameras_c)
        save_info_contaminant = np.append(save_info_contaminant, SPR_tot_sec)

        # Now we save the important metrics w.r.t the extended mask
        save_info_ext = np.array([ID_t, m_t, extended_mask_key, extended_mask_size, NSR_ext_1h_24_cameras, sprk_ext[ind_sprk], SPR_crit_ext, eta_ext, delta_obs_ext, abs_cob_ext, eta_cob_ext, 
                                  sigma_1_24_ext, n_bad_ext, SPR_tot_ext])
        
        save_info_ext = np.append(save_info_ext, SPRK_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_ext_10first_6_cameras)
        save_info_ext = np.append(save_info_ext, NSR_ext_1h_6_cameras)
        save_info_ext = np.append(save_info_ext, eta_cob_ext_10first)
        save_info_ext = np.append(save_info_ext, sigma_cob_ext_10first)
        save_info_ext = np.append(save_info_ext, abs_cob_shift_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_cob_ext_10first_6_cameras)
        save_info_ext = np.append(save_info_ext, sigma_cob_ext_10first_6_cameras)
        save_info_ext = np.append(save_info_ext, abs_cob_shift_ext_10first_6_cameras)


        # Now we save the import metrics w.r.t the 4 pixel mask described by Bray et al. 2023
        save_info_bray = np.array([ID_t, m_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray, sprk_bray[ind_sprk], SPR_tot_bray])
        
        file_out.write('%.2f %i %.2f\n'% (ID_t, n_bad, sprk[ind_sprk]))
        
        return save_info, save_info_contaminant, save_info_ext, save_info_bray     

    for k in range(len(ID_target)):
        result = process_target(k)
        save_info[counter] = result[0]
        save_info_sec[counter] = result[1]
        save_info_ext[counter] = result[2]
        save_info_bray[counter] = result[3]
        counter = counter + 1
        
       
        

    print("NUMBER OF PROCESSED TARGETS: %i" % counter)

save_info = save_info[0:counter]
save_info_sec = save_info_sec[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]

# Now it is time to save the metrics into several .npy files, each file for each mask
np.save(DIRout + 'targets_P5.npy', save_info)
np.save(DIRout + 'targets_P5_secondary.npy', save_info_sec)
np.save(DIRout + 'targets_P5_extended.npy', save_info_ext)
np.save(DIRout + 'targets_P5_bray.npy', save_info_bray)
file_out.close()

fout = open(DIRout + 'star_count.txt', 'w')
for i in range(nP):
    Pi = Pmin + i * binsize
    fout.write('%.2f %i\n' % (Pi, n_star_p_bin[i]))
fout.close()
