import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import spline2dbase
import scipy.signal
import h5py as h5py
import numpy as np
import sys
from fitting_psf import from_mm_2_pix, from_pix_2_mm, closest_psf, contaminants, reference_flux_target, \
    reference_flux_contaminant
from imagette import catalogue, list_psf, barycenter, gauss, window, ran_unique_int, ploting_imagettes, ploting_nsr, \
    ploting_nsr_s, \
    ploting_initial,psf_gauss_int
from NSR import spr_crit, aperture, NSRn, nsr_AGG, SPR, mask_to_bitmask, bitmask_to_mask, extended_binary_mask
from pylab import *
from tqdm import tqdm # only used for convenience, to monitor progress in our  a calculation
import multiprocessing
import math
import time
import traceback

#--------------------------------------------------
# CONFIGURATION PARAMETERS

verbose = True
CatalogueDIR = '/home/reza/plato/share/catalogues/'
# CatalogueDIR = '/volumes/astro/sismo/general/plato/web/grids/catalogues/'
DIRout = 'test/'
PSFFileName = 'PSF.npz'
# CatalogueFileName = 'SFP_DR3_20220831.npy'
CatalogueFileName = 'SFP_DR3_20230101.npy'

gauss_psf = False
gauss_width_x = 0.5
gauss_width_y = 0.5

bsres = 20 # resolution adopted for the b-spline decomposition

mp =  False # multiprocessing mode
nproc_max = 4 # number of processors

Pmin = 8
Pmax = 13
binsize = 0.5
nStar = 300
Delta_P_max = 15.
distance_max = 7
n_c_max = 300

extend_2ndmask = 0

# Parameters for the NSR
sb = (45. * 21)  # Background noise form zodiacal light in units of e-/px after multiplying by the integration time (poisson noise)
sd = 50.2  # Overall detector noise (including readout at beginning of life, smearing and dark current) in units of e- rms/px
sq = 7.2  # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries
dback = 85000  # transit depth in ppm

td = 4  # transit duration in hours
ntr = 3  # number of transits

seed = 300

#--------------------------------------------------
data = catalogue(CatalogueDIR+CatalogueFileName)

psfdata = np.load(PSFFileName)

psfbs = psfdata['psfbs']
pxc = psfdata['pxc']
pyc = psfdata['pyc']
xpsf_pix = psfdata['xpsf_pix']
ypsf_pix = psfdata['ypsf_pix']

# Define an ID for every target
ID = np.arange(0, data.shape[0])

# Now we save the x and y coordinates on the focal plane of all the stars in the catalogue
x_star = data[:, 3]
y_star = data[:, 4]

nP = int(round ( (Pmax-Pmin)/binsize + 1 ))

# Define a numpy array for saving the metrics of interest (Target ID, magnitude, N_bad, etc.) (Is hard-coded now)
save_info = np.zeros((nStar * nP, 52))

# for the secondary mask
save_info_2ndmask =  np.zeros((nStar * nP,64))

# The same for the extended mask
save_info_ext = np.zeros((nStar * nP,51))

# The same for bray's et al. assu#mption of using 2 x 2 masks
save_info_bray = np.zeros((nStar * nP, 6))

# Now we choose the random targets using Réza's function
np.random.seed(seed)
n_star_p_bin = np.zeros((nP))

def cob_shift(Itot,Ic,w,dback):
    # compute the COB shift and the associated uncertainty

    db = dback*1e-6
    x = (np.arange(0, Ic.shape[1]) + 0.5)
    y = (np.arange(0, Ic.shape[0]) + 0.5)
    x, y = np.meshgrid(x, y)

    Icw= Ic*w
    Itotw = Itot*w
    Ftot = np.sum(Itotw)
    SPRk = np.sum(Icw)/Ftot
    Cx = np.sum(x*Itotw)/Ftot
    Cy = np.sum(y*Itotw)/Ftot

    VarItot = Itot + sb + sd**2 + sq**2
    VarItotw = VarItot*w

    # centroid shift calculation
    Gammax = np.sum(x*Icw)/Ftot-Cx*SPRk
    s = db/(1-db*SPRk)
    delta_Cx = s*Gammax
    Gammay = np.sum(y*Icw)/Ftot-Cy*SPRk
    delta_Cy = s*Gammay
    Gamma = np.sqrt(Gammax**2+Gammay**2)
    delta_C =  s * Gamma

    # uncertainty calculation
    delta_Cx_var = (np.sum(x**2*VarItotw))/Ftot**2 + (Cx/Ftot)**2*np.sum(VarItotw)
    delta_Cy_var = (np.sum(y**2*VarItotw))/Ftot**2 + (Cy/Ftot)**2*np.sum(VarItotw)

    # Lambdax = np.sum(x**2*Icw)/Ftot**2 + (Cx**2/Ftot)*SPRk - 2*SPRk*(Cx/Ftot)**2*np.sum(VarItotw)
    # # Lambday = np.sum(y**2*Icw)/Ftot**2 + (Cy**2/Ftot)*SPRk - 2*SPRk*(Cy/Ftot)**2*np.sum(VarItotw)
    # # Lambda = Lambdax + Lambday
    # delta_Cx_var_in = (delta_Cx_var -db*Lambdax)/(1.-db*SPRk)**2
    # # delta_Cy_var_in = (delta_Cy_var -db*Lambday)/(1.-db*SPRk)**2
    # Itot_in = Itot - db*Ic
    # delta_Cx_var_in2,delta_Cy_var_in2 =  barycenter_var(Itot_in,sb,sd,sq,mask=w)
    # print( (delta_Cx_var_in-delta_Cx_var_in2)/(delta_Cx_var_in+delta_Cx_var_in2)*0.5)
    # print( (delta_Cx_var-delta_Cx_var_in2)/(delta_Cx_var+delta_Cx_var_in2)*0.5,SPRk)
    #

    # we assume that the transit does not significantly change the uncertainty
    # delta_Cx_var_in = delta_Cx_var
    # delta_Cy_var_in = delta_Cy_var

    delta_C_sig = np.sqrt( 2* ( Gammax**2*delta_Cx_var + Gammay**2*delta_Cy_var) )/ Gamma

    delta_C_sig_1h_24c  = delta_C_sig / (12 * np.sqrt(24)) # random error re-scaled to 1h and 24 cameras

    ## print('delta_COB = ',delta_C)
    ## print('delta_COB_sig =',delta_C_sig_1h_24c)
    return delta_C,delta_C_sig_1h_24c,Gamma

#ftrack = open(DIRout+'track.log','w')
# We define a counter to store our data
counter = 0
# Now we can create the mask for getting only stars from P5 sample magnitude range
for i in range(nP):
    Pi = Pmin + i*binsize
    mask = (data[:, 2] >= Pi -binsize/2.) & (data[:, 2] < Pi + binsize/2.)
    n_star_p_bin[i] = mask.sum()
    print('P magnitude range: [%f - %f] number of stars in this range: %i' % (Pi -binsize/2.,Pi + binsize/2.,n_star_p_bin[i]))
    targets_P5 = data[mask, :]
    ID_target = ID[mask]

    j = ran_unique_int(n=nStar, interval=[0, targets_P5.shape[0] - 1])
    targets_P5 = targets_P5[j]
    ID_target = ID_target[j]
    # Now we obtain the x and y coordinates of the targets on the focal plane
    x_tar = targets_P5[:, 3]
    y_tar = targets_P5[:, 4]

    # We convert the coordinates of the randomly chosen targets to mm for obtaining the vignetting afterwards
    x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)

    n_t = len(j) # number of targets in the current range of magnitude

    def process_target(k): # main processing routine (work on a given target)

        # Define target ID
        ID_t = ID_target[k]
        if(verbose):
            print('Processing star ID %i' % ID_t)
        # First we compute the angle for obtaining the vignetting
        alpha = np.arctan(np.sqrt(x_tar_mm[k] ** 2 + y_tar_mm[k] ** 2) / 247.732)
        # Now we found the closest psf to every target (we will use this psf for the contaminants as well)
        s_d = (xpsf_pix - x_tar[k]) ** 2 + (ypsf_pix - y_tar[k]) ** 2
        # We store as a string the psf index of the closest psf
        psf_idx = np.argmin(s_d)
        psf_id = str(psf_idx + 1)

        # Now we we define the window (imagette) and find the coordinates of the target inside of it
        x_t_im, y_t_im, i0, j0 = window(x_tar[k], y_tar[k], 6, 6)
        # Then we obtain the offset between the center of the imagette and the center of the PSF
        offx = x_t_im - pxc[psf_idx]
        offy = y_t_im - pyc[psf_idx]
        if(gauss_psf):
           imagette = psf_gauss_int(x_t_im,y_t_im,gauss_width_x,gauss_width_y,6,6)
        else:
            psfbs_k = np.ascontiguousarray(psfbs[psf_idx])
        # Then we finally compute the imagette for the target by integrating the b-spline decomposition of the PSF
            imagette = spline2dbase.Spline2Imagette(psfbs_k, bsres, 6, 6, offx=offx, offy=offy)
        # ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
        # Then we can print the coordinates of the C.O.B.
        ## COBx, COBy = barycenter(imagette, subres=1)
        # Let's obtain the value of the reference flux after the integration time for the target star including the vignetting
        P_t = targets_P5[:,2][k]
        f_ref_t = reference_flux_target(P_t) * (np.cos(alpha) ** 2)
        # Let's obtain the flux per pixel of the target
        It = f_ref_t * imagette
        # Now it is time to find the contaminants sorrounding each target. We write the distance condition (10 pixels)
        dist = np.sqrt((x_star - x_tar[k]) ** 2 + (y_star - y_tar[k]) ** 2)
        # We define a useful mask now
        Delta_P = data[:, 2] - P_t
        m = (dist > 0) & (dist < distance_max)   & (Delta_P< Delta_P_max) & (data[:, 2]>0)
        n_c = m.sum()
        if(n_c> n_c_max):
            # too many contaminant stars, we keep only those for which Delta_P is smaller than
            # Delta_P_sorted[n_c_max], i.e.  Delta_P of the n_c_max-th contaminant star
             Delta_P_sorted = sort(Delta_P[m])
             m = m & (Delta_P <  Delta_P_sorted[n_c_max])
        ID_contaminants = ID[m] # IDs of the contaminant stars
        # We get the the index of all the contaminants now with the following line
        n = np.where(m)[0]
        # Now we find the magnitude of each contaminant
        m_contaminants = data[:, 2][n]
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
        offx_c = x_c_im - pxc[psf_idx]
        offy_c = y_c_im - pyc[psf_idx]
        # We define an array that will contain the 'imagettes' of every contaminant
        Ic = np.zeros((n_c, 6, 6))

        if(verbose):
            print('n_c = %i' % n_c)

        if(n_c==0):
            # no contaminant stars: skipping this non-relevant situation
            return None

        for o in range(0, n_c):
            if(gauss_psf):
                Ic[o] = psf_gauss_int(x_c_im[o],y_c_im[o],gauss_width_x,gauss_width_y,6,6)
            else:
                Ic[o] = spline2dbase.Spline2Imagette(psfbs_k, bsres , 6, 6, offx=offx_c[o], offy=offy_c[o])
            f_ref_c = reference_flux_contaminant(f_ref_t,m_contaminants[o],P_t)
            Ic[o] *= f_ref_c
            ## Ic = spline2dbase.Spline2ImagetteMulti(psfbs_k, bsres, 6, 6, offx_c, offy_c)
            # t3 = time.time()
            # print(t2-t1)
            # print(t3-t2)
        # Now we define an array with the contribution from all the stars to each pixel
        Ic_acc = np.sum(Ic, axis=0)

        # Let's compute the aperture of the target
        NSR1h, w_t = aperture(ft=It, fc=Ic_acc, sb=sb, sd=sd, sq=sq)

        # Now we store the nominal mask into a mask_key
        w_t_key = mask_to_bitmask(w_t)
        w_t_size = w_t.sum() # mask size

        # Now in this part of the code we present the calculations for the sprk of every contaminant as well as the
        # calculation of the SPR_crit.

        sprk, sprk_max, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)

        # We compute the critical SPR now
        SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr)

        n_bad = np.sum(sprk > SPR_crit)

        n_bad_wrong = np.sum(sprk > (SPR_crit/(1.-SPR_tot)) )

        # Now we have to know for which contaminant corresponds the highest sprk value we just found
        ind_sprk = np.argmax(sprk)
        Ic_max = Ic[ind_sprk]

        ID_c = ID_contaminants[ind_sprk]

        # computing the COB shift and its associated uncertainty
        Itot = It + Ic_acc
        delta_COB, delta_COB_sig_1h_24c,Gamma = cob_shift(Itot,Ic_max,w_t,dback)
        eta_cob = delta_COB*np.sqrt(td*ntr)/delta_COB_sig_1h_24c

        SPRk_10first = np.zeros(10)
        Gamma_10first = np.zeros(10)
        delta_COB_sig_10first = np.zeros(10)
        nsprmax = min(10,n_c)
        # sorting the SPRk by decreasing order and taking the 10 first values
        sprk_sorted_index  = (np.argsort(sprk)[::-1])
        SPRk_10first[0:nsprmax] = sprk[sprk_sorted_index[0:nsprmax]]
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            _, delta_COB_sig_l,Gamma_l = cob_shift(Itot,Ic[m],w_t,dback)
            Gamma_10first[l] = Gamma_l
            delta_COB_sig_10first[l] = delta_COB_sig_l

        if(verbose):
            print('SPR_tot=',SPR_tot)
            print('Mask size=',w_t_size)
            print('delta_COB = ',delta_COB)
            print('delta_COB_sig=',delta_COB_sig_1h_24c)
            print('eta_cob=',eta_cob)



########################################################################################################################
#                                   NOW THE EXTENDED MASK METHOD                                                       #
########################################################################################################################

        # First we create the extended mask given the nominal mask
        w_ext = extended_binary_mask(w_t, W=1)
        w_ext_key = mask_to_bitmask(w_ext)
        w_ext_size = w_ext.sum() # mask size

        # computing the COB shift and associated uncertainty for the extended mask
        delta_COB_ext, delta_COB_sig_1h_24c_ext,Gamma_ext  = cob_shift(Itot,Ic_max,w_ext,dback)
        eta_cob_ext = delta_COB_ext*np.sqrt(td*ntr)/delta_COB_sig_1h_24c_ext
        if(verbose):
            print('delta_COB_ext = ',delta_COB_ext)
            print('delta_COB_sig_ext=',delta_COB_sig_1h_24c_ext)
            print('eta_cob_ext=',eta_cob_ext)


        # Now we compute all the metrics associated with this mask. Let's begin with the NSR
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext)) / np.sum(It * w_ext)
        NSR_ext_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext


        # Then we compute the sprk over the extended mask for all the contaminants for a this target
        sprk_ext, sprk_max_ext, SPR_tot_ext= SPR(n_c=n_c, f_contaminant=Ic,
                                                             f_tot=(It + Ic_acc), w=w_ext)
        # Then compute the critical SPR
        SPR_crit_ext = spr_crit(dback=dback,SPR_tot=SPR_tot_ext,nsr=NSR_ext_1h, td=td, ntr=ntr)

        n_bad_ext = np.sum(sprk_ext>SPR_crit_ext)

        # observed transit depth given the maximum value of sprk with the extended mask
        delta_obs_ext_max = sprk_max_ext * dback

        # And now we compute the statistical significance over the extended mask
        eta_ext_max = sprk_max_ext * np.sqrt(td * ntr) * dback / NSR_ext_1h/(1.-SPR_tot_ext)

        # observed transit depth given the maximum value of sprk with the nominal mask
        delta_obs_ext = sprk_ext[ind_sprk] * dback

        # associated statistical significance
        eta_ext = sprk_ext[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h/(1.-SPR_tot_ext)

        # nb of positive detection with the extended mask
        n_det_ext = ( (sprk_ext > SPR_crit_ext) & (sprk > SPR_crit) & (sprk_ext > sprk)).sum()
        if(verbose):
            print('n_det_ext=',n_det_ext)

        SPRk_ext_10first = np.zeros(10)
        Gamma_ext_10first = np.zeros(10)
        delta_COB_sig_ext_10first = np.zeros(10)
        # taking the 10 first SPRk_ext, Gamma_ext and delta_COB_sig_ext (sorted by decreasing SPRk)
        SPRk_ext_10first[0:nsprmax] = sprk_ext[sprk_sorted_index[0:nsprmax]]
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            _, delta_COB_sig_ext_l,Gamma_ext_l = cob_shift(Itot,Ic[m],w_ext,dback)
            Gamma_ext_10first[l] = Gamma_ext_l
            delta_COB_sig_ext_10first[l] = delta_COB_sig_ext_l

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


        sprk_bray, sprk_max_bray, SPR_tot_bray  = SPR(n_c=n_c, f_contaminant=Ic,
                                                                 f_tot=(It + Ic_acc), w=w_bray)

        # We compute the critical SPR now
        SPR_crit_bray = spr_crit(dback=dback,SPR_tot=SPR_tot_bray, nsr=NSR_bray_1h, td=td, ntr=ntr)

        n_bad_bray = np.sum(sprk_bray>SPR_crit_bray)

########################################################################################################################
#                                          END OF TESTING Bray et al's ASSUMPTION                                      #
########################################################################################################################



        # Now we compute the secondary aperture for this contamninant with the highets sprk

        # We define the term that englobes the sigma of the target and the accumulated flux of the contaminants without
        # the contaminant of interest
        Itc_acc = Itot - Ic_max

        # Then we procedd to compute the secondary aperture
        NSR1h_c, w_c = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq)

        if(extend_2ndmask>0):
            w_c = extended_binary_mask(w_c, W=extend_2ndmask)
            NSR1h_c = ((10 ** 6) / (12 * np.sqrt(24))) * \
                      np.sqrt(np.sum((Itot + sb + sd ** 2 + sq ** 2) * w_c)) / np.sum(Ic_max * w_c)

        # Now we store this secondary mask in a mask key
        w_c_key = mask_to_bitmask(w_c)
        w_c_size = w_c.sum() # mask size

        # calculating the COG shift and the associated uncertainty
        delta_COB_c, delta_COB_sig_1h_24c_c,Gamma_c = cob_shift(Itot,Ic_max,w_c,dback)
        eta_cob_c = delta_COB_c*np.sqrt(td*ntr)/delta_COB_sig_1h_24c_c
        if(verbose):
            print('delta_COB_c = ',delta_COB_c)
            print('delta_COB_sig_c=',delta_COB_sig_1h_24c_c)
            print('eta_cob_c=',eta_cob_c)

        # We compute the flux over the secondary mask
        f_beb = Ic_acc * w_c
        f_t_c = It * w_c

        # We define the denominator of the spr_tot_c calculation
        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))

        # We compute spr_tot_c
        spr_tot_c = np.sum(Itc_acc * w_c) / f_tot_c

        # We compute now the delta_obs for the two apertures
        delta_obs_t = sprk_max * dback
        delta_obs_c = (1 - spr_tot_c) * dback

        # We compute now the statistical significances for a given transit event
        eta_t = sprk_max * np.sqrt(td * ntr) * dback / NSR1h/ (1. - SPR_tot)
        eta_c = np.sqrt(td * ntr) * dback / NSR1h_c

        d_c = np.sqrt(  (x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        m_c = m_contaminants[ind_sprk]
        if(verbose):
            print('Delta P is:', m_c - P_t)
            print('Distance between the target and contaminant',d_c)

            # if eta_t > eta_c:
            print('NSR_T is:', NSR1h)
            print('NSR_c is:', NSR1h_c)
            print('eta_t is:', eta_t)
            print('eta_c is:', eta_c)
            print('spr_t', sprk_max)
            print('spr_tot_c', spr_tot_c)

        SPRtot_c_10first = np.zeros(10)
        Gamma_c_10first = np.zeros(10)
        delta_COB_sig_c_10first = np.zeros(10)
        NSR1h_c_10first = np.zeros(10)
        w_c_key_10first = np.zeros(10)

        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            Itc_acc_l = Itot - Ic[m]

            # Then we procedd to compute the secondary aperture
            NSR1h_c_l, w_c_l = aperture(ft=Ic[m], fc=Itc_acc_l, sb=sb, sd=sd, sq=sq)
            NSR1h_c_10first[l] = NSR1h_c_l
            w_c_key_10first[l] = mask_to_bitmask(w_c_l)

            SPRtot_c_10first[l] = np.sum(Itc_acc_l * w_c_l) / np.sum(Itot*w_c_l)

            _, delta_COB_sig_c_l,Gamma_c_l = cob_shift(Itot,Ic[m],w_c_l,dback)
            Gamma_c_10first[l] = Gamma_c_l
            delta_COB_sig_c_10first[l] = delta_COB_sig_c_l


        save_info = np.array([ID_t,P_t,psf_idx,n_c,w_t_key,w_t_size,NSR1h,n_bad,SPR_crit,SPR_tot,ID_c,
                                 m_c,sprk_max,eta_t,delta_obs_t,delta_COB,delta_COB_sig_1h_24c,
                                 eta_cob,d_c,Gamma,0.,n_bad_wrong])
        save_info = np.append(save_info,SPRk_10first)
        save_info = np.append(save_info,Gamma_10first)
        save_info = np.append(save_info,delta_COB_sig_10first)

        save_info_2ndmask = np.array([ID_t,P_t,ID_c,m_c,NSR1h_c,w_c_key,
                                         w_c_size,eta_c,delta_obs_c,delta_COB_c,delta_COB_sig_1h_24c_c,eta_cob_c,spr_tot_c
                                         ,Gamma_c])
        save_info_2ndmask = np.append(save_info_2ndmask,SPRtot_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask,Gamma_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask,delta_COB_sig_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask,NSR1h_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask,w_c_key_10first)

        save_info_ext = np.array([ID_t,P_t,w_ext_key,w_ext_size,NSR_ext_1h,n_bad_ext,SPR_crit_ext,
                                     SPR_tot_ext,ID_c,m_c,eta_ext,delta_obs_ext,delta_COB_ext,
                                     delta_COB_sig_1h_24c_ext,eta_cob_ext,delta_obs_ext,eta_ext_max,delta_obs_ext_max,
                                     n_det_ext,Gamma_ext,0.])
        save_info_ext = np.append(save_info_ext,SPRk_ext_10first)
        save_info_ext = np.append(save_info_ext,Gamma_ext_10first)
        save_info_ext = np.append(save_info_ext,delta_COB_sig_ext_10first)

        save_info_bray =np.array([ID_t,P_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray])

        return save_info,save_info_2ndmask,save_info_ext,save_info_bray

    def process_target_wrapper(k):
        try:
            result = process_target(k)
        except:
            result = None
            traceback.print_exc()
            print(('an unkown error occur while processing star ID=%i (index=%i, magnitude interval %f-%f)')
                  % (ID_target[k],k,Pi-binsize/2.,Pi+binsize/2.))
        return result
    if(mp):  # multiprocessing mode
        nproc = min(n_t,nproc_max)
        pool = multiprocessing.Pool(nproc)
        results = pool.map(process_target_wrapper, [(t) for t in range(n_t)])
        pool.close()
        for k in range(n_t):
            if(results[k] is not None):
                save_info[counter] = results[k][0]
                save_info_2ndmask[counter] = results[k][1]
                save_info_ext[counter] = results[k][2]
                save_info_bray[counter] = results[k][3]
                counter += 1
    else: # sequential mode
        for k in range(n_t):
            result = process_target_wrapper(k)
            if(result is not None):
                save_info[counter] = result[0]
                save_info_2ndmask[counter] = result[1]
                save_info_ext[counter] = result[2]
                save_info_bray[counter] = result[3]
                counter += 1
            sys.stdout.flush()


    #     ftrack.write('%i %i %i %i\n' % (i,k,counter,ID_t))
    #     ftrack.flush()
    #
    # np.save(DIRout+'targets_P5.npy', save_info)
    # np.save(DIRout+'targets_P5_2ndmask.npy', save_info_2ndmask)
    # np.save(DIRout+'targets_P5_bray.npy', save_info_bray)
    # np.save(DIRout+'targets_P5_extended.npy', save_info_ext)
    print("%i targets processed so far" % counter)
    sys.stdout.flush()




# ftrack.close()
save_info = save_info[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]
# Now it is time to save the metrics of interest into a.npy file
np.save(DIRout+'targets_P5.npy', save_info)
# 0: ID_t
# 1: P_t
# 2: PSF index
# 3: n_c
# 4: w_t_key
# 5: w_t_size
# 6: NSR1h
# 7: n_bad
# 8: SPR_crit
# 9: SPR_tot
# 10: ID_c
# 11: P_c
# 12: spr max
# 13: eta_t
# 14: delta_obs_t
# 15: delta_COB
# 16: delta_COB_sig_1h_24c
# 17: eta_cob
# 18: d_c
# 19: Gamma
# 20:
# 21: n_bad_wrong (Victor's formula)
# 22-31: 10 first SPRk values
# 32-41: 10 first Gamma values
# 42-51: 10 first delta_COB_sig_1h_24c

np.save(DIRout+'targets_P5_2ndmask.npy', save_info_2ndmask)
# 0: ID_t
# 1: P_t
# 2: ID_c
# 3: P_c
# 4: NSR1h_c
# 5: w_c_key
# 6: w_c_size
# 7: eta_c
# 8: delta_obs_c
# 9: delta_COB_c
# 10: delta_COB_sig_1h_24c_c
# 11: eta_cob_c
# 12: spr_tot_c
# 13: Gamma_c
# 14-23: 10 first SPRtot_c values
# 24-33: 10 first Gamma_c values
# 34-43: 10 first delta_COB_sig_1h_24c
# 44-53: 10 first NSR1h_c
# 54-63: 10 first w_c_key



np.save(DIRout+'targets_P5_extended.npy', save_info_ext)
# 0: ID_t
# 1: P_t
# 2: w_ext_key
# 3: w_ext_size
# 4: NSR_ext_1h
# 5: n_bad_ext
# 6: SPR_crit_ext
# 7: SPR_tot_ext
# 8: ID_c
# 9: P_c
# 10: eta_ext
# 11: delta_obs_ext
# 12: delta_COB_ext
# 13: delta_COB_sig_1h_24c_ext
# 14: eta_cob_ext
# 15: delta_obs_ext
# 16: eta_ext_max
# 17: delta_obs_ext_max
# 18: n_det_ext
# 19: Gamma
# 20:
# 21-30: 10 first SPRk values
# 31-40: 10 first Gamma values
# 41-50: 10 first delta_COB_sig_1h_24c

np.save(DIRout+'targets_P5_bray.npy', save_info_bray)

fout = open(DIRout+'star_count.txt','w')
for i in range(nP):
    Pi = Pmin + i*binsize
    fout.write('%.2f %i\n' % (Pi,n_star_p_bin[i]))
fout.close()



