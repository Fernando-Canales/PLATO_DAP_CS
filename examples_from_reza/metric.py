import matplotlib.pyplot as plt # type: ignore
import spline2dbase # type: ignore
import scipy.signal # type: ignore
import h5py as h5py # type: ignore
import numpy as np # type: ignore
import sys
from fitting_psf import from_mm_2_pix, from_pix_2_mm, closest_psf, contaminants, reference_flux_target, \
    reference_flux_contaminant
from imagette import catalogue, list_psf, barycenter, gauss, window, ran_unique_int, ploting_imagettes, ploting_nsr, \
    ploting_nsr_s, \
    ploting_initial, psf_gauss_int
from NSR import spr_crit, aperture, NSRn, nsr_AGG, SPR, mask_to_bitmask, bitmask_to_mask, extended_binary_mask
from pylab import * # type: ignore
from tqdm import tqdm  # type: ignore # only used for convenience, to monitor progress in our  a calculation
import multiprocessing
import math
import time
import traceback
import os
# --------------------------------------------------
# CONFIGURATION PARAMETERS

verbose =  False
doplot = False
mp = True  # multiprocessing mode

CatalogueDIR = '/home/samadi/plato/share/catalogues/'
# CatalogueDIR = '/volumes/astro/sismo/general/plato/web/grids/catalogues/'
# CatalogueDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/catalogues_stars/'
# CatalogueDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues
DIRout = 'test40/'
# test 6 : v2 optimal extended mask (union of the secondary mask + nominal)
# test7 # v1 optimal extended mask: for each contaminant, take the extra pixels that maximize the significance
# test8 as test4 with a maximum distance of 10 pixel
# test 9   PSF_Focus_0mu_0.2pxdif.npz
# test 10  PSF_Focus_0mu_0.1pxdif.npz
# test 11 as 9 with optimal extended mask v1
# test 12 as 9 no RON no background no quantification noise
# test 13 factor sqrt(2) removed in equation for the centroid error
# test 14 optimal extended mask v1
# test 15 optimal extended mask v2 and extend_2ndmask=1
# test 16 optimal extended mask v1  : extension by 2 pixels
# test 17 low background level: 25
# test 18  PSF_Focus_0mu_0.1pxdif.npz
# test 19  PSF_Focus_0mu_0.2pxdif.npz with bsres = 10
# test 20 Pmax = 16
# test 21 SFP_DR3_20230926_gr1.npy
# test 22 SFP_DR3_20230926_gr2.npy
# test 23 SFP_DR3_20230926_gr3.npy
# test 24 SFP_DR3_20230926_gr4.npy
# test 25 PSF_Focus_0mu_0.2pxdif_N4000K.npz
# test 26 PSF_Focus_0mu_0.2pxdif_N6500K.npz
# test 27: target: SF_Focus_0mu_0.2pxdif.npz contaminant: PSF_Focus_0mu_0.2pxdif_N4000K.npz
# test 28: FM_Leopold7_02470_IAS_inversion_240927_iter2.npz warning  bsres = 5
# test 29: FM_Joup_02182_IAS_inversion_240927_iter2   bsres = 5
# test 30: gaussian PSF 0.5 px width
# test 31: gaussian PSF 1 px width
# test 32: optimal extended mask v3 and extension by 2 pixels
# test 33: optimal extended mask v3 and extension by 1 pixel
# test 34: optimal extended mask v3 and extension by 4  pixels (full imagette)
# test 35: NFP1
# test 36: RON 120e-
# test 37: RON 120e-  gaussian sigma = 1.5. px  win_size = 12
# test 38 RON 120e-  gaussian sigma = 1.5. px  win_size = 12 Kepler field
# test 39: as 13 with rho = 5% (sub-optimal mask)
# test 40: as 13, with additionnal informations


PSFFileName = 'PSF_Focus_0mu_0.2pxdif.npz'
# PSFFileName = 'PSF_Focus_0mu_0.1pxdif.npz'
# PSFFileName = 'PSF_Focus_0mu_0.2pxdif_N4000K.npz'
# PSFFileName = 'PSF_Focus_0mu_0.2pxdif_N6500K.npz'
# PSFFileName = 'FM_Leopold7_02470_IAS_inversion_240927_iter2.npz'
# PSFFileName = '/home/samadi/plato/share/test_house/FM/Joup/02182_IAS/inversion_240927_iter2/FM_Joup_02182_IAS_inversion_240927_iter2.npz'

PSFFileNameC = None #  'PSF_Focus_0mu_0.2pxdif_N4000K.npz' # PSF used for the contaminant stars


# CatalogueFileName = 'SFP_DR3_20220831.npy'
# CatalogueFileName = 'LOPN1_DR3_20241011_gr0.npy' #
CatalogueFileName = 'SFP_DR3_20230101.npy'
# CatalogueFileName = 'gaia_dr3_kepler_241022.npy'
# CatalogueFileName =  'SFP_DR3_20230926_gr1.npy'
# CatalogueFileName =  'SFP_DR3_20230926_gr2.npy'
# CatalogueFileName =  'SFP_DR3_20230926_gr3.npy'
# CatalogueFileName =  'SFP_DR3_20230926_gr4.npy'

gauss_psf = False
gauss_width_x =  0.5 #   1.5 px -> encircle diameter 95% ~ in 2.45 sigma
gauss_width_y =   0.5
# encircle energy (EE) 95% in  2.45 sigma

bsres =  20  # resolution adopted for the b-spline decomposition

win_size =  6 # window size

nproc_max = 10  # number of processors

Pmin = 10
Pmax = 13
binsize = 0.5
nStar = 1000
Delta_P_max = 15.
Pc_max = 99 #  16.
distance_max = 7
n_c_max = 300
mask_extend =  1  # the number by which the nominal mask is extended (typically: 1)
extend_2ndmask = 0  # the number by which the secondary mask is extended (default: 0)
opt_ext_mask = 0 #  version of extended mask calculation: 0-> no optimization, >0 -> optimization
rho_mask = 0. # parameter to compte a sub-optimal mask

# Parameters for the NSR
sb = (45. * 21)  # Background noise form zodiacal light in units of e-/px after multiplying by the integration time (poisson noise)
# sb = (25. * 21)  # Background noise form zodiacal light in units of e-/px after multiplying by the integration time (poisson noise)

sd = 50.2  # Overall detector noise (including readout at beginning of life, smearing and dark current) in units of e- rms/px
sq = 7.2  # Quantization noise in units of e-rms/px

# Parameters for the eclipsing binaries
dback = 132000   # mean transit depth in ppm (Victor's value : 85 000)
td = 6.72 *0.46**2  # mean transit duration in hours  changed to have gamma = 0.46  assuming <td>=6.72h and <dback>= 132 000 ppm
ntr = 3  # number of transits

seed = 300

# --------------------------------------------------
#
_, fileExtension = os.path.splitext(CatalogueFileName)
if(fileExtension == '.coo'):
    data = np.genfromtxt(CatalogueDIR + CatalogueFileName,usecols=(0,1,2,3,4))
else:
    data = np.load(CatalogueDIR + CatalogueFileName)
psfdata = np.load(PSFFileName)

psfbs = psfdata['psfbs']
pxc = psfdata['pxc']
pyc = psfdata['pyc']
xpsf_pix = psfdata['xpsf_pix']
ypsf_pix = psfdata['ypsf_pix']

if( (PSFFileNameC is not None) and (PSFFileNameC != PSFFileName)): psfdata_c = np.load(PSFFileNameC) # contaminants have a PSF different to the target PSF
else: psfdata_c = psfdata # same PSF for the target and the contaminants

psfbs_c = psfdata_c['psfbs']
pxc_c = psfdata_c['pxc']
pyc_c = psfdata_c['pyc']
xpsf_pix_c = psfdata_c['xpsf_pix']
ypsf_pix_c = psfdata_c['ypsf_pix']

# Define an ID for every star
ID = np.arange(0, data.shape[0])

# Now we save the x and y coordinates on the focal plane of all the stars in the catalogue
x_star = data[:, 3]
y_star = data[:, 4]

nP = int(round((Pmax - Pmin) / binsize + 1))


file_out = open(DIRout + 'metrics_reza.txt', 'w')
# Define a numpy array for saving the metrics of interest (Target ID, magnitude, N_bad, etc.) (Is hard-coded now)
save_info = np.zeros((nStar * nP, 106))

# for the secondary mask
save_info_2ndmask = np.zeros((nStar * nP, 85))

# The same for the extended mask
save_info_ext = np.zeros((nStar * nP, 81))

# The same for bray's et al. assu#mption of using 2 x 2 masks
save_info_bray = np.zeros((nStar * nP, 6))

# Now we choose the random targets using Réza's function
np.random.seed(seed)
n_star_p_bin = np.zeros((nP))


def cob_shift(Itot, Ic, w, dback):
    # compute the COB shift and the associated uncertainty
    db = dback * 1e-6
    x = (np.arange(0, Ic.shape[1]) + 0.5)
    y = (np.arange(0, Ic.shape[0]) + 0.5)
    x, y = np.meshgrid(x, y)

    Icw = Ic * w
    Itotw = Itot * w
    Ftot = np.sum(Itotw)
    SPRk = np.sum(Icw) / Ftot
    Cx = np.sum(x * Itotw) / Ftot
    Cy = np.sum(y * Itotw) / Ftot

    VarItot = Itot + sb + sd ** 2 + sq ** 2
    VarItotw = VarItot * w

    # centroid shift calculation
#    Gammax = np.sum(x * Icw) / Ftot - Cx * SPRk
    Gammax = np.sum( (x - Cx) * Icw) / Ftot  # centered formula
    s = db / (1 - db * SPRk)
#    Gammay = np.sum(y * Icw) / Ftot - Cy * SPRk
    Gammay = np.sum( (y-Cy) * Icw) / Ftot # centered formula
    Gamma = np.sqrt(Gammax ** 2 + Gammay ** 2)
    delta_C = s * Gamma

    # uncertainty calculation
    # wrong Bryson 's forumula
    # delta_Cx_var = (np.sum(x ** 2 * VarItotw)) / Ftot ** 2 + (Cx / Ftot) ** 2 * np.sum(VarItotw)
    # delta_Cy_var = (np.sum(y ** 2 * VarItotw)) / Ftot ** 2 + (Cy / Ftot) ** 2 * np.sum(VarItotw)
    #

    delta_Cx_var = np.sum( (x-Cx) ** 2 * VarItotw) / Ftot ** 2
    delta_Cy_var = np.sum( (y-Cy) ** 2 * VarItotw) / Ftot ** 2


    if(Gamma>0.): delta_C_sig = np.sqrt( (Gammax ** 2 * delta_Cx_var + Gammay ** 2 * delta_Cy_var)) / Gamma
    else: delta_C_sig = 1e99
    delta_C_sig_1h_24c = delta_C_sig / (12 * np.sqrt(24))  # random error re-scaled to 1h and 24 cameras

    return delta_C, delta_C_sig_1h_24c, Gamma


def cal_opt_extended_mask_1(nmask,It,Ic,flag,background,ron,verbose=False,doplot=False):
    '''
    build an optimal extended mask by joining together (union)
    the optimal extended mask of each contaminant that can (potentially) generate a FP

    the optimal extended mask of each contaminant is calculating by gathering together the pixels in the ring
    until we reach an optimal significance , prior these pixels are sorted in decreasing SPR

    '''
    Npix = nmask.shape[0]
    emask = np.zeros((Npix,Npix),dtype=bool)
    emask[:,:] = nmask
    we = np.array(extended_binary_mask(nmask,mask_extend),dtype=bool)
    extra_pixel = we & (emask==False) # the pixels in the ring (called 'extra' pixel in the following)
    extra_pixel = extra_pixel.flatten()
    extra_pixel_idx = np.where(extra_pixel)[0]
    ne = len(extra_pixel_idx)
    idx = np.where(flag)[0]
    Ic_acc = np.sum(Ic, axis=0).flatten()
    It_f = It.flatten()
    wt = np.array(nmask,dtype=bool).flatten()
    for i in idx: # loop over the flagged contaminant stars
        w = np.zeros(Npix*Npix,dtype=bool)
        # for each extra pixel, compute its SPR
        spr = np.zeros(ne)
        k = 0
        for j in extra_pixel_idx:
            spr[k] = Ic[i].flatten()[j]/(Ic_acc[j] + It_f[j])
            k+= 1
        # sorting the extra pixel in decreasing SPR and storing their indexes
        j = np.argsort(spr)[::-1]
        eta_agg = np.zeros(ne)
        for k in range(ne):
            l = j[0:k+1]
            w[extra_pixel_idx[l]] = True
            w = w | wt
            ftot_k = np.sum( (Ic_acc + It_f)*w)
            spr_k = np.sum(  Ic[i].flatten()*w )/ftot_k
            sprtot_k = np.sum(Ic_acc*w)/ftot_k
            ft_k = np.sum(It_f*w)
            nsr_k = math.sqrt( np.sum( w* ( Ic_acc +  It_f  + ron**2 * background) ) )/ft_k
            eta_agg[k] = spr_k/((1-sprtot_k)*nsr_k) # aggregated significance (arbitrary unit)
        m = np.argmax(eta_agg)
        w = np.zeros(Npix * Npix, dtype=bool)
        w[extra_pixel_idx[j][0:m+1]] = True
        w = np.reshape(w,(Npix,Npix))
        emask = emask | w
        if(verbose):
            print('number of extra pixel kept:', (m + 1) )
            print('fraction of the extra pixesl:', i, float(m + 1) / float(ne))
        if(doplot):
            figure(0)
            clf()
            plot(spr[j])

            figure(1)
            clf()
            plot(eta_agg)

            figure(2)
            clf()
            imshow(w,origin='lower')
            title('Extra pixels kept')

            figure(3)
            clf()
            imshow(emask,origin='lower')
            title('emask')

            figure(4)
            clf()
            imshow(nmask,origin='lower')
            title('nmaxk')

            figure(5)
            clf()
            imshow(we & (wt==False).reshape((Npix,Npix)) ,origin='lower')
            title('Extra pixels')
            show(block=True)


    return emask


def cal_opt_extended_mask_2(W_t,It,Ic,flag):
    '''
    build an optimal extended mask by joining together (union) the secondary mask of each contaminant star that can create a FP
    '''
    W_ext = np.zeros(W_t.shape,dtype=bool)
    W_ext[:,:] = W_t

    idx = np.where(flag)[0]
    for i in idx: # loop over the flagged contaminant stars
        It_i = Ic[i]
        Ic_acc_i = np.sum(Ic, axis=0) - It_i + It
        # Then we proceed to compute the secondary aperture for the
        _, W_c = aperture(ft=It_i, fc=Ic_acc_i, sb=sb, sd=sd, sq=sq,window_size=win_size,rho=rho_mask)
        W_c = np.array(W_c,dtype=bool)
        W_ext = W_ext | W_c
    return W_ext




def cal_opt_extended_mask_3(nmask, W=1):
    '''
    build an extended mask where only the pixel external to the nominal mask are kept

    '''
    ny, nx = nmask.shape
    emask = np.zeros((ny, nx))
    w_ext =  extended_binary_mask(nmask, W=W)
    m = (np.array(w_ext,dtype=bool)==True ) & (np.array(nmask,dtype=bool)==False) # the pixels in the ring (called 'extra' pixel in the following)
    emask[m] = 1.

    return emask



# ftrack = open(DIRout+'track.log','w')
# We define a counter to store our data
counter = 0
# Now we can create the mask for getting only stars from P5 sample magnitude range
for i in range(nP):
    Pi = Pmin + i * binsize
    mask = (data[:, 2] >= Pi - binsize / 2.) & (data[:, 2] < Pi + binsize / 2.)
    n_star_p_bin[i] = mask.sum()
    print('P magnitude range: [%f - %f] number of stars in this range: %i' % (
    Pi - binsize / 2., Pi + binsize / 2., n_star_p_bin[i]))
    targets_P5 = data[mask, :]
    ID_target = ID[mask]

    j = ran_unique_int(n=nStar, interval=[0, targets_P5.shape[0] - 1])
    print('number of stars processed in this range: %i' % (len(j)))
    targets_P5 = targets_P5[j]
    ID_target = ID_target[j]
    # Now we obtain the x and y coordinates of the targets on the focal plane
    x_tar = targets_P5[:, 3]
    y_tar = targets_P5[:, 4]

    # We convert the coordinates of the randomly chosen targets to mm for obtaining the vignetting afterwards
    x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar, y_tar)

    n_t = len(j)  # number of targets in the current range of magnitude
    #n_t = 1

    def process_target(k):  # main processing routine (work on a given target)

        # Define target ID
        ID_t = ID_target[k]
        if verbose:
            print('Processing star ID %i' % ID_t)
        # First we compute the angle for obtaining the vignetting
        alpha = np.arctan(np.sqrt(x_tar_mm[k] ** 2 + y_tar_mm[k] ** 2) / 247.732)
        # Now we found the closest psf to every target (we will use this psf for the contaminants as well)
        s_d = (xpsf_pix - x_tar[k]) ** 2 + (ypsf_pix - y_tar[k]) ** 2
        # We store as a string the psf index of the closest psf
        psf_idx = np.argmin(s_d)
        psf_id = str(psf_idx + 1)

        s_d_c = (xpsf_pix_c - x_tar[k]) ** 2 + (ypsf_pix_c - y_tar[k]) ** 2
        psf_idx_c = np.argmin(s_d_c)

        # Now we we define the window (imagette) and find the coordinates of the target inside of it
        x_t_im, y_t_im, i0, j0 = window(x_tar[k], y_tar[k], win_size, win_size)
        # Then we obtain the offset between the center of the imagette and the center of the PSF
        offx = x_t_im - pxc[psf_idx]
        offy = y_t_im - pyc[psf_idx]
        if gauss_psf:
            imagette = psf_gauss_int(x_t_im, y_t_im, gauss_width_x, gauss_width_y, win_size, win_size)
        else:
            psfbs_k = np.ascontiguousarray(psfbs[psf_idx])
            psfbs_k_c = np.ascontiguousarray(psfbs_c[psf_idx_c])

            # Then we finally compute the imagette for the target by integrating the b-spline decomposition of the PSF
            imagette = spline2dbase.Spline2Imagette(psfbs_k, bsres, win_size, win_size, offx=offx, offy=offy)

        # B = barycenter(imagette)
        # print(psf_idx,pxc[psf_idx],pyc[psf_idx])
        # print(xpsf_pix[psf_idx],ypsf_pix[psf_idx])
        # print(i0,j0,x_tar[k], y_tar[k],B,x_t_im,y_t_im,x_t_im-B[0] ,y_t_im - B[1])
        # print(offx,offy)
        # ploting_initial(2, 1, psf, imagette, i='PSF', j='Target')
        # Then we can print the coordinates of the C.O.B.
        ## COBx, COBy = barycenter(imagette, subres=1)
        # Let's obtain the value of the reference flux after the integration time for the target star including the vignetting
        P_t = targets_P5[:, 2][k]
        f_ref_t = reference_flux_target(P_t) * (np.cos(alpha) ** 2)
        # Let's obtain the flux per pixel of the target
        It = f_ref_t * imagette
        # Now it is time to find the contaminants sorrounding each target. We write the distance condition (10 pixels)
        dist = np.sqrt((x_star - x_tar[k]) ** 2 + (y_star - y_tar[k]) ** 2)
        # We define a useful mask now
        Delta_P = data[:, 2] - P_t

        m = (dist > 0) & (dist < distance_max) & (Delta_P < Delta_P_max) & (data[:, 2] > 0) & (data[:, 2]<Pc_max)
        n_c = m.sum()
        if n_c > n_c_max:
            # too many contaminant stars, we keep only those for which Delta_P is smaller than
            # Delta_P_sorted[n_c_max], i.e.  Delta_P of the n_c_max-th contaminant star
            Delta_P_sorted = sort(Delta_P[m]) # type: ignore
            m = m & (Delta_P < Delta_P_sorted[n_c_max])
        ID_contaminants = ID[m]  # IDs of the contaminant stars
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
        offx_c = x_c_im - pxc_c[psf_idx_c]
        offy_c = y_c_im - pyc_c[psf_idx_c]
        # We define an array that will contain the 'imagettes' of every contaminant
        Ic = np.zeros((n_c, win_size, win_size)) # image associated with each contaminant, computed with the contaminant PSF
        if(psfdata_c is not psfdata):
            Ic_ap = np.zeros((n_c, win_size, win_size)) #  image associated with each contaminant, computed with the PSF of the target, used for the aperture calculation
            offx_c_ap = x_c_im - pxc_c[psf_idx]
            offy_c_ap = y_c_im - pyc_c[psf_idx]
        else:
            Ic_ap = Ic
            offx_c_ap = offx_c
            offy_c_ap = offy_c

        if verbose:
            print('n_c = %i' % n_c)

        if (n_c == 0):
            # no contaminant stars: skipping this non-relevant situation
            return None

        for o in range(0, n_c):
            if gauss_psf:
                Ic[o] = psf_gauss_int(x_c_im[o], y_c_im[o], gauss_width_x, gauss_width_y, win_size, win_size)
                if (psfdata_c is not psfdata):
                    Ic_ap[o] = Ic[o]
            else:
                Ic[o] = spline2dbase.Spline2Imagette(psfbs_k_c, bsres, win_size, win_size, offx=offx_c[o], offy=offy_c[o])
                if (psfdata_c is not psfdata): # contaminant with the target PSF, only for the aperture calculation
                    Ic_ap[o] = spline2dbase.Spline2Imagette(psfbs_k, bsres, win_size, win_size, offx=offx_c_ap[o], offy=offy_c_ap[o])
            f_ref_c = reference_flux_contaminant(f_ref_t, m_contaminants[o], P_t)
            Ic[o] *= f_ref_c
            if (psfdata_c is not psfdata):
                Ic_ap[o] *= f_ref_c
            ## Ic = spline2dbase.Spline2ImagetteMulti(psfbs_k, bsres, win_size, win_size, offx_c, offy_c)
            # t3 = time.time()
            # print(t2-t1)
            # print(t3-t2)
        # Now we define an array with the contribution from all the stars to each pixel
        Ic_acc = np.sum(Ic, axis=0)
        Ic_acc_ap = np.sum(Ic_ap, axis=0) # for the aperture calculation only

        # Now we define the total flux (target and all contaminants)
        f_tot = It + Ic_acc

        # Let's compute the aperture of the target
        # calculation based on the target PSF assuming the same PSF for the contaminants
        NSR1h_0, w_t = aperture(ft=It, fc=Ic_acc_ap, sb=sb, sd=sd, sq=sq,window_size=win_size,rho=rho_mask)
        # NSR1h_0 : theoretical NSR ->  calculation based on the target PSF
        NSR1h = ((10 ** 6) / (12 * np.sqrt(24))) * \
                  np.sqrt(np.sum((f_tot + sb + sd ** 2 + sq ** 2) * w_t)) / np.sum(It * w_t)
        if(verbose): print("compare n-mask: theoretical NSR = %f ; effective NSR = %f" % (NSR1h_0,NSR1h))
        # Now we store the nominal mask into a mask_key
        if(win_size<=8): w_t_key = mask_to_bitmask(w_t)
        else: w_t_key = 0
        w_t_size = w_t.sum()  # mask size

        # Now in this part of the code we present the calculations for the sprk of every contaminant as well as the
        # calculation of the SPR_crit.

        sprk, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_t)
        # sprk, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_t)

        # We compute the critical SPR now
        SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr)

        n_bad = np.sum(sprk > SPR_crit)
        if(verbose): print('number of bad = %i' % n_bad)

        n_bad_wrong = np.sum(sprk > (SPR_crit / (1. - SPR_tot)))

        # Now we get the index of the contaminant star with the highest sprk value with respect to the nominal mask
        ind_sprk = np.argmax(sprk)
        Ic_max = Ic[ind_sprk]
        Ic_max_ap = Ic_ap[ind_sprk]

        sprk_max = sprk[ind_sprk]


        ID_c = ID_contaminants[ind_sprk]

        # computing the COB shift and its associated uncertainty
        Itot = It + Ic_acc
        delta_COB, delta_COB_sig_1h_24c, Gamma = cob_shift(Itot, Ic_max, w_t, dback)
        eta_cob = delta_COB * np.sqrt(td * ntr) / delta_COB_sig_1h_24c

        SPRk_10first = np.zeros(10)
        Gamma_10first = np.zeros(10)
        delta_COB_sig_10first = np.zeros(10)
        eta_COB_10first = np.zeros(10)
        eta_10first = np.zeros(10)
        nsprmax = min(10, n_c)
        # sorting the SPRk by decreasing order and taking the 10 first values
        sprk_sorted_index = (np.argsort(sprk)[::-1])
        SPRk_10first[0:nsprmax] = sprk[sprk_sorted_index[0:nsprmax]]
        IDs_10first = np.zeros(10)
        Delta_x_10first = np.zeros(10)
        Delta_y_10first = np.zeros(10)
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            eta_10first[l] = sprk[m] * np.sqrt(td * ntr) * dback / NSR1h / (1. - SPR_tot)
            _, delta_COB_sig_l, Gamma_l = cob_shift(Itot, Ic[m], w_t, dback)
            Gamma_10first[l] = Gamma_l
            delta_COB_sig_10first[l] = delta_COB_sig_l
            eta_COB_10first[l] = _ * np.sqrt(td * ntr) / delta_COB_sig_l
            IDs_10first[l] = ID_contaminants[m]
            Delta_x_10first[l] = x_c[m] -  x_tar[k]
            Delta_y_10first[l] = y_c[m] -  y_tar[k]

        if verbose:
            print('SPR_tot=', SPR_tot)
            print('Mask size=', w_t_size)
            print('delta_COB = ', delta_COB)
            print('delta_COB_sig=', delta_COB_sig_1h_24c)
            print('eta_cob=', eta_cob)
        # figure(0)
        # clf()
        # imshow(w_t)
        # figure(1)
        # clf()
        # imshow(Itot)
        # stop
        ########################################################################################################################
        #                                   NOW THE EXTENDED MASK METHOD                                                       #
        ########################################################################################################################
        # First we create the extended mask given the nominal mask
        if(opt_ext_mask ==1 ): w_ext = cal_opt_extended_mask_1(w_t,It,Ic,sprk > SPR_crit,sb,sd,verbose=verbose,doplot=doplot)
        elif(opt_ext_mask == 2 ): w_ext = cal_opt_extended_mask_2(w_t,It,Ic,sprk > SPR_crit)
        elif(opt_ext_mask == 3 ): w_ext = cal_opt_extended_mask_3(w_t, W=mask_extend)
        else: w_ext = extended_binary_mask(w_t, W=mask_extend)
        if(verbose):
            print('n_bad = %i'  % n_bad)
            print('extended maks size = %i'  % w_ext.sum())
            # figure(0)
            # clf()
            # imshow(w_t,origin='lower')
            # figure(1)
            # clf()
            # imshow(w_ext,origin='lower')
            # show(block=True)

        if(win_size<=8): w_ext_key = mask_to_bitmask(w_ext)
        else: w_ext_key = 0
        w_ext_size = w_ext.sum()  # mask size

        # computing the COB shift and associated uncertainty for the extended mask
        delta_COB_ext, delta_COB_sig_1h_24c_ext, Gamma_ext = cob_shift(Itot, Ic_max, w_ext, dback)
        eta_cob_ext = delta_COB_ext * np.sqrt(td * ntr) / delta_COB_sig_1h_24c_ext
        if verbose:
            print('delta_COB_ext = ', delta_COB_ext)
            print('delta_COB_sig_ext=', delta_COB_sig_1h_24c_ext)
            print('eta_cob_ext=', eta_cob_ext)

        # Now we compute all the metrics associated with this mask. Let's begin with the NSR
        NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_ext)) / np.sum(It * w_ext)
        NSR_ext_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext

        # mask with only the extra pixels of the extended mask
        Delta_W_ext =   np.array(w_ext,dtype=bool)  &  (np.array(w_t,dtype=bool) == False)
        NSR_chi = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * Delta_W_ext)) / np.sum(It * w_t)
        NSR_chi = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_chi
        rho = np.sum((It + Ic_acc) * w_t) / np.sum((It + Ic_acc) * w_ext) # rho = Ftot_nom / Ftot_ext


        # Then we compute the sprk over the extended mask for all the contaminants for a this target
        sprk_ext, SPR_tot_ext = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_ext)

        sprk_max_ext = np.max(sprk_ext)

        # Then compute the critical SPR
        SPR_crit_ext = spr_crit(dback=dback, SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h, td=td, ntr=ntr)

        n_bad_ext = np.sum(sprk_ext > SPR_crit_ext)

        # observed transit depth given the maximum value of sprk with the extended mask
        delta_obs_ext_max = sprk_max_ext * dback

        # And now we compute the statistical significance over the extended mask
        eta_ext_max = sprk_max_ext * np.sqrt(td * ntr) * dback / NSR_ext_1h / (1. - SPR_tot_ext)

        # observed transit depth given the maximum value of sprk with the nominal mask
        delta_obs_ext = sprk_ext[ind_sprk] * dback

        # associated statistical significance
        eta_ext = sprk_ext[ind_sprk] * np.sqrt(td * ntr) * dback / NSR_ext_1h / (1. - SPR_tot_ext)

        # nb of positive detection with the extended mask
        n_det_ext = ((sprk_ext > SPR_crit_ext) & (sprk > SPR_crit) & (sprk_ext > sprk)).sum()
        if verbose:
            print('n_det_ext=', n_det_ext)

        SPRk_ext_10first = np.zeros(10)
        Gamma_ext_10first = np.zeros(10)
        delta_COB_sig_ext_10first = np.zeros(10)
        eta_COB_ext_10first = np.zeros(10)
        eta_ext_10first = np.zeros(10)
        eta_dtd_10first = np.zeros(10) # significance of the differential transit depth
        # taking the 10 first SPRk_ext, Gamma_ext and delta_COB_sig_ext (sorted by decreasing SPRk)
        SPRk_ext_10first[0:nsprmax] = sprk_ext[sprk_sorted_index[0:nsprmax]]
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            eta_ext_10first[l] = sprk_ext[m] * np.sqrt(td * ntr) * dback / NSR_ext_1h / (1. - SPR_tot_ext)
            _, delta_COB_sig_ext_l, Gamma_ext_l = cob_shift(Itot, Ic[m], w_ext, dback)
            Gamma_ext_10first[l] = Gamma_ext_l
            delta_COB_sig_ext_10first[l] = delta_COB_sig_ext_l
            eta_COB_ext_10first[l] = _ * np.sqrt(td * ntr) / delta_COB_sig_ext_l
            eta_dtd_10first[l] = dback * (sprk_ext[m] - sprk[m]) * np.sqrt(td * ntr) / (
               (1-SPR_tot) * math.sqrt (  (rho-1.)**2 * NSR1h**2  + rho**2 * NSR_chi**2 )           )

        ########################################################################################################################
        #                                     END OF THE EXTENDED MASK METHOD                                                  #
        ########################################################################################################################

        ########################################################################################################################
        #                                   TESTING  J.C. Bray et al's ASSUMPTION OF A 2 x 2 MASK                              #
        ########################################################################################################################

        # The mask has to contain the 4 pixels around the center,
        w_bray = np.zeros((win_size, win_size))
        w_bray[2:4, 2:4] = 1

        NSR_bray = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * w_bray)) / np.sum(It * w_bray)
        NSR_bray_1h = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_bray

        sprk_bray,  SPR_tot_bray = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_bray)

        sprk_max_bray = np.max(sprk_bray)

        # We compute the critical SPR now
        SPR_crit_bray = spr_crit(dback=dback, SPR_tot=SPR_tot_bray, nsr=NSR_bray_1h, td=td, ntr=ntr)

        n_bad_bray = np.sum(sprk_bray > SPR_crit_bray)

        ########################################################################################################################
        #                                          END OF TESTING Bray et al's ASSUMPTION                                      #
        ########################################################################################################################


        ########################################################################################################################
        #                                          SECONDARY MASK
        ########################################################################################################################

        # Now we compute the secondary aperture for this contamninant with the highets sprk

        # We define the term that englobes the sigma of the target and the accumulated flux of the contaminants without
        # the contaminant of interest
        Itc_acc = Itot - Ic_max
        Itc_acc_ap = Itot - Ic_max_ap # for the aperture calculation only

        # Then we procedd to compute the secondary aperture
        # calculation based on the target PSF assuming the same PSF for the contaminants
        NSR1h_c_0, w_c = aperture(ft=Ic_max_ap, fc=Itc_acc_ap, sb=sb, sd=sd, sq=sq,window_size=win_size,rho=rho_mask)
        # NSR1h_c_0:  theoretical NSR -> calculation based on the target PSF

        if extend_2ndmask > 0:
            w_c = extended_binary_mask(w_c, W=extend_2ndmask)
        NSR1h_c = ((10 ** 6) / (12 * np.sqrt(24))) * \
                  np.sqrt(np.sum((Itot + sb + sd ** 2 + sq ** 2) * w_c)) / np.sum(Ic_max * w_c)
        if(verbose): print("compare s-mask: theoretical NSR = %f ; effective NSR = %f" % (NSR1h_c_0,NSR1h_c))

        # Now we store this secondary mask in a mask key
        if(win_size<=8): w_c_key = mask_to_bitmask(w_c)
        else: w_c_key = 0
        w_c_size = w_c.sum()  # mask size

        # calculating the COG shift and the associated uncertainty
        delta_COB_c, delta_COB_sig_1h_24c_c, Gamma_c = cob_shift(Itot, Ic_max, w_c, dback)
        if(w_c_size==1): delta_COB_c = 0. # no centroid shift in that cas, we impose zero
        eta_cob_c = delta_COB_c * np.sqrt(td * ntr) / delta_COB_sig_1h_24c_c
        if verbose:
            print('delta_COB_c = ', delta_COB_c)
            print('delta_COB_sig_c=', delta_COB_sig_1h_24c_c)
            print('eta_cob_c=', eta_cob_c)

        # We compute the flux over the secondary mask
        f_beb = Ic_acc * w_c
        f_t_c = It * w_c

        # We define the denominator of the spr_tot_c calculation
        f_tot_c = (np.sum(f_t_c) + np.sum(f_beb))

        # We compute spr_tot_c
        spr_tot_c = np.sum(Itc_acc * w_c) / f_tot_c
        # sprk_c, _ = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=w_c)

        # We compute now the delta_obs for the two apertures
        delta_obs_t = sprk_max * dback
        delta_obs_c = (1 - spr_tot_c) * dback

        # We compute now the statistical significances for a given transit event
        eta_t = sprk_max * np.sqrt(td * ntr) * dback / NSR1h / (1. - SPR_tot)
        eta_c = np.sqrt(td * ntr) * dback / NSR1h_c

        d_c = np.sqrt((x_t_im - x_c_im[ind_sprk]) ** 2 + (y_t_im - y_c_im[ind_sprk]) ** 2)
        m_c = m_contaminants[ind_sprk]
        if verbose:
            print('Delta P is:', m_c - P_t)
            print('Distance between the target and contaminant', d_c)

            # if eta_t > eta_c:
            print('Nbad for this target =', n_bad)
            print('NSR_T is:', NSR1h)
            print('NSR_c is:', NSR1h_c)
            print('eta_t is:', eta_t)
            print('eta_c is:', eta_c)
            print('spr_t is:', sprk[ind_sprk])
            # print('spr_c is:', sprk_c[ind_sprk])
            print('spr_tot_c is:', spr_tot_c)
            print('abs_cob', delta_COB)
            print('abs_cob_ext:', delta_COB_ext)
            print('Sprk_10_first_ext are:', SPRk_ext_10first)
            print('eta_10_first_ext are:', eta_ext_10first)
            print('eta_10_first are:', eta_10first)
            print('++++++++++++++++++++++++++++++++++++++++++++++++')

        SPRtot_c_10first = np.zeros(10)
        Gamma_c_10first = np.zeros(10)
        delta_COB_sig_c_10first = np.zeros(10)
        eta_COB_c_10first = np.zeros(10)
        eta_c_10first = np.zeros(10)
        NSR1h_c_10first = np.zeros(10)
        w_c_key_10first = np.zeros(10)
        for l in range(nsprmax):
            m = sprk_sorted_index[l]
            Itc_acc_l = Itot - Ic[m]

            # Then we procedd to compute the secondary aperture
            NSR1h_c_l, w_c_l = aperture(ft=Ic[m], fc=Itc_acc_l, sb=sb, sd=sd, sq=sq,window_size=win_size,rho=rho_mask)
            w_c_size_l = w_c_l.sum()  # mask size
            NSR1h_c_10first[l] = NSR1h_c_l
            if(win_size<=8): w_c_key_10first[l] = mask_to_bitmask(w_c_l)
            SPRtot_c_10first[l] = np.sum(Itc_acc_l * w_c_l) / np.sum(Itot * w_c_l)
            eta_c_10first[l] = np.sqrt(td * ntr) * dback / NSR1h_c_l
            delta_COB_c_l, delta_COB_sig_c_l, Gamma_c_l = cob_shift(Itot, Ic[m], w_c_l, dback)
            if(w_c_size_l==1): delta_COB_c_l = 0.
            Gamma_c_10first[l] = Gamma_c_l
            delta_COB_sig_c_10first[l] = delta_COB_sig_c_l
            eta_COB_c_10first[l] =  delta_COB_c_l * np.sqrt(td * ntr) / delta_COB_sig_c_l

        save_info = np.array([ID_t, P_t, psf_idx, n_c, w_t_key, w_t_size, NSR1h, n_bad, SPR_crit, SPR_tot, ID_c,
                              m_c, sprk[ind_sprk], eta_t, delta_obs_t, delta_COB, delta_COB_sig_1h_24c,
                              eta_cob, d_c, Gamma, 0., n_bad_wrong])
        save_info = np.append(save_info, SPRk_10first)
        save_info = np.append(save_info, Gamma_10first)
        save_info = np.append(save_info, delta_COB_sig_10first)
        save_info = np.append(save_info, eta_COB_10first)
        save_info = np.append(save_info, eta_10first)
        save_info = np.append(save_info, IDs_10first)
        save_info = np.append(save_info, Delta_x_10first)
        save_info = np.append(save_info, Delta_y_10first)
        save_info = np.append(save_info, x_t_im)
        save_info = np.append(save_info, y_t_im)
        save_info = np.append(save_info, psf_idx_c)
        save_info = np.append(save_info, NSR1h_0)

        save_info_2ndmask = np.array([ID_t, P_t, ID_c, m_c, NSR1h_c, w_c_key, w_c_size, eta_c, delta_obs_c, delta_COB_c,
                                      delta_COB_sig_1h_24c_c, eta_cob_c, spr_tot_c, Gamma_c])
        save_info_2ndmask = np.append(save_info_2ndmask, SPRtot_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, Gamma_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, delta_COB_sig_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, NSR1h_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, w_c_key_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, eta_COB_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, eta_c_10first)
        save_info_2ndmask = np.append(save_info_2ndmask, NSR1h_c_0)



        save_info_ext = np.array([ID_t, P_t, w_ext_key, w_ext_size, NSR_ext_1h, n_bad_ext, SPR_crit_ext, SPR_tot_ext,
                                  ID_c, m_c, eta_ext, delta_obs_ext, delta_COB_ext, delta_COB_sig_1h_24c_ext,
                                  eta_cob_ext, delta_obs_ext, eta_ext_max, delta_obs_ext_max, n_det_ext, Gamma_ext, NSR_chi])
        save_info_ext = np.append(save_info_ext, SPRk_ext_10first)
        save_info_ext = np.append(save_info_ext, Gamma_ext_10first)
        save_info_ext = np.append(save_info_ext, delta_COB_sig_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_COB_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_ext_10first)
        save_info_ext = np.append(save_info_ext, eta_dtd_10first)

        save_info_bray = np.array([ID_t, P_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray])
        
        file_out.write('%.2f %.2f %i %.2f\n'% (ID_t, P_t, n_bad, sprk[ind_sprk]))

        return save_info, save_info_2ndmask, save_info_ext, save_info_bray


    def process_target_wrapper(k):
        try:
            result = process_target(k)
        except:
            result = None
            traceback.print_exc()
            print('an unkown error occur while processing star ID=%i (index=%i, magnitude interval %f-%f)'
                  % (ID_target[k], k, Pi - binsize / 2., Pi + binsize / 2.))
        return result


    if mp:  # multiprocessing mode
        nproc = min(n_t, nproc_max)
        pool = multiprocessing.Pool(nproc)
        results = pool.map(process_target_wrapper, [(t) for t in range(n_t)])
        pool.close()
        for k in range(n_t):
            if results[k] is not None:
                save_info[counter] = results[k][0]
                save_info_2ndmask[counter] = results[k][1]
                save_info_ext[counter] = results[k][2]
                save_info_bray[counter] = results[k][3]
                counter += 1
        #file_out.close()
    else:  # sequential mode
        for k in range(n_t):
            result = process_target(k)
            if result is not None:
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
file_out.close()
# ftrack.close()
save_info = save_info[0:counter]
save_info_ext = save_info_ext[0:counter]
save_info_bray = save_info_bray[0:counter]
# Now it is time to save the metrics of interest into a.npy file
np.save(DIRout + 'targets_P5.npy', save_info)
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
# 52-61: 10 first eta_COB
# 62-71: 10 first eta_10first
# 72-81: IDs of the 10 first contaminants
# 82-91: delta x
# 92-101: delta y
# 102: x target position in the imagette
# 103: Y target position in the imagette
# 104:  PSF index for the contaminant stars
# 105: theoretical n-mask NSR (calculation based on the target PSF assuming the same PSF for the contaminants)

np.save(DIRout + 'targets_P5_2ndmask.npy', save_info_2ndmask)
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
# 64-73: 10 first eta_COB_c
# 74-83: 10 first eta_c_10first
# 84: theoretical s-mask NSR (calculation based on the target PSF assuming the same PSF for the contaminants)


np.save(DIRout + 'targets_P5_extended.npy', save_info_ext)
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
# 20: NSRchi
# 21-30: 10 first SPRk values
# 31-40: 10 first Gamma values
# 41-50: 10 first delta_COB_sig_1h_24c
# 51-60: 10 first eta_COB_ext
# 61-70: 10 first eta_ext_10first
# 71-80: 10 first eta_dtd  (significance of the differental transit depth

np.save(DIRout + 'targets_P5_bray.npy', save_info_bray)

fout = open(DIRout + 'star_count.txt', 'w')
for i in range(nP):
    Pi = Pmin + i * binsize
    fout.write('%.2f %i\n' % (Pi, n_star_p_bin[i]))
fout.close()
