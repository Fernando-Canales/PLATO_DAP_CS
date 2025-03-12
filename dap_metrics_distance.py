"""
Script that computes the efficiency in the
distance space

Fernando 05.03.2025
"""
import numpy as np #type:ignore
import spline2dbase #type:ignore
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant
from imagette import window, centroid_shift
from NSR import spr_crit, aperture_computation, SPR, mask_to_bitmask, extended_binary_mask
from tqdm import tqdm # type:ignore

#CONFIGURATION PARAMETERS
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_4_22_hr_CONDITION_1-PIXEL_SEC_MASK/' 
PSFfile = 'PSF_Focus_0mu_0.2pxdif.npz'
DIRout =  '/home/fercho/double-aperture-photometry/simulation_results/rings/distances/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues
data_nominal_mask = np.load(dataDIR + 'targets_P5.npy')
# Parameters for the imagette and PSF decomposition
size_im_x = 6  # size of the imagette (x-direction)
size_im_y = 6  # size of the imagette (y-direction)
subres = 128   # resolution of the PSF
bsres = 20     # resolution of the b-spline decomposition of the PSF
# Parameters for the NSR
sb = (45. * 21) # Background noise from zodiacal light in e-/px (poisson noise)times integration time (21 sec.)
sd = 50.2      # Overall detector noise(includ. readout at beginning of life,smearing and dark current)in units of e-rms/px
sq = 7.2       # Quantization noise in units of e-rms/px

# Define ring numbers for x-axis (1 to number_of_circles)
number_of_circles = 7 
ring_numbers = np.arange(1, number_of_circles + 1)  # [1, 2, ..., number_of_circles]

# Define parameters
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other possible value: 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1

# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
transit_depth = 132000  # transit depth in ppm
transit_duration = 6.72*0.46**2   # transit duration in hours
ntr = 3        # number of transits in one hour
distance_max = 7 # maximum distance in pixels, from the target, to a star in the window in order to be considered a contaminant
#n_c_max = 300 # maximum number of contaminants in each window

Delta_P_max = 15.
data = np.load(cataDIR + 'SFP_DR3_20230101.npy') # star catalogue from GAIA
psfdata = np.load(PSFfile)                       # processed PSFs 
del_back, tr_dur = np.loadtxt(cataDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt', unpack=True, usecols=[0, 1]) # transit depth and duration from Kepler Eclipsing Binary Catalogue
number_of_circles = 7   # Number of circles                          # binsize around every magn

# First values
ID_cts = np.arange(0, data.shape[0]) # ID for every star in the catalogue
ID = data_nominal_mask[:, 0]

x_star = data[:, 3]              # x-coordinate on the focal plane for every star in the catalogue
y_star = data[:, 4]              # y-coordinate on the focal plane for every star in the catalogue
x_tar_fp_pix = data_nominal_mask[:, 219]
y_tar_fp_pix = data_nominal_mask[:, 220]
x_tar_mm, y_tar_mm = from_pix_2_mm(x_tar_fp_pix, y_tar_fp_pix)
distance_from_target_to_10first_contaminants = data_nominal_mask[:, 179:189] # within the imagette, in pixels
# Compute a representative distance for each target (e.g. the minimum distance among the 10)
rep_distance = np.median(distance_from_target_to_10first_contaminants, axis=1)
distance_from_target_to_10first_contaminants_1D = distance_from_target_to_10first_contaminants.flatten()

psfbs = psfdata['psfbs']         # PSF b-spline decomposition
pxc = psfdata['pxc']             # x-coordinate of the PSF b-spline decomposition centroid
pyc = psfdata['pyc']             # y-cooridnate of the PSF b-spline decomposition centroid
xpsf_pix = psfdata['xpsf_pix']   # x-coordinate of the PSF in pixel
ypsf_pix = psfdata['ypsf_pix']   # y-coordinate of the PSF in pixel

# Now we use a seed for obtaining the same number of targets and dback and td values
np.random.seed(300)
# Lists to store results for each ring
ring_results_nominal = []
ring_results_secondary = []
ring_results_extended = []
ring_results_bray = []

def process_target(target_index):
    ID_t = ID[target_index]
    m_t = data_nominal_mask[:, 1][target_index]     
    alpha = np.arctan(np.sqrt(x_tar_mm[target_index] ** 2 + y_tar_mm[target_index] ** 2) / 247.732)

    # Compute distance to PSF and find closest PSF
    s_d = (xpsf_pix - x_tar_fp_pix[target_index]) ** 2 + (ypsf_pix - y_tar_fp_pix[target_index]) ** 2
    psf_idx = np.argmin(s_d)
    x_t_im, y_t_im, i0, j0 = window(x_tar_fp_pix[target_index], y_tar_fp_pix[target_index], size_im_x, size_im_y)
    offx, offy = x_t_im - pxc[psf_idx], y_t_im - pyc[psf_idx]

    # Imagette and intensity calculations
    imagette = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y, offx=offx, offy=offy)
    Delta_P = data[:, 2] - m_t
    f_ref_t = reference_flux_target(m_t) * (np.cos(alpha) ** 2)
    It = f_ref_t * imagette

    # Now we find all the contaminants surrounding each TARGET. We put a distance condition (10 pixels)
    dist = np.sqrt((x_star - x_tar_fp_pix[target_index]) ** 2 + (y_star - y_tar_fp_pix[target_index]) ** 2)
    m = (dist > 0) & (dist < distance_max) & (Delta_P < Delta_P_max)  & (data[:, 2] > 0)# mask containing the distance condition for contaminants
    print(len(m))
    n = np.where(m)[0]             # index of the contaminants
    m_c = data[:, 2][n]            # magnitude of each contaminant star for a given TARGET
    x_c = x_star[n]                # x-coordinate of a given contaminant in the focal plane
    y_c = y_star[n]                # y-coordinate of a given contaminant in the focal plane
    n_c = len(x_c)                 # number of contaminants for a given TARGET
    ID_contaminants = ID_cts[m] # IDs of the contaminant stars
    x_c_im = x_c - i0              # x-coordinate of a given contaminant inside the imagette (window)
    y_c_im = y_c - j0              # y-coordinate of a given contaminant inside the imagette (window)
    offx_c = x_c_im - pxc[psf_idx] # x-coordinate of the offset between the center of each contaminant and each PSF
    offy_c = y_c_im - pyc[psf_idx] # y-coordinate of the offse between the center of each contaminant and each PSF

        # ---------------------------------------------------------------------
    # Calcular la distancia entre el target y cada contaminante (en píxeles)
    # ---------------------------------------------------------------------
    # Filtrar contaminantes: dist > 0, dist < distance_max, y condición adicional (magnitud positiva)
    mask_cont = (dist > 0) & (dist < distance_max) & (data[:, 2] > 0)
    idx_cont = np.where(mask_cont)[0]
    n_c = len(idx_cont)
    if n_c == 0:
        return np.empty((0, 4))
    
    # Obtaining magnitude_difference
    mag_cont = data[idx_cont, 2]
    delta_mag = m_t - mag_cont
    
    # ---------------------------------------------------------------------


    # Now we have to build an 'imagette' for every contaminant
    Ic = np.zeros((n_c, 6, 6))
    for o in range(1, n_c + 1):
    # Now we make sure to deal only with stars with positive magnitudes
        if m_c[o - 1] > 0:
        # Then we compute the imagette for each contaminant by integrating the b-spline decomposition of the PSF
            Ic[o - 1, :, :] = spline2dbase.Spline2Imagette(psfbs[psf_idx], bsres, size_im_x, size_im_y, offx=offx_c[o - 1], offy=offy_c[o - 1])
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
    nsr_1h_24_cameras_nominal_mask = ((10 ** 6) / (12 * np.sqrt(24))) * nsr_nominal_mask
    nsr_1h_6_cameras_nominal_mask = ((10 ** 6) / (12 * np.sqrt(6))) * nsr_nominal_mask

    # sprk, sprk_max, SPR_tot, n_bad = SPR(SPR_crit=SPR_crit, n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=w_t)
    sprk, SPR_tot = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=nominal_mask) # sprk for every contaminant (as well as sprkmax and spr-tot)
    sprk_6_cameras,  SPR_tot_6_cameras = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=nominal_mask)
    # We give to each contaminant a transit depth and duration
    dback = np.zeros(n_c)
    td = np.zeros(n_c)

    # Define whether to use fixed values or random values from the catalogue
    use_fixed_values = True  # Change to False if you want to use random values from the catalogue
    if use_fixed_values:
    # Case 1: Use fixed values
        td.fill(transit_duration)  # Fill the array with the fixed value
        dback.fill(transit_depth)  # Fill the array with the fixed value
    else:
    # We obtain 10 random values for transit depth and transit duration
        for n in range(n_c):
            dback_index = np.random.randint(0, len(del_back))
        #td_index = np.random.randint(0, len(tr_dur))

            dback[n] = del_back[dback_index]
            td[n] = tr_dur[dback_index]

    # We compute the critical SPR now
    #SPR_crit = spr_crit(dback=dback, SPR_tot=SPR_tot, nsr=NSR1h, td=td, ntr=ntr) # spr_crit w.r.t. the nominal mask
    SPR_crit_24_cameras = spr_crit(dback=dback[0], SPR_tot=SPR_tot, nsr=nsr_1h_24_cameras_nominal_mask, td=td[0], ntr=ntr) # spr_crit w.r.t. the nominal mask
    SPR_crit_6_cameras = spr_crit(dback=dback[0], SPR_tot=SPR_tot, nsr=nsr_1h_6_cameras_nominal_mask, td=td[0], ntr=ntr)

    n_bad = np.sum(sprk > SPR_crit_24_cameras) # N_bad
    index_contaminant_highest_sprk = np.argmax(sprk)      # index of the contaminant with the highest sprk w.r.t. nominal mask
    m_c_bad = m_c[index_contaminant_highest_sprk]         # magnitude for the contaminant with the highest sprk value w.r.t the nominal mask
    ID_contaminant_highest_sprk = ID_contaminants[index_contaminant_highest_sprk]

    # Now we select the (intensity per pixel array) 'imagette' of the contaminant star with the highest sprk value
    # w.r.t. the nominal mask
    Ic_max = Ic[index_contaminant_highest_sprk]
    dback_contaminant_highest_spr = dback[index_contaminant_highest_sprk]
    td_contaminant_highest_spr = td[index_contaminant_highest_sprk]

    # Now we find the distance between the target and the contaminant star with the highest value of sprk w.r.t the
    # nominal mask
    dist_bad = np.sqrt((x_t_im - x_c_im[index_contaminant_highest_sprk]) ** 2 + (y_t_im - y_c_im[index_contaminant_highest_sprk]) ** 2)
    # Now we define a term that contains the flux of the targets and their respective contaminants except the
    # one with the highest value of SPRk (i.e. the contaminant of interest)
    Itc_acc = It + Ic_acc - Ic_max

    #NSR1h_c, NSR1h_c_6, w_c, w_c_6 = aperture(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq) # secondary mask for the most problematic contaminant
    secondary_mask = aperture_computation(ft=Ic_max, fc=Itc_acc, sb=sb, sd=sd, sq=sq)
    secondary_mask_key = mask_to_bitmask(secondary_mask)
    secondary_mask_size = np.count_nonzero(secondary_mask)

    nsr_secondary_mask = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * (secondary_mask**2))) / np.sum(Ic_max * secondary_mask) # E1. 11 in Marchiori paper
    nsr_1h_24_cameras_secondary_mask = ((10 ** 6) / (12 * np.sqrt(24))) * nsr_secondary_mask 
    nsr_1h_6_cameras_secondary_mask = ((10 ** 6) / (12 * np.sqrt(6))) * nsr_secondary_mask 

    sprk_sec, SPR_tot_sec = SPR(n_c=n_c, f_contaminant=Ic, f_tot=f_tot, w=secondary_mask) # sprk for every contaminant (as well as sprkmax and spr-tot) w.r.t. secondary mask
    # We compute the flux over the secondary mask and the extended secondary mask
    f_beb = Ic_acc * secondary_mask
    f_t_c = It * secondary_mask
    f_tot_c = np.sum(f_t_c + f_beb)
    # We compute spr_tot_c, that is, the expression given in Marchiori presentation for PLATO week #8
    spr_tot_secondary_mask = np.sum(Itc_acc * secondary_mask) / f_tot_c

    # We compute now the delta_obs and etas for the two apertures w.r.t. the most significant contaminant star
    delta_obs_nominal_mask = sprk[index_contaminant_highest_sprk] * dback_contaminant_highest_spr                                         # observed transit depth in the nominal mask
    delta_obs_secondary_mask = (1 - spr_tot_secondary_mask) * dback_contaminant_highest_spr                                        # observed transit depth in the secondary mask
    eta_t = sprk[index_contaminant_highest_sprk] * np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / (nsr_1h_24_cameras_nominal_mask * (1 - SPR_tot)) # signal statistical significance in the nominal mask
    eta_t_6_cameras = sprk[index_contaminant_highest_sprk] * np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr/ (nsr_1h_6_cameras_nominal_mask * (1 - SPR_tot)) # signal statistical significance in the nominal mask
    eta_c = np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / nsr_1h_24_cameras_secondary_mask # signal statistical significance in the secondary mask
    eta_c_6_cameras = np.sqrt(td_contaminant_highest_spr * ntr) * dback_contaminant_highest_spr / nsr_1h_6_cameras_secondary_mask

    eta_true_positive_24_cameras_earth_like = 84 * np.sqrt(13 * 3) / nsr_1h_24_cameras_nominal_mask
    eta_true_positive_6_cameras_jovian_planet =  10100 * np.sqrt(29.6 * 3) / nsr_1h_6_cameras_nominal_mask
    eta_true_positive_24_cameras_super_earth = 522 * np.sqrt(42 * 3) / nsr_1h_24_cameras_nominal_mask
    eta_true_positive_6_cameras_super_earth = 522 * np.sqrt(42 * 3) / nsr_1h_6_cameras_nominal_mask


    nsprmax = min(10, n_c)
    eta_10first = np.zeros(10)
    sprk_10first = np.zeros(10)
    eta_cob_10first = np.zeros(10)
    sigma_cob_10first = np.zeros(10)
    gamma_nom_10first = np.zeros(10)
    abs_cob_shift_10first = np.zeros(10)
    eta_cob_10first_6_cameras = np.zeros(10)
    gamma_nom_10first_6cameras = np.zeros(10)
    sigma_cob_10first_6_cameras = np.zeros(10)
    abs_cob_shift_10first_6_cameras = np.zeros(10)
    IDs_from_the_10first_contaminants = np.zeros(10)
    dist_from_target_to_10first_contaminants = np.zeros(10)
    delta_x_from_target_to_10first_contaminants = np.zeros(10)
    delta_y_from_target_to_10first_contaminants = np.zeros(10)
    dback_10first = np.zeros(10)
    td_10first = np.zeros(10)

    magnitude_contaminant_10first = np.zeros(10)
    x_coordinate_in_the_imagette_for_a_contaminant_10first = np.zeros(10)
    y_coordinate_in_the_imagette_for_a_contaminant_10first = np.zeros(10)
    #dback_10first_index = np.random.randint(0, len(del_back), size=10)
    # sorting the SPRk by decreasing order and taking the 10 first values
    sprk_sorted_index = (np.argsort(sprk)[::-1])
        # Ordenar contaminantes por SPR (descendente) y tomar los 10 primeros
    sprk_sorted_index = (np.argsort(sprk)[::-1])
    top_n = min(10, n_c)
    top_indices = idx_cont[sprk_sorted_index[:top_n]]
    
    # Extraer distancias y delta_mag de los top contaminantes
    dist_top = dist[top_indices]
    delta_mag_top = m_t - data[top_indices, 2]
    
    # Formar el array de pares con la información básica: 
    # [distance, ID_target, ID_contaminant, delta_mag]
    pairs = np.column_stack((dist_top,
                             np.full(top_n, ID_t),
                             top_indices, delta_mag_top))
    td_10first[0:nsprmax] = td[sprk_sorted_index[0:nsprmax]]
    sprk_10first[0:nsprmax] = sprk[sprk_sorted_index[0:nsprmax]]
    dback_10first[0:nsprmax] = dback[sprk_sorted_index[0:nsprmax]]
    IDs_from_the_10first_contaminants[0:nsprmax] = ID_contaminants[sprk_sorted_index[0:nsprmax]]
    x_coordinate_in_the_imagette_for_a_contaminant_10first[0:nsprmax] = x_c_im[sprk_sorted_index[0:nsprmax]]
    y_coordinate_in_the_imagette_for_a_contaminant_10first[0:nsprmax] = y_c_im[sprk_sorted_index[0:nsprmax]]
    delta_x_from_target_to_10first_contaminants[0:nsprmax] = x_c_im[sprk_sorted_index[0:nsprmax]] - x_t_im 
    delta_y_from_target_to_10first_contaminants[0:nsprmax] = y_c_im[sprk_sorted_index[0:nsprmax]] - y_t_im 
    dist_from_target_to_10first_contaminants[0:nsprmax] = np.sqrt((x_t_im - x_c_im[sprk_sorted_index[0:nsprmax]]) ** 2 + (y_t_im - y_c_im[sprk_sorted_index[0:nsprmax]]) ** 2)
    magnitude_contaminant_10first[0:nsprmax] = m_c[sprk_sorted_index[0:nsprmax]]
    for l in range(nsprmax):
        m = sprk_sorted_index[l]
        eta_10first[l] = sprk[m] * np.sqrt(td[m] * ntr) * dback[m] / nsr_1h_24_cameras_nominal_mask / (1. - SPR_tot)
        eta_cob_10first[l], sigma_cob_10first[l], abs_cob_shift_10first[l], gamma_nom_10first[l] = centroid_shift(w=nominal_mask, Ik=Ic[m],
        n_cam=24, I_t=It, I_contaminants=Ic_acc, sprk=sprk[m], dback=dback_10first[l], sb=sb, sd=sd, sq=sq, td=td_10first[l], ntr=ntr)
        eta_cob_10first_6_cameras[l], sigma_cob_10first_6_cameras[l], abs_cob_shift_10first_6_cameras[l], gamma_nom_10first_6cameras[l] = centroid_shift(w=nominal_mask, Ik=Ic[m],
        n_cam=6, I_t=It, I_contaminants=Ic_acc, sprk=sprk[m], dback=dback_10first[l], sb=sb, sd=sd, sq=sq, td=td_10first[l], ntr=ntr)
    # -------------------------------------------NOMINAL COB-------------------------------------------------------#
    ## 24 cameras
    eta_cob, sigma_1_24, abs_cob, gamma_nom = centroid_shift(w=nominal_mask, Ik=Ic_max, n_cam=24, I_t=It, I_contaminants=Ic_acc, 
                                    sprk=sprk[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)

    ## 6 cameras 
    eta_cob_6_cameras, sigma_1_6_cameras, abs_cob_6_cameras, gamma_nom_6cameras = centroid_shift(w=nominal_mask, Ik=Ic_max, n_cam=6, I_t=It,
    I_contaminants=Ic_acc, sprk=sprk[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr) 
    # -------------------------------------------NOMINAL COB-------------------------------------------------------#

    # ------------------------------------------SECONDARY COB------------------------------------------------------#
    ## 24 cameras
    eta_cob_c, sigma_1_24_c, abs_cob_c, gamma_cob_c = centroid_shift(w=secondary_mask, Ik=Ic_max, n_cam=24, I_t=It, I_contaminants=Ic_acc, 
                                    sprk=sprk_sec[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)

    sprk_cob_c = np.sum(Ic_max * secondary_mask)/ np.sum((It + Ic_acc) * secondary_mask)
    ## 6 cameras
    eta_cob_c_6_cameras, sigma_1_6_cameras_c, abs_cob_c_6_cameras, gamma_cob_c_6_cameras = centroid_shift(w=secondary_mask, Ik=Ic_max, n_cam=6, I_t=It, I_contaminants=Ic_acc, 
                                    sprk=sprk_sec[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
    # ------------------------------------------SECONDARY COB------------------------------------------------------#

    ################################################################################################################
    #                                               EXTENDED MASK                                                  #
    ################################################################################################################
    extended_mask = extended_binary_mask(nominal_mask, W=1) # extended mask (1 pixel ring around the nominal mask)
    extended_mask_key = mask_to_bitmask(extended_mask)     # mask key for the extended mask
    extended_mask_size = np.count_nonzero(extended_mask)   # mask size for the extended mask

    extended_mask_2_pix = extended_binary_mask(nominal_mask, W=2) # extended mask (1 pixel ring around the nominal mask)
    extended_mask_key_2_pix = mask_to_bitmask(extended_mask_2_pix)     # mask key for the extended mask
    extended_mask_size_2_pix = np.count_nonzero(extended_mask_2_pix)   # mask size for the extended mask

    # Now we compute all the metrics associated with every extended mask.
    NSR_ext = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * extended_mask)) / np.sum(It * extended_mask) # E1. 11 in Marchiori paper
    NSR_ext_1h_24_cameras = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext
    NSR_ext_1h_6_cameras = ((10 ** 6) / (12 * np.sqrt(6))) * NSR_ext
    sprk_ext, SPR_tot_ext = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=extended_mask)
    SPR_crit_ext = spr_crit(dback=dback_10first[0], SPR_tot=SPR_tot_ext, nsr=NSR_ext_1h_24_cameras, td=td_10first[0], ntr=ntr)
    n_bad_ext = np.sum(sprk_ext > SPR_crit_ext)                           # N_bad for the extended mask
    delta_obs_ext = sprk_ext[index_contaminant_highest_sprk] * dback_10first[0]                            # w.r.t. the nominal mask
    eta_ext = sprk_ext[index_contaminant_highest_sprk] * np.sqrt(td_10first[0] * ntr) * dback_10first[0] / (NSR_ext_1h_24_cameras *(1 - SPR_tot_ext)) # w.r.t the nominal mask
    #eta_ext_6_cameras = sprk_ext[index_contaminant_highest_sprk] * np.sqrt(td * ntr) * dback / (NSR_ext_1h_6_cameras *(1 - SPR_tot_ext))
    NSR_ext_2_pix = np.sqrt(np.sum((It + Ic_acc + sb + sd ** 2 + sq ** 2) * extended_mask_2_pix)) / np.sum(It * extended_mask_2_pix) # E1. 11 in Marchiori paper
    NSR_ext_2_pix_1h_24_cameras = ((10 ** 6) / (12 * np.sqrt(24))) * NSR_ext_2_pix
    NSR_ext_2_pix_1h_6_cameras = ((10 ** 6) / (12 * np.sqrt(6))) * NSR_ext_2_pix
    sprk_ext_2_pix, SPR_tot_ext_2_pix = SPR(n_c=n_c, f_contaminant=Ic, f_tot=(It + Ic_acc), w=extended_mask_2_pix)
    SPR_crit_ext_2_pix = spr_crit(dback=dback_10first[0], SPR_tot=SPR_tot_ext_2_pix, nsr=NSR_ext_2_pix_1h_24_cameras, td=td_10first[0], ntr=ntr)
    n_bad_ext_2_pix = np.sum(sprk_ext_2_pix > SPR_crit_ext_2_pix)                           # N_bad for the extended mask
    delta_obs_ext_2_pix = sprk_ext_2_pix[index_contaminant_highest_sprk] * dback_10first[0]                            # w.r.t. the nominal mask
    eta_ext_2_pix = sprk_ext_2_pix[index_contaminant_highest_sprk] * np.sqrt(td_10first[0] * ntr) * dback_10first[0] / (NSR_ext_2_pix_1h_24_cameras *(1 - SPR_tot_ext_2_pix)) # w.r.t the nominal mask
    # ------------------------------------------EXTENDED COB-------------------------------------------------------#
    ## 24 cameras w.r.t. the most significant contaminant
    eta_cob_ext, sigma_1_24_ext, abs_cob_ext, gamma_ext = centroid_shift(w=extended_mask, Ik=Ic_max, n_cam=24, 
    I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
    eta_cob_ext_2_pix, sigma_1_24_ext_2_pix, abs_cob_ext_2_pix, gamma_ext_2_pix = centroid_shift(w=extended_mask_2_pix, Ik=Ic_max, n_cam=24, 
    I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext_2_pix[index_contaminant_highest_sprk], dback=dback_10first[0], sb=sb, sd=sd, sq=sq, td=td_10first[0], ntr=ntr)
    # ------------------------------------------EXTENDED COB-------------------------------------------------------#

    # Now we will take the first 10 contaminants in terms of SPR

    # We take the 10 first contaminants in terms of SPRk to use them for the extended mask
    nsprmax = min(10, n_c)
    sprk_sorted_index = (np.argsort(sprk)[::-1])

    SPRK_ext_10first = np.zeros(10)
    eta_ext_10first = np.zeros(10)
    gamma_ext_10first = np.zeros(10)
    eta_ext_10first_6_cameras = np.zeros(10)
    eta_cob_ext_10first = np.zeros(10)
    eta_cob_ext_10first_6_cameras = np.zeros(10)
    sigma_cob_ext_10first = np.zeros(10)
    sigma_cob_ext_10first_6_cameras = np.zeros(10)
    abs_cob_shift_ext_10first = np.zeros(10)
    abs_cob_shift_ext_10first_6_cameras = np.zeros(10)
    gamma_ext_10first_6_cameras = np.zeros(10)
    dback_ext_10first = np.zeros(10)
    td_ext_10first = np.zeros(10)

    SPRK_ext_2_pix_10first = np.zeros(10)
    eta_ext_2_pix_10first = np.zeros(10)
    gamma_ext_2_pix_10first = np.zeros(10)
    eta_ext_2_pix_10first_6_cameras = np.zeros(10)
    eta_cob_ext_2_pix_10first = np.zeros(10)
    eta_cob_ext_2_pix_10first_6_cameras = np.zeros(10)
    sigma_cob_ext_2_pix_10first = np.zeros(10)
    sigma_cob_ext_2_pix_10first_6_cameras = np.zeros(10)
    abs_cob_shift_ext_2_pix_10first = np.zeros(10)
    abs_cob_shift_ext_2_pix_10first_6_cameras = np.zeros(10)
    gamma_ext_2_pix_10first_6_cameras = np.zeros(10)


    SPRK_ext_10first[0:nsprmax] = sprk_ext[sprk_sorted_index[0:nsprmax]]
    dback_ext_10first[0:nsprmax] = dback[sprk_sorted_index[0:nsprmax]]
    td_ext_10first[0:nsprmax] = td[sprk_sorted_index[0:nsprmax]]
    SPRK_ext_2_pix_10first[0:nsprmax] = sprk_ext_2_pix[sprk_sorted_index[0:nsprmax]]
    for l in range(nsprmax):
        m = sprk_sorted_index[l]
        eta_ext_10first[l] = sprk_ext[m] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_1h_24_cameras * (1 - SPR_tot_ext))       
        eta_ext_10first_6_cameras[l] = sprk_ext[l] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_1h_6_cameras * (1 - SPR_tot_ext))
        eta_cob_ext_10first[l], sigma_cob_ext_10first[l], abs_cob_shift_ext_10first[l], gamma_ext_10first[l] = centroid_shift(w=extended_mask, Ik=Ic[m],  
        n_cam=24, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)
        eta_cob_ext_10first_6_cameras[l], sigma_cob_ext_10first_6_cameras[l], abs_cob_shift_ext_10first_6_cameras[l], gamma_ext_10first_6_cameras[l] = centroid_shift(w=extended_mask, Ik=Ic[m],
        n_cam=6, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)
        eta_ext_2_pix_10first[l] = sprk_ext_2_pix[m] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_2_pix_1h_24_cameras * (1 - SPR_tot_ext_2_pix))       
        eta_ext_2_pix_10first_6_cameras[l] = sprk_ext_2_pix[l] * np.sqrt(td_ext_10first[l] * ntr) * dback_ext_10first[l] / (NSR_ext_2_pix_1h_6_cameras * (1 - SPR_tot_ext_2_pix))
        eta_cob_ext_2_pix_10first[l], sigma_cob_ext_2_pix_10first[l], abs_cob_shift_ext_2_pix_10first[l], gamma_ext_2_pix_10first[l] = centroid_shift(w=extended_mask_2_pix, Ik=Ic[m],  
        n_cam=24, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)
        eta_cob_ext_2_pix_10first_6_cameras[l], sigma_cob_ext_2_pix_10first_6_cameras[l], abs_cob_shift_ext_2_pix_10first_6_cameras[l], gamma_ext_2_pix_10first_6_cameras[l] = centroid_shift(w=extended_mask_2_pix, Ik=Ic[m],
        n_cam=6, I_t=It, I_contaminants=Ic_acc, sprk=sprk_ext_2_pix[m], dback=dback_ext_10first[l], sb=sb, sd=sd, sq=sq, td=td_ext_10first[l], ntr=ntr)       
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
    # n_eff_ext = len(np.where((sprk_ext > SPR_crit_ext) & (sprk[index_contaminant_highest_sprk] > SPR_crit) & (sprk_ext > sprk[index_contaminant_highest_sprk])
    # )[0])
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('P magnitude of the target=', m_t)
    print('Total number of contaminant stars for this target = ', n_c)
    print('Nbad for this target = ', n_bad)
    print('Most significant contaminant ID =', ID_contaminant_highest_sprk)
    print('Transit depth of the contaminant with the highest spr =', dback_contaminant_highest_spr, 'ppm')
    print('Magnitude difference between the target and the contaminant with the highest sprk = ',
        m_c[index_contaminant_highest_sprk] - m_t)
    print('Distance between the target and the contaminant with the highest sprk = ',
        (x_t_im - x_c_im[index_contaminant_highest_sprk]) ** 2 + (y_t_im - y_c_im[index_contaminant_highest_sprk]) ** 2)
    print('IDs from the 10 first contaminants:', IDs_from_the_10first_contaminants)
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

    # Now we save the important metrics for every target w.r.t the nominal mask
    save_info = np.array([ID_t, m_t, n_c, m_c_bad, nominal_mask_key, nominal_mask_size, nsr_1h_24_cameras_nominal_mask, n_bad, SPR_crit_24_cameras, sprk[index_contaminant_highest_sprk], SPR_tot, eta_t, delta_obs_nominal_mask, abs_cob, eta_cob, 
                        sigma_1_24])

    save_info = np.append(save_info, sprk_10first)
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
    save_info = np.append(save_info, gamma_nom_10first)
    save_info = np.append(save_info, gamma_nom_10first_6cameras)
    save_info = np.append(save_info, td_10first)
    save_info = np.append(save_info, dback_10first)
    save_info = np.append(save_info, gamma_nom)
    save_info = np.append(save_info, gamma_nom_6cameras)
    save_info = np.append(save_info, nsr_1h_6_cameras_nominal_mask)
    save_info = np.append(save_info, x_coordinate_in_the_imagette_for_a_contaminant_10first)
    save_info = np.append(save_info, y_coordinate_in_the_imagette_for_a_contaminant_10first)
    save_info = np.append(save_info, magnitude_contaminant_10first)
    save_info = np.append(save_info, dist_from_target_to_10first_contaminants)
    save_info = np.append(save_info, IDs_from_the_10first_contaminants)
    save_info = np.append(save_info, delta_x_from_target_to_10first_contaminants)
    save_info = np.append(save_info, delta_y_from_target_to_10first_contaminants)


    # Now we save the important metrics w.r.t the secondary mask
    save_info_contaminant = np.array([ID_t, m_t, secondary_mask_key, secondary_mask_size, nsr_1h_24_cameras_secondary_mask, spr_tot_secondary_mask, eta_c, delta_obs_secondary_mask, abs_cob_c, eta_cob_c, sigma_1_24_c])
    save_info_contaminant = np.append(save_info_contaminant, eta_c_6_cameras)
    save_info_contaminant = np.append(save_info_contaminant, eta_cob_c_6_cameras)
    save_info_contaminant = np.append(save_info_contaminant, abs_cob_c_6_cameras)
    save_info_contaminant = np.append(save_info_contaminant, sigma_1_6_cameras_c)
    save_info_contaminant = np.append(save_info_contaminant, gamma_cob_c_6_cameras)
    save_info_contaminant = np.append(save_info_contaminant, SPR_tot_sec)
    save_info_contaminant = np.append(save_info_contaminant, gamma_cob_c)
    save_info_contaminant = np.append(save_info_contaminant, nsr_1h_6_cameras_secondary_mask)
    save_info_contaminant = np.append(save_info_contaminant, ID_contaminant_highest_sprk)

    # Now we save the important metrics w.r.t the extended mask
    save_info_ext = np.array([ID_t, m_t, extended_mask_key, extended_mask_size, NSR_ext_1h_24_cameras, sprk_ext[index_contaminant_highest_sprk], SPR_crit_ext, eta_ext, delta_obs_ext, abs_cob_ext, eta_cob_ext, 
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
    save_info_ext = np.append(save_info_ext, gamma_ext_10first)
    save_info_ext = np.append(save_info_ext, gamma_ext_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, gamma_ext)
    save_info_ext = np.append(save_info_ext, SPRK_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, eta_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, eta_ext_2_pix_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, NSR_ext_2_pix_1h_6_cameras)
    save_info_ext = np.append(save_info_ext, eta_cob_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, sigma_cob_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, abs_cob_shift_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, eta_cob_ext_2_pix_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, sigma_cob_ext_2_pix_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, abs_cob_shift_ext_2_pix_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, gamma_ext_2_pix_10first)
    save_info_ext = np.append(save_info_ext, gamma_ext_2_pix_10first_6_cameras)
    save_info_ext = np.append(save_info_ext, gamma_ext_2_pix)
    save_info_ext = np.append(save_info_ext, NSR_ext_2_pix_1h_24_cameras)
    save_info_ext = np.append(save_info_ext, extended_mask_key_2_pix)
    save_info_ext = np.append(save_info_ext, extended_mask_size_2_pix)
    save_info_ext = np.append(save_info_ext, delta_obs_ext_2_pix)
    save_info_ext = np.append(save_info_ext, eta_ext_2_pix)
    save_info_ext = np.append(save_info_ext, n_bad_ext_2_pix)
    save_info_ext = np.append(save_info_ext, eta_cob_ext_2_pix)
    save_info_ext = np.append(save_info_ext, sigma_1_24_ext_2_pix)
    save_info_ext = np.append(save_info_ext, abs_cob_ext_2_pix)
    save_info_ext = np.append(save_info_ext, SPR_crit_ext_2_pix)
    save_info_ext = np.append(save_info_ext, SPR_tot_ext_2_pix)

    # Now we save the import metrics w.r.t the 4 pixel mask described by Bray et al. 2023
    save_info_bray = np.array([ID_t, m_t, n_c, NSR_bray_1h, n_bad_bray, SPR_crit_bray, sprk_bray[index_contaminant_highest_sprk], SPR_tot_bray])
    return save_info, save_info_contaminant, save_info_ext, save_info_bray


# =============================================================================
# Global accumulation: process all targets and store outputs for each mask and pairs.
# =============================================================================
all_nominal = []
all_secondary = []
all_extended = []
all_bray = []

# Instead of radius_focal_plane, we use our fixed maximum distance for the contaminants:
max_distance = 7.0  # maximum distance in pixels
number_of_shells = 7
shell_edges = np.sqrt(np.linspace(0, max_distance**2, number_of_shells + 1))
print("Shell edges (pixels):", shell_edges)


for i in range(number_of_shells):
    r_lo = shell_edges[i]
    r_hi = shell_edges[i+1]
    # Select targets whose representative distance is in the current shell
    targets_in_shell = np.where((rep_distance >= r_lo) & (rep_distance < r_hi))[0]
    print(f"Shell {i}: Found {targets_in_shell.size} targets with rep. distance in [{r_lo:.2f}, {r_hi:.2f})")

    ring_nominal_mask_results = []
    ring_secondary_mask_results = []
    ring_extended_mask_results = []
    ring_bray_mask_results = []
    
    for t in tqdm(targets_in_shell, desc=f"Processing Shell {i+1}", unit="target"):
         results = process_target(t)
         # Save or accumulate these results for later analysis.

        # Append results to the current ring's result lists
         ring_nominal_mask_results.append(results[0])
         ring_secondary_mask_results.append(results[1])
         ring_extended_mask_results.append(results[2])
         ring_bray_mask_results.append(results[3])
        # counter = counter +1
        # Convert the lists to arrays and store them in the main list for each ring
    ring_results_nominal.append(np.array(ring_nominal_mask_results))
    ring_results_secondary.append(np.array(ring_secondary_mask_results))
    ring_results_extended.append(np.array(ring_extended_mask_results))
    ring_results_bray.append(np.array(ring_bray_mask_results))

    # Save each ring's results as separate .npy files
for i, (nominal, secondary, extended, bray) in enumerate(zip(ring_results_nominal, ring_results_secondary, ring_results_extended, ring_results_bray)):
    np.save(f"{DIRout}ring_{i}_nominal.npy", nominal)
    np.save(f"{DIRout}ring_{i}_secondary.npy", secondary)
    np.save(f"{DIRout}ring_{i}_extended.npy", extended)
    np.save(f"{DIRout}ring_{i}_bray.npy", bray)

print("Ring-based processing and saving complete.")