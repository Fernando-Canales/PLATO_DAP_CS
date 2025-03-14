"""
Script that computes the efficiency in the
distance space

Fernando 05.03.2025
"""
import numpy as np #type:ignore
import matplotlib.pyplot as plt # type: ignore
from fitting_psf import from_pix_2_mm, reference_flux_target, reference_flux_contaminant

from tqdm import tqdm # type:ignore

#CONFIGURATION PARAMETERS
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_4_22_hr_CONDITION_1-PIXEL_SEC_MASK/' 
PSFfile = 'PSF_Focus_0mu_0.2pxdif.npz'
DIRout =  '/home/fercho/double-aperture-photometry/simulation_results/rings/distances/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/' # directory with all star catalogues
data_nominal_mask = np.load(dataDIR + 'targets_P5.npy')
data_extended_mask = np.load(dataDIR + 'targets_P5_extended.npy')
data_secondary_mask = np.load(dataDIR + 'targets_P5_secondary.npy')
# Parameters for the imagette and PSF decomposition

# Define parameters
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other possible value: 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1

eta_ext_threshold = 3
# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
transit_depth = 132000  # transit depth in ppm
transit_duration = 6.72*0.46**2   # transit duration in hours
ntr = 3        # number of transits in one hour
distance_max = 7 # maximum distance in pixels, from the target, to a star in the window in order to be considered a contaminant
#n_c_max = 300 # maximum number of contaminants in each window                      # processed PSFs 
del_back, tr_dur = np.loadtxt(cataDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt', unpack=True, usecols=[0, 1]) # transit depth and duration from Kepler Eclipsing Binary Catalogue

distance_from_target_to_10first_contaminants = data_nominal_mask[:, 179:189] # within the imagette, in pixels
flat_distance = distance_from_target_to_10first_contaminants.flatten()


# Instead of radius_focal_plane, we use our fixed maximum distance for the contaminants:
max_distance = 7.0  # maximum distance in pixels
number_of_shells = 6
shell_edges = np.sqrt(np.linspace(0, max_distance**2, number_of_shells + 1))
print("Shell edges (pixels):", shell_edges)
td = data_nominal_mask[:, 126:136]
dback = data_nominal_mask[:, 136:146]


 #Precompute flattened arrays for the extended mask metrics outside the shell loop:
flat_Sprk_ext    = data_extended_mask[:, 14:24].flatten()  # shape: (n_targets*10,)
flat_dback_ext   = data_nominal_mask[:, 136:146].flatten()  # assuming dback is taken from the nominal mask
flat_val_ext     = data_extended_mask[:, 4].flatten()       # used in the denominator
flat_frac_ext    = (1 - data_extended_mask[:, 13]).flatten()  # also used in the denominator
# Compute the eta metrics for all targets (for 24 cameras)
#eta_nom = (dback * data_nominal_mask[:, 17:27] * np.sqrt(td * ntr) /
#           (data_nominal_mask[:, 7][:, None] * (1 - data_nominal_mask[:, 11][:, None])))

# Flatten the eta arrays so that each target’s 10 values become one long vector:
#flat_eta_nom = eta_nom.flatten()  # shape: (n_targets * 10,)
# For each shell, create a boolean mask using the flattened distances.
for i in range(number_of_shells):
    r_lo = shell_edges[i]
    r_hi = shell_edges[i+1]
    mask = (flat_distance >= r_lo) & (flat_distance < r_hi)
    print(f"Shell {i}: Found {np.sum(mask)} contaminants with distance in [{r_lo:.2f}, {r_hi:.2f})")


    # Now flatten the other arrays in the same way:
    flat_Sprk_ext = data_extended_mask[:, 14:24].flatten()  # shape: (n_targets*10,)
    flat_Sprk_nom = data_nominal_mask[:, 17:27].flatten()  # shape: (n_targets*10,)
    flat_dback = data_nominal_mask[:, 136:146].flatten()
    flat_td = data_nominal_mask[:, 16:136].flatten()
    flat_eta_cob_nom = data_nominal_mask[:, 46:56].flatten()
    flat_eta_cob_ext = data_extended_mask[:, 45:55].flatten()
    #flat_val_ext     = np.repeat(data_extended_mask[:, 4], 10)     # used in the denominator
    #flat_frac_ext    = np.repeat((1 - data_extended_mask[:, 13]), 10) # also used in the denominator
    flat_eta_secondary = np.repeat(data_secondary_mask[:, 6], 10)
    flat_delta_obs_sec = np.repeat(data_secondary_mask[:, 7], 10)
    flat_delta_obs_sec_target = np.repeat(data_secondary_mask[:, 13], 10)
    flat_td = data_nominal_mask[:, 126:136].flatten()
    flat_eta_nom = data_nominal_mask[:, 27:37].flatten()
    flat_eta_ext = data_extended_mask[:, 25:35].flatten()

    # Filter the flattened eta values using the mask:
    eta_nom_in_shell = flat_eta_nom[mask]
    eta_ext_in_shell = flat_eta_ext[mask]

    Sprk_ext_in_shell = flat_Sprk_ext[mask]
    Sprk_nom_in_shell = flat_Sprk_nom[mask]
    eta_cob_nom_in_shell = flat_eta_cob_nom[mask]
    eta_cob_ext_in_shell = flat_eta_cob_ext[mask]
    eta_sec_in_shell = flat_eta_secondary[mask]
    dback = flat_dback[mask]
    
    # For example, compute a derived metric:
    delta_obs_ext_in_shell = dback * Sprk_ext_in_shell
    delta_obs_nom_in_shell = dback * Sprk_nom_in_shell
    delta_obs_sec_in_shell = flat_delta_obs_sec[mask]
    delta_obs_sec_in_shell_target = flat_delta_obs_sec_target[mask]
    #val_ext_in_shell = flat_val_ext[mask]
    #frac_ext_in_shell = flat_frac_ext[mask]
    #td = flat_td[mask]

    # Now filter the flattened extended mask arrays:
    #eta_ext_in_shell = (dback * Sprk_ext_in_shell * np.sqrt(td * ntr)) / (val_ext_in_shell * frac_ext_in_shell)

    # Suppose you want to compute an efficiency-like quantity:
    # Let’s say you define a boolean condition based on eta_nom_in_shell:
    nfp = (eta_nom_in_shell > eta_nom_threshold)
    # And another condition for the extended metric:
    nfp_ext = (eta_ext_in_shell > eta_ext_threshold) & (delta_obs_ext_in_shell > delta_obs_nom_in_shell) & (eta_nom_in_shell > eta_nom_threshold)
    nfp_sec = (eta_sec_in_shell > eta_ext_threshold) & (delta_obs_sec_in_shell > delta_obs_sec_in_shell_target) & (eta_nom_in_shell > eta_nom_threshold)
    nfp_ext_cob = (eta_cob_ext_in_shell > 3) & (eta_nom_in_shell > eta_nom_threshold)
    nfp_nom_cob = (eta_cob_nom_in_shell > 3) & (eta_nom_in_shell > eta_nom_threshold)
    count_nom = np.sum(nfp_ext)
    count_ext = np.sum(nfp)
    print(f"Shell {i}: count_nom = {count_nom}, count_ext = {count_ext}")
    # Then compute efficiency:
    eff_ext_flux = np.sum(nfp_ext) / np.sum(nfp) * 100.
    eff_ext_cob = np.sum(nfp_ext_cob) / np.sum(nfp) *100.
    eff_nom_cob = np.sum(nfp_nom_cob) / np.sum(nfp) *100.
    eff_sec_flux = np.sum(nfp_sec) / np.sum(nfp) * 100.


    
    # You can then plot or store the efficiency for each shell.
    plt.plot(r_lo, eff_ext_flux, 'ro')  # Example: a point per shell
    plt.plot(r_lo, eff_ext_cob, '^')
    plt.plot(r_lo, eff_nom_cob, 'P')
    plt.plot(r_lo, eff_sec_flux, 'd')

plt.xlabel("Distance (pixels)")
plt.ylabel("Efficiency (%)")
plt.show()
    

