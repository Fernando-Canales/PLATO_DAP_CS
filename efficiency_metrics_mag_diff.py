"""
Script that computes the efficiency
for the mag difference parameter space

Fernando 17.03.2025
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
flux_thresh_nom_mask, cob_thresh = 7.1, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1
eta_ext_threshold = 3

# Parameters for the eclipsing binaries (comment if you don't want fixed values for every contaminant)
transit_depth = 132000  # transit depth in ppm
td_ref = 6.72*0.46**2
ntr = 3        # number of transits in one hour
distance_max = 5 # maximum distance in pixels, from the target, to a star in the window in order to be considered a contaminant
#n_c_max = 300 # maximum number of contaminants in each window                      # processed PSFs 


mag = data_nominal_mask[:, 1] # mag of the targets
mag_2d = np.repeat(mag[:,np.newaxis], 10, axis=1)
mag_diff_from_target_to_10first_contaminants = data_nominal_mask[:, 169:179] - mag[:, None] # mag of the contaminants - mag of the target
flat_mag_diff = mag_diff_from_target_to_10first_contaminants.flatten()


nsr1h_sec = data_secondary_mask[:, 3]    # adjust index if needed
sig_depth_secondary_mask = nsr1h_sec * (1 - data_secondary_mask[:, 5]) / np.sqrt(td_ref * ntr)

sig_depth_nominal = data_nominal_mask[:, 7] * (1 - data_nominal_mask[:, 11]) / np.sqrt(td_ref * ntr)
# For extended:
sig_depth_extended = data_extended_mask[:, 4] * (1 - data_extended_mask[:, 13]) / np.sqrt(td_ref * ntr)
# Combined uncertainty:
sig_depth_total = np.sqrt(sig_depth_nominal**2 + sig_depth_extended**2)
sig_depth_total_secondary = np.sqrt(sig_depth_nominal**2 + sig_depth_secondary_mask**2)

# Now, we need to “repeat” these per–target values 10 times (one per contaminant).
# Assume the number of targets is:
sig_depth_total_reshaped = np.repeat(sig_depth_total[:, None], 10, axis=1)  # shape: (n_targets, 10)
# Then flatten this array to match the flattened version of your other metrics:
flat_sig_depth_total = sig_depth_total_reshaped.flatten()
flat_sig_depth_secondary = np.repeat(sig_depth_secondary_mask, 10)

# Initialize accumulators for overall counts:
total_nfp = 0
total_nfp_ext = 0
total_nfp_ext_cob = 0
total_nfp_nom_cob = 0


flat_mag = np.repeat(mag, 10)  # shape (N*10,)

# 2) Define uniform bins from flat_mag itself:
mdiff_min = np.min(flat_mag_diff)
mdiff_max = np.max(flat_mag_diff)
binsize = 0.5

# Number of bins (integer) that covers the entire range from mdiff_min to mdiff_max
n_bins = int((mdiff_max - mdiff_min) / binsize + 1)

# 3) QUANTILE BIN EDGES
num_bins = 5
clipped_diff = np.clip(flat_mag_diff, -2, None)
bin_edges = np.quantile(clipped_diff, np.linspace(0, 1, num_bins + 1))
print("Bin edges in Δm:", bin_edges)

# 3) Loop over bins:
for i in range(num_bins):

    lo = bin_edges[i]
    hi = bin_edges[i+1]
    
    # Build mask over flat_mag_diff
    mask = (flat_mag_diff >= lo) & (flat_mag_diff < hi)
    
    # Build mask for [left, right):
    #mask = (flat_mag >= left) & (flat_mag < right)
    # Now flatten the other arrays in the same way:
    flat_Sprk_ext = data_extended_mask[:, 14:24].flatten()  # shape: (n_targets*10,)
    flat_Sprk_nom = data_nominal_mask[:, 17:27].flatten()  # shape: (n_targets*10,)
    flat_dback = data_nominal_mask[:, 136:146].flatten()
    flat_td = data_nominal_mask[:, 126:136].flatten()
    flat_eta_cob_nom = data_nominal_mask[:, 46:56].flatten()
    flat_eta_cob_ext = data_extended_mask[:, 45:55].flatten()
    flat_eta_secondary = np.repeat(data_secondary_mask[:, 6], 10)
    flat_delta_obs_sec = np.repeat(data_secondary_mask[:, 7], 10)
    flat_eta_target_secondary = np.repeat(data_nominal_mask[:, 12], 10)
    flat_delta_obs_target_sec = np.repeat(data_nominal_mask[:, 13], 10)

    flat_eta_nom = data_nominal_mask[:, 27:37].flatten()
    flat_eta_ext = data_extended_mask[:, 24:34].flatten()

    # Filter the flattened eta values using the mask:
    eta_nom_in_bin = flat_eta_nom[mask]
    eta_ext_in_bin = flat_eta_ext[mask]

    Sprk_ext_in_bin = flat_Sprk_ext[mask]
    Sprk_nom_in_bin = flat_Sprk_nom[mask]
    eta_cob_nom_in_bin = flat_eta_cob_nom[mask]
    eta_cob_ext_in_bin = flat_eta_cob_ext[mask]
    eta_sec_in_bin = flat_eta_secondary[mask]
    dback = flat_dback[mask]
    
    # For example, compute a derived metric:
    delta_obs_ext_in_shell = dback * Sprk_ext_in_bin
    delta_obs_nom_in_shell = dback * Sprk_nom_in_bin
    delta_obs_sec_in_shell = flat_delta_obs_sec[mask]
    eta_target_secondary_in_shell = flat_eta_target_secondary[mask]
    delta_obst_target_sec_in_shell = flat_delta_obs_target_sec[mask]

    # Also filter the flattened sig_depth_total:
    sig_depth_in_shell = flat_sig_depth_total[mask]
    sig_depth_secondary_in_shell = flat_sig_depth_secondary[mask]

    # Let’s say you define a boolean condition based on eta_nom_in_bin:
    nfp = (eta_nom_in_bin > eta_nom_threshold)
    # And another condition for the extended metric:
    nfp_ext = (eta_ext_in_bin > eta_ext_threshold) & (delta_obs_ext_in_shell > delta_obs_nom_in_shell + 3*sig_depth_in_shell) & (eta_nom_in_bin > eta_nom_threshold)
    nfp_sec = (eta_sec_in_bin > eta_ext_threshold) & (delta_obs_sec_in_shell > delta_obst_target_sec_in_shell + 3*sig_depth_secondary_in_shell) & (eta_target_secondary_in_shell > eta_nom_threshold)
    nfp_ext_cob = (eta_cob_ext_in_bin > 3) & (eta_nom_in_bin > eta_nom_threshold)
    nfp_nom_cob = (eta_cob_nom_in_bin > 3) & (eta_nom_in_bin > eta_nom_threshold)
    count_nfp = np.sum(nfp)
    count_nfp_ext = np.sum(nfp_ext)
    count_nfp_ext_cob = np.sum(nfp_ext_cob)
    count_nfp_nom_cob = np.sum(nfp_nom_cob)

   
    # Then compute efficiency:
    eff_ext_flux = np.sum(nfp_ext) / np.sum(nfp) * 100.
    eff_ext_cob = np.sum(nfp_ext_cob) / np.sum(nfp) *100.
    eff_nom_cob = np.sum(nfp_nom_cob) / np.sum(nfp) *100.
    eff_sec_flux = np.sum(nfp_sec) / np.sum((eta_target_secondary_in_shell > eta_nom_threshold)) * 100.

      # (A) Bin midpoint for plotting:
    bin_center = (lo + hi)/2

    
    # Plot the efficiency points. Add label only for the first iteration.
    if i == 0:
        plt.plot(bin_center, eff_ext_flux, 'o-', color='blue', label='Extended flux')
        plt.plot(bin_center, eff_ext_cob, '^-', color='orange', label='Extended COB')
        plt.plot(bin_center, eff_nom_cob, 'P-', color='red', label='Nominal COB')
        plt.plot(bin_center, eff_sec_flux, 'd-', color='green', label='Sec flux')
    else:
        plt.plot(bin_center, eff_ext_flux, 'o-', color='blue')
        plt.plot(bin_center, eff_ext_cob, '^-', color='orange')
        plt.plot(bin_center, eff_nom_cob, 'P-', color='red')
        plt.plot(bin_center, eff_sec_flux, 'd-', color='green')

            # Sum counts for this shell:
    bin_nfp = np.sum(nfp)
    bin_nfp_ext = np.sum(nfp_ext)
    bin_nfp_ext_cob = np.sum(nfp_ext_cob)
    bin_nfp_nom_cob = np.sum(nfp_nom_cob)
    
    # Accumulate overall counts:
    total_nfp += bin_nfp # type: ignore
    total_nfp_ext += bin_nfp_ext # type: ignore
    total_nfp_ext_cob += bin_nfp_ext_cob
    total_nfp_nom_cob += bin_nfp_nom_cob

# After processing all shells, compute overall efficiencies:
if total_nfp > 0: # type: ignore
    overall_eff_ext_flux = total_nfp_ext / total_nfp * 100.
    overall_eff_ext_cob = total_nfp_ext_cob / total_nfp * 100.
    overall_eff_nom_cob = total_nfp_nom_cob / total_nfp * 100.
else:
    overall_eff_ext_flux = overall_eff_ext_cob = overall_eff_nom_cob = 0

print("Overall Extended Flux Efficiency: {:.2f}%".format(overall_eff_ext_flux))
print("Overall Extended COB Efficiency: {:.2f}%".format(overall_eff_ext_cob))
print("Overall Nominal COB Efficiency: {:.2f}%".format(overall_eff_nom_cob))


plt.xlabel("Contaminant mag − Target mag (Δm)")
plt.ylabel("Efficiency (%)")
plt.legend()
plt.savefig('/home/fercho/Documents/' +'efficiency_rings_mag_diff.pdf', dpi=300, format='pdf')
plt.show()
