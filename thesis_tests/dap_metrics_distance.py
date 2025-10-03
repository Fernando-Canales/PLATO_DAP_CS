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

distance_from_target_to_10first_contaminants = data_nominal_mask[:, 179:189] # within the imagette, in pixels
flat_distance = distance_from_target_to_10first_contaminants.flatten()


# Instead of radius_focal_plane, we use our fixed maximum distance for the contaminants:
max_distance = 5  # maximum distance in pixels
number_of_shells = 5
shell_edges = np.sqrt(np.linspace(0, max_distance**2, number_of_shells + 1))
print("Shell edges (pixels):", shell_edges)

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

#ext_eff = [66.51, 70.22, 70.69, 70.76, 70.80]

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
    eta_target_secondary_in_shell = flat_eta_target_secondary[mask]
    delta_obst_target_sec_in_shell = flat_delta_obs_target_sec[mask]

    # Also filter the flattened sig_depth_total:
    sig_depth_in_shell = flat_sig_depth_total[mask]
    sig_depth_secondary_in_shell = flat_sig_depth_secondary[mask]

    # Let’s say you define a boolean condition based on eta_nom_in_shell:
    nfp = (eta_nom_in_shell > eta_nom_threshold)
    # And another condition for the extended metric:
    nfp_ext = (eta_ext_in_shell > eta_ext_threshold) & (delta_obs_ext_in_shell > delta_obs_nom_in_shell + 3*sig_depth_in_shell) & (eta_nom_in_shell > eta_nom_threshold)
    nfp_sec = (eta_sec_in_shell > eta_ext_threshold) & (delta_obs_sec_in_shell > delta_obst_target_sec_in_shell + 3*sig_depth_secondary_in_shell) & (eta_target_secondary_in_shell > eta_nom_threshold)
    nfp_ext_cob = (eta_cob_ext_in_shell > 3) & (eta_nom_in_shell > eta_nom_threshold)
    nfp_nom_cob = (eta_cob_nom_in_shell > 3) & (eta_nom_in_shell > eta_nom_threshold)
    count_nfp = np.sum(nfp)
    count_nfp_ext = np.sum(nfp_ext)
    count_nfp_ext_cob = np.sum(nfp_ext_cob)
    count_nfp_nom_cob = np.sum(nfp_nom_cob)

    print(f"Shell {i}: nfp = {count_nfp}, nfp_ext = {count_nfp_ext} nfp_ext_cob = {count_nfp_ext_cob}, nfp_nom_cob = {count_nfp_nom_cob}")
    # Then compute efficiency:
    eff_ext_flux = np.sum(nfp_ext) / np.sum(nfp) * 100.
    eff_ext_cob = np.sum(nfp_ext_cob) / np.sum(nfp) *100.
    eff_nom_cob = np.sum(nfp_nom_cob) / np.sum(nfp) *100.
    eff_sec_flux = np.sum(nfp_sec) / np.sum((eta_target_secondary_in_shell > eta_nom_threshold)) * 100.

    r = (r_lo + r_hi) / 2
    
    # Plot the efficiency points. Add label only for the first iteration.
    if i == 0:
        #plt.plot(r_lo, ext_eff[i], 'o-', color='blue', label='Extended flux')
        plt.plot(r_lo, eff_ext_cob, '^-', color='orange', label='Extended COB')
        plt.plot(r_lo, eff_nom_cob, 'P-', color='red', label='Nominal COB')
        plt.plot(r, eff_sec_flux, 'd-', color='green', label='Sec flux')
    else:
        #plt.plot(r_lo, ext_eff[i], 'o-', color='blue')
        plt.plot(r_lo, eff_ext_cob, '^-', color='orange')
        plt.plot(r_lo, eff_nom_cob, 'P-', color='red')
        plt.plot(r_lo, eff_sec_flux, 'd-', color='green')

            # Sum counts for this shell:
    shell_nfp = np.sum(nfp)
    shell_nfp_ext = np.sum(nfp_ext)
    shell_nfp_ext_cob = np.sum(nfp_ext_cob)
    shell_nfp_nom_cob = np.sum(nfp_nom_cob)
    
    print(f"Shell {i}: nfp = {shell_nfp}, nfp_ext = {shell_nfp_ext}, "
          f"nfp_ext_cob = {shell_nfp_ext_cob}, nfp_nom_cob = {shell_nfp_nom_cob}")
    
    # Accumulate overall counts:
    total_nfp += shell_nfp # type: ignore
    total_nfp_ext += shell_nfp_ext # type: ignore
    print(total_nfp_ext/total_nfp)
    total_nfp_ext_cob += shell_nfp_ext_cob
    total_nfp_nom_cob += shell_nfp_nom_cob

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


plt.xlabel("Shells")
plt.ylabel("Efficiency (%)")
plt.legend()
plt.savefig('/home/fercho/Documents/' +'efficiency_rings_distance.pdf', dpi=300, format='pdf')
plt.show()
