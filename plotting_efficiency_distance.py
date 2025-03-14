"""
Script that computes the efficieincy of
the different metrics (dap and centroid shift)
in terms of the distance from the target to
the 10 first contaminants in terms of SPR

Fernando 12.03.2025
"""
import numpy as np              # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.colors import LogNorm # type: ignore
import matplotlib.patheffects as PathEffects # type: ignore

# Import parameters
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/rings/distances/'

dataDIRnominal = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_4_22_hr_CONDITION_1-PIXEL_SEC_MASK/'
data_nominal_mask = np.load(dataDIRnominal + 'targets_P5.npy')

distance_from_target_to_10first_contaminants = data_nominal_mask[:, 179:189] # within the imagette, in pixels
# Compute a representative distance for each target (e.g. the minimum distance among the 10)
rep_distance = np.median(distance_from_target_to_10first_contaminants, axis=1)
distance_from_target_to_10first_contaminants_1D = distance_from_target_to_10first_contaminants.flatten()

# Define ring numbers for x-axis (1 to number_of_circles)
number_of_circles = 7 
ring_numbers = np.arange(1, number_of_circles + 1)  # [1, 2, ..., number_of_circles]

# Initialize lists to store efficiencies for each ring
# Initialize lists to store efficiencies for each ring
eff_extended_ring = []
eff_secondary_ring = []
eff_nom_cob_ring = []
eff_ext_cob_ring = []
eff_sec_cob_ring = []

# Initialize lists to store uncertainties for efficiencies
eff_nom_cob_ring_error = []
eff_ext_cob_ring_error = []
eff_sec_cob_ring_error = []
eff_secondary_ring_error = []
eff_extended_ring_error = []

# Define parameters
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Reference transit depth in ppm
td_ref = 6.72 * 0.46**2  # Transit duration in hours
depth_sig_scaling = 3
gamma_factor_significance = 1  # (Alternate value: 0.46)
flux_thresh_nom_mask, cob_thresh = 7.1, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1

# Function to compute binomial confidence interval (standard error)
def binomial_error(k, n):
    """
    Computes the binomial confidence interval (standard error) for a proportion.
    Args:
        k (int): Number of successes
        n (int): Number of trials
    Returns:
        float: Standard error in percentage points
    """
    if n > 0:
        p = k / n
        return np.sqrt(p * (1 - p) / n) * 100
    else:
        return 0.0

# Loop through each ring, load corresponding .npy files, and compute metrics
for i in range(number_of_circles):
    # Load .npy files for the current ring (assume these files exist)
    ring_nominal = np.load(f"{dataDIR}ring_{i}_nominal.npy")
    ring_secondary = np.load(f"{dataDIR}ring_{i}_secondary.npy")
    ring_extended = np.load(f"{dataDIR}ring_{i}_extended.npy")
    
    # Number of targets in this ring
    n_targets_ring = ring_nominal.shape[0]
    
    # Initialize arrays for this ring (for example, for dap and centroid metrics)
    eta_ext_bt_24_cameras = np.zeros((n_targets_ring, 10))
    eta_ext_bt_6_cameras = np.zeros((n_targets_ring, 10))
    eta_nom_bt_24_cameras = np.zeros((n_targets_ring, 10))
    eta_nom_bt_6_cameras = np.zeros((n_targets_ring, 10))
    delta_obs = np.zeros((n_targets_ring, 10))
    delta_obs_ext = np.zeros((n_targets_ring, 10))
    sig_depth_24_cameras = np.zeros((n_targets_ring, 10))
    sig_depth_6_cameras = np.zeros((n_targets_ring, 10))
    sig_depth_sec_24_cameras = np.zeros((n_targets_ring, 10))
    
    # Initialize boolean arrays for conditions
    nfp = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_ext_mask = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_nom_cob = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_ext_cob = np.zeros((n_targets_ring, 10), dtype=bool)
    fp_single_contaminant_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)
    secondary_mask_conditions_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)
    secondary_mask_cob_conditions_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)
    
    # Compute metrics for each target in the current ring
    for j in range(n_targets_ring):
        dback = np.ones(10) * dback_ref
        td = np.ones(10) * td_ref
        
        # Extract necessary values from the loaded ring arrays (adjust indices as needed)
        sprk_10first = ring_nominal[j, 17:27]
        nsr_1h_24_cameras_nominal_mask = ring_nominal[j, 7]
        SPR_tot = ring_nominal[j, 11]
        
        eta_nom_bt_24_cameras[j, :] = (gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                                       / (ring_nominal[j, 7] * (1 - ring_nominal[j, 11])))
        eta_nom_bt_6_cameras[j, :] = (gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                                      / (ring_nominal[j, 148] * (1 - ring_nominal[j, 11])))
        eta_ext_bt_24_cameras[j, :] = (gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                                       / (ring_extended[j, 4] * (1 - ring_extended[j, 13])))
        eta_ext_bt_6_cameras[j, :] = (gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                                      / (ring_extended[j, 44] * (1 - ring_extended[j, 13])))
        delta_obs[j, :] = dback * ring_nominal[j, 17:27]
        delta_obs_ext[j, :] = dback * ring_extended[j, 14:24]
        
        # COB-related variables
        eta_cob_nom_10first_24_cameras = ring_nominal[j, 46:56]
        eta_cob_ext_10first_24_cameras = ring_extended[j, 45:55]
        eta_cob_sec = ring_secondary[j, 9]

        
        
        nfp[j, :] = eta_nom_bt_24_cameras[j, :] > flux_thresh_nom_mask
        
        # Compute significant transit depth values:
        sig_depth_extended_mask_24_cameras = ring_extended[j, 4] * (1 - ring_extended[j, 13]) / np.sqrt(td_ref * ntr)
        sig_depth_extended_mask_6_cameras = ring_extended[j, 44] * (1 - ring_extended[j, 13]) / np.sqrt(td_ref * ntr)
        sig_depth_nominal_mask_24_cameras = ring_nominal[j, 7] * (1 - ring_extended[j, 11]) / np.sqrt(td_ref * ntr)
        sig_depth_nominal_mask_6_cameras = ring_nominal[j, 148] * (1 - ring_extended[j, 11]) / np.sqrt(td_ref * ntr)
        
        # Here we assume that sig_depth_nominal_mask_24_cameras and sig_depth_extended_mask_24_cameras
        # are scalars for the given target; if they are arrays, remove the [j, :] indexing.
        sig_depth_24_cameras[j, :] = np.sqrt(sig_depth_nominal_mask_24_cameras**2 +
                                              sig_depth_extended_mask_24_cameras**2)
        sig_depth_6_cameras[j, :] = np.sqrt(sig_depth_nominal_mask_6_cameras**2 +
                                             sig_depth_extended_mask_6_cameras**2)
        
        # Extended mask detection condition:
        nfp_ext_mask[j, :] = ((eta_ext_bt_24_cameras[j, :] > flux_thresh_ext_mask) &
                              (delta_obs_ext[j, :] > delta_obs[j, :] + depth_sig_scaling * sig_depth_24_cameras[j, :]))
        
        nfp_nom_cob[j, :] = eta_cob_nom_10first_24_cameras > cob_thresh
        nfp_ext_cob[j, :] = eta_cob_ext_10first_24_cameras > cob_thresh
        
        fp_single_contaminant_24_cameras[j, :] = nfp[j, :]
        
        eta_c = ring_secondary[j, 6]
        delta_obs_c = ring_secondary[j, 7]
        sig_depth_sec_24_cameras[j, :] = ring_secondary[j, 4] * (1 - ring_secondary[j, 5]) / np.sqrt(td_ref * ntr)
        
        secondary_mask_conditions_24_cameras[j, :] = ((eta_c > flux_thresh_sec_mask) &
                                                      (delta_obs_c > delta_obs[j, :] + depth_sig_scaling * sig_depth_sec_24_cameras[j, :]) &
                                                      fp_single_contaminant_24_cameras[j, :])
        
        eta_cob_sec = ring_secondary[j, 9]  # assuming shape[1] > 9
        secondary_mask_cob_conditions_24_cameras[j, :] = ((eta_cob_sec > cob_thresh) & fp_single_contaminant_24_cameras[j, :])
    
    
    # Compute total number of false positives (FP) for all targets in the ring
    nfp_total_all = nfp.sum()
    if nfp_total_all > 0:
        eff_ext_flux_all = (nfp & nfp_ext_mask).sum() / nfp_total_all * 100.
        eff_nom_cob_all = (nfp & nfp_nom_cob).sum() / nfp_total_all * 100.
        eff_ext_cob_all = (nfp & nfp_ext_cob).sum() / nfp_total_all * 100.
    else:
        eff_ext_flux_all = 0
        eff_nom_cob_all = 0
        eff_ext_cob_all = 0

    fp_single_total_all = fp_single_contaminant_24_cameras.sum()
    if fp_single_total_all > 0:
        eff_secondary_all = secondary_mask_conditions_24_cameras.sum() / fp_single_total_all * 100.
        eff_sec_cob_all = secondary_mask_cob_conditions_24_cameras.sum() / fp_single_total_all * 100.
    else:
        eff_secondary_all = 0
        eff_sec_cob_all = 0

    # Store efficiencies for the ring
    eff_extended_ring.append(eff_ext_flux_all)
    eff_secondary_ring.append(eff_secondary_all)
    eff_nom_cob_ring.append(eff_nom_cob_all)
    eff_ext_cob_ring.append(eff_ext_cob_all)
    eff_sec_cob_ring.append(eff_sec_cob_all)
    
    # Compute uncertainties using the binomial error function
    delta_eff_ext_flux_all = binomial_error(int(np.logical_and(nfp, nfp_ext_mask).sum()), int(nfp_total_all))
    delta_eff_nom_cob_all = binomial_error(int(np.logical_and(nfp, nfp_nom_cob).sum()), int(nfp_total_all))
    delta_eff_ext_cob_all = binomial_error(int(np.logical_and(nfp, nfp_ext_cob).sum()), int(nfp_total_all))
    delta_eff_secondary_all = binomial_error(int(secondary_mask_conditions_24_cameras.sum()), int(fp_single_total_all))
    delta_eff_sec_cob_all = binomial_error(int(secondary_mask_cob_conditions_24_cameras.sum()), int(fp_single_total_all))
    
    eff_extended_ring_error.append(delta_eff_ext_flux_all)
    eff_nom_cob_ring_error.append(delta_eff_nom_cob_all)
    eff_ext_cob_ring_error.append(delta_eff_ext_cob_all)
    eff_secondary_ring_error.append(delta_eff_secondary_all)
    eff_sec_cob_ring_error.append(delta_eff_sec_cob_all)
    
    print(f"Ring {i} Efficiencies for All Targets:")
    print(f"  Extended Flux Efficiency: {eff_ext_flux_all:.2f}% ± {delta_eff_ext_flux_all:.2f}%")
    print(f"  Secondary Mask Efficiency: {eff_secondary_all:.2f}% ± {delta_eff_secondary_all:.2f}%")
    print(f"  Nominal COB Efficiency: {eff_nom_cob_all:.2f}% ± {delta_eff_nom_cob_all:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob_all:.2f}% ± {delta_eff_ext_cob_all:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_sec_cob_all:.2f}% ± {delta_eff_sec_cob_all:.2f}%")
    
# Now, you can plot the efficiencies vs. ring number.
plt.figure(1, figsize=(10, 6))
plt.errorbar(ring_numbers, eff_extended_ring, yerr=eff_extended_ring_error, fmt='o-', label='Extended Flux Efficiency', capsize=5)
plt.errorbar(ring_numbers, eff_secondary_ring, yerr=eff_secondary_ring_error, fmt='s-', label='Secondary Mask Efficiency', capsize=5)
plt.errorbar(ring_numbers, eff_nom_cob_ring, yerr=eff_nom_cob_ring_error, fmt='^-', label='Nominal COB Efficiency', capsize=5)
plt.errorbar(ring_numbers, eff_ext_cob_ring, yerr=eff_ext_cob_ring_error, fmt='d-', label='Extended COB Efficiency', capsize=5)
plt.errorbar(ring_numbers, eff_sec_cob_ring, yerr=eff_sec_cob_ring_error, fmt='v-', label='Secondary COB Efficiency', capsize=5)

plt.xlabel("Ring Number")
plt.ylabel("Efficiency (%)")
plt.title("Efficiency vs. Distance Shell (Ring)")
plt.legend()
plt.xlim(1,6)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()


distance_from_target_to_10first_contaminants = data_nominal_mask[:, 179:189] # within the imagette, in pixels

# Instead of radius_focal_plane, we use our fixed maximum distance for the contaminants:
max_distance = 7.0  # maximum distance in pixels
number_of_shells = 7
shell_edges = np.sqrt(np.linspace(0, max_distance**2, number_of_shells + 1))
print("Shell edges (pixels):", shell_edges)

plt.figure(2)
plt.hist(distance_from_target_to_10first_contaminants_1D, bins=shell_edges, edgecolor='black', color='skyblue')
for edge in shell_edges:
    plt.axvline(edge, color='red', linestyle='--', linewidth=1)
    plt.xlabel("Distance (pixels)")
    plt.ylabel("Count")
    plt.title("Histogram of Distances with Shell Boundaries")

plt.show()