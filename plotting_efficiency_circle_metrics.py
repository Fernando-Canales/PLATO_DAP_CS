"""
Script for computing the efficiency for detecting FPs
with each method and plotting median SPR_tot as a function of ring number.
Fernando 28th October 2024
"""
import numpy as np              # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.colors import LogNorm # type: ignore
import matplotlib.patheffects as PathEffects # type: ignore

# If you have a custom module 'fitting_psf' with 'from_pix_2_mm', uncomment the following line
# from fitting_psf import from_pix_2_mm 

# Some file parameters
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/rings/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'

# Initialize lists to store efficiencies for each ring
eff_extended_ring = []
eff_secondary_ring = []
eff_nom_cob_ring = []
eff_ext_cob_ring = []
eff_sec_cob_ring = []

# Lists to store efficiencies for masked targets
eff_extended_ring_masked = []
eff_secondary_ring_masked = []
eff_nom_cob_ring_masked = []
eff_ext_cob_ring_masked = []
eff_sec_cob_ring_masked = []

# Initialize lists to store SPR_tot metrics for each ring
spr_tot_nom_all = []        # For all targets
spr_tot_ext_all = []
spr_tot_sec_all = []

spr_tot_nom_masked = []     # For masked targets (eta_nom > 7.1)
spr_tot_ext_masked = []
spr_tot_sec_masked = []

# Initialize lists to store uncertainties for efficiencies
eff_nom_cob_ring_error = []
eff_ext_cob_ring_error = []
eff_sec_cob_ring_error = []
eff_secondary_ring_error = []
eff_extended_ring_error = []

# Initialize lists to store uncertainties for masked efficiencies
eff_nom_cob_ring_masked_error = []
eff_ext_cob_ring_masked_error = []
eff_sec_cob_ring_masked_error = []
eff_secondary_ring_masked_error = []
eff_extended_ring_masked_error = []

# Define parameters
number_of_circles = 7  # Adjust as needed
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Reference value for transit depth (ppm)
td_ref = 6.72 * 0.46**2  # Transit duration in hours
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other possible value: 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1

# Function to convert pixels to mm
def from_pix_2_mm(x_star, y_star):
    """
    Placeholder function for converting pixel coordinates to millimeters.
    Replace with actual implementation as needed.
    """
    # Example conversion factors (these should be defined based on your setup)
    pix_to_mm_x = 0.1  # mm per pixel in x-direction
    pix_to_mm_y = 0.1  # mm per pixel in y-direction
    x_mm = x_star * pix_to_mm_x
    y_mm = y_star * pix_to_mm_y
    return x_mm, y_mm

# Function to add concentric circles to the plot
def add_concentric_circles(ax, R, N):
    for i in range(N):
        r_i = (i + 1) / N * R  # Linear spacing for circle radii
        circle = plt.Circle((0, 0), r_i, color='darkorange', linestyle='--', fill=False, linewidth=2, zorder=3)
        ax.add_artist(circle)

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
    # Load .npy files for the current ring
    try:
        ring_nominal = np.load(f"{dataDIR}ring_{i}_nominal.npy")
        ring_secondary = np.load(f"{dataDIR}ring_{i}_secondary.npy")
        ring_extended = np.load(f"{dataDIR}ring_{i}_extended.npy")
    except FileNotFoundError as e:
        print(f"Error loading data for Ring {i}: {e}")
        continue  # Skip to the next ring if files are missing

    # Number of targets in this ring
    n_targets_ring = ring_nominal.shape[0]

    # Initialize arrays for this ring
    eta_ext_bt_24_cameras = np.zeros((n_targets_ring, 10))
    eta_ext_bt_6_cameras = np.zeros((n_targets_ring, 10))
    eta_nom_bt_24_cameras = np.zeros((n_targets_ring, 10))
    eta_nom_bt_6_cameras = np.zeros((n_targets_ring, 10))
    delta_obs = np.zeros((n_targets_ring, 10))
    delta_obs_ext = np.zeros((n_targets_ring, 10))
    sig_depth_24_cameras = np.zeros((n_targets_ring, 10))
    sig_depth_sec_nom_quad = np.zeros((n_targets_ring, 10))

    # Initialize boolean arrays for conditions
    nfp = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_ext_mask = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_nom_cob = np.zeros((n_targets_ring, 10), dtype=bool)
    nfp_ext_cob = np.zeros((n_targets_ring, 10), dtype=bool)
    fp_single_contaminant_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)
    secondary_mask_conditions_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)
    secondary_mask_cob_conditions_24_cameras = np.zeros((n_targets_ring, 10), dtype=bool)

    # Compute metrics for each target in the ring
    for j in range(n_targets_ring):
        dback = np.ones(10) * dback_ref  # Reference transit depth
        td = np.ones(10) * td_ref  # Reference transit duration

        # Ensure indices correspond to your data structure
        sprk_10first = ring_nominal[j, 17:27]   # Adjust indices if necessary
        try:
            nsr_1h_24_cameras_nominal_mask = ring_nominal[j, 7]  # Adjust index if necessary
            SPR_tot = ring_nominal[j, 11]  # Adjust index if necessary
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in SPR_tot extraction: {e}")
            SPR_tot = 0

        # Compute eta_nom_bt_24_cameras and eta_nom_bt_6_cameras
        try:
            eta_nom_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                / (ring_nominal[j, 7] * (1 - ring_nominal[j, 11]))
            )
            eta_nom_bt_6_cameras[j, :] = (
                gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                / (ring_nominal[j, 148] * (1 - ring_nominal[j, 11]))  # Adjust index if necessary
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_nom calculations: {e}")
            eta_nom_bt_24_cameras[j, :] = 0
            eta_nom_bt_6_cameras[j, :] = 0

        # Compute eta_ext_bt_24_cameras and eta_ext_bt_6_cameras
        try:
            eta_ext_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                / (ring_extended[j, 4] * (1 - ring_extended[j, 13]))
            )
            eta_ext_bt_6_cameras[j, :] = (
                gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                / (ring_extended[j, 44] * (1 - ring_extended[j, 13]))  # Adjust index if necessary
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_ext calculations: {e}")
            eta_ext_bt_24_cameras[j, :] = 0
            eta_ext_bt_6_cameras[j, :] = 0

        # Observed transit depth
        try:
            delta_obs[j, :] = dback * ring_nominal[j, 17:27]
            delta_obs_ext[j, :] = dback * ring_extended[j, 14:24]
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in observed transit depth extraction: {e}")
            delta_obs[j, :] = 0
            delta_obs_ext[j, :] = 0

        # COB-related variables
        try:
            eta_cob_nom_10first_24_cameras = ring_nominal[j, 46:56]
            eta_cob_ext_10first_24_cameras = ring_extended[j, 45:55]
            eta_cob_sec = ring_secondary[j, 9]
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in COB extraction: {e}")
            eta_cob_nom_10first_24_cameras = np.zeros(10)
            eta_cob_ext_10first_24_cameras = np.zeros(10)
            eta_cob_sec = 0

        # False Positive detection conditions
        nfp[j, :] = eta_nom_bt_24_cameras[j, :] > flux_thresh_nom_mask

        # Extended mask detection conditions
        # Assuming sig_depth_24_cameras is supposed to be computed as uncertainties
        # Placeholder: If sig_depth_24_cameras should be derived from delta_obs or other metrics, compute it here
        # For demonstration, let's assume it's proportional to delta_obs
        sig_depth_24_cameras[j, :] = delta_obs[j, :] * 0.1  # Example: 10% uncertainty

        nfp_ext_mask[j, :] = (
            (eta_ext_bt_24_cameras[j, :] > flux_thresh_ext_mask) &
            (delta_obs_ext[j, :] > delta_obs[j, :] + depth_sig_scaling * sig_depth_24_cameras[j, :])
        )

        # COB detection conditions
        nfp_nom_cob[j, :] = eta_cob_nom_10first_24_cameras > cob_thresh
        nfp_ext_cob[j, :] = eta_cob_ext_10first_24_cameras > cob_thresh

        # Single contaminant false positives
        fp_single_contaminant_24_cameras[j, :] = nfp[j, :]

        # Secondary mask conditions
        try:
            eta_c = ring_secondary[j, 6]
            delta_obs_c = ring_secondary[j, 7]
            # Assuming sig_depth_sec_nom_quad should be computed similarly
            sig_depth_sec_nom_quad[j, :] = delta_obs_c * 0.1  # Example: 10% uncertainty
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in secondary mask conditions: {e}")
            eta_c = 0
            delta_obs_c = 0
            sig_depth_sec_nom_quad[j, :] = 0

        secondary_mask_conditions_24_cameras[j, :] = (
            (eta_c > flux_thresh_sec_mask) &
            (delta_obs_c > delta_obs[j, :] + depth_sig_scaling * sig_depth_sec_nom_quad[j, :]) &
            fp_single_contaminant_24_cameras[j, :]
        )

        # Secondary mask COB conditions
        eta_cob_sec = ring_secondary[j, 9] if ring_secondary.shape[1] > 9 else 0
        secondary_mask_cob_conditions_24_cameras[j, :] = (
            (eta_cob_sec > cob_thresh) & fp_single_contaminant_24_cameras[j, :]
        )

    # --- New Code for SPR_tot Accumulation Starts Here ---

    # Extract SPR_tot_nom, SPR_tot_ext, and SPR_tot_sec for all targets in this ring
    try:
        SPR_tot_nom = ring_nominal[:, 11]   # Adjust index if necessary
        SPR_tot_ext = ring_extended[:, 13]  # Adjust index if necessary
        SPR_tot_sec = ring_secondary[:, 9]  # Adjust index if necessary
    except IndexError as e:
        print(f"Index error while extracting SPR_tot metrics for Ring {i}: {e}")
        SPR_tot_nom = np.array([])
        SPR_tot_ext = np.array([])
        SPR_tot_sec = np.array([])

    # Append to all targets lists
    spr_tot_nom_all.append(SPR_tot_nom)
    spr_tot_ext_all.append(SPR_tot_ext)
    spr_tot_sec_all.append(SPR_tot_sec)

    # Identify masked targets where any eta_nom_bt_24_cameras > threshold
    eta_nom_above_threshold = np.any(eta_nom_bt_24_cameras > eta_nom_threshold, axis=1)  # Shape (n_targets_ring,)
    SPR_tot_nom_masked = SPR_tot_nom[eta_nom_above_threshold]
    SPR_tot_ext_masked = SPR_tot_ext[eta_nom_above_threshold]
    SPR_tot_sec_masked = SPR_tot_sec[eta_nom_above_threshold]

    # Append to masked targets lists
    spr_tot_nom_masked.append(SPR_tot_nom_masked)
    spr_tot_ext_masked.append(SPR_tot_ext_masked)
    spr_tot_sec_masked.append(SPR_tot_sec_masked)

    # --- New Code for SPR_tot Accumulation Ends Here ---

    # --- Existing Code for Efficiency Calculations ---

    # Total number of False Positives (FP) for all targets
    nfp_total_all = nfp.sum()

    if nfp_total_all > 0:
        eff_ext_flux_all = (nfp & nfp_ext_mask).sum() / nfp_total_all * 100.
        eff_nom_cob_all = (nfp & nfp_nom_cob).sum() / nfp_total_all * 100.
        eff_ext_cob_all = (nfp & nfp_ext_cob).sum() / nfp_total_all * 100.
    else:
        eff_ext_flux_all = 0
        eff_nom_cob_all = 0
        eff_ext_cob_all = 0

    # Secondary mask efficiency for all targets
    fp_single_total_all = fp_single_contaminant_24_cameras.sum()
    if fp_single_total_all > 0:
        eff_secondary_all = secondary_mask_conditions_24_cameras.sum() / fp_single_total_all * 100.
        eff_sec_cob_all = secondary_mask_cob_conditions_24_cameras.sum() / fp_single_total_all * 100.
    else:
        eff_secondary_all = 0
        eff_sec_cob_all = 0

    # Store efficiencies for all targets
    eff_extended_ring.append(eff_ext_flux_all)
    eff_secondary_ring.append(eff_secondary_all)
    eff_nom_cob_ring.append(eff_nom_cob_all)
    eff_ext_cob_ring.append(eff_ext_cob_all)
    eff_sec_cob_ring.append(eff_sec_cob_all)

    # Compute uncertainties for all targets using binomial error
    delta_eff_ext_flux_all = binomial_error(int((nfp & nfp_ext_mask).sum()), int(nfp_total_all))
    eff_extended_ring_error.append(delta_eff_ext_flux_all)

    delta_eff_nom_cob_all = binomial_error(int((nfp & nfp_nom_cob).sum()), int(nfp_total_all))
    eff_nom_cob_ring_error.append(delta_eff_nom_cob_all)

    delta_eff_ext_cob_all = binomial_error(int((nfp & nfp_ext_cob).sum()), int(nfp_total_all))
    eff_ext_cob_ring_error.append(delta_eff_ext_cob_all)

    delta_eff_secondary_all = binomial_error(int(secondary_mask_conditions_24_cameras.sum()), int(fp_single_total_all))
    eff_secondary_ring_error.append(delta_eff_secondary_all)

    delta_eff_sec_cob_all = binomial_error(int(secondary_mask_cob_conditions_24_cameras.sum()), int(fp_single_total_all))
    eff_sec_cob_ring_error.append(delta_eff_sec_cob_all)

    # Now compute efficiencies only for targets where eta_nom_bt_24_cameras > threshold

    if np.any(eta_nom_above_threshold):
        nfp_masked = nfp[eta_nom_above_threshold, :]
        nfp_ext_mask_masked = nfp_ext_mask[eta_nom_above_threshold, :]
        nfp_nom_cob_masked = nfp_nom_cob[eta_nom_above_threshold, :]
        nfp_ext_cob_masked = nfp_ext_cob[eta_nom_above_threshold, :]
        fp_single_contaminant_masked = fp_single_contaminant_24_cameras[eta_nom_above_threshold, :]
        secondary_mask_conditions_masked = secondary_mask_conditions_24_cameras[eta_nom_above_threshold, :]
        secondary_mask_cob_conditions_masked = secondary_mask_cob_conditions_24_cameras[eta_nom_above_threshold, :]
    else:
        nfp_masked = np.array([])
        nfp_ext_mask_masked = np.array([])
        nfp_nom_cob_masked = np.array([])
        nfp_ext_cob_masked = np.array([])
        fp_single_contaminant_masked = np.array([])
        secondary_mask_conditions_masked = np.array([])
        secondary_mask_cob_conditions_masked = np.array([])

    # Total number of False Positives (FP) for masked targets
    nfp_total_masked = nfp_masked.sum() if nfp_masked.size > 0 else 0

    if nfp_total_masked > 0:
        eff_ext_flux_masked = (nfp_masked & nfp_ext_mask_masked).sum() / nfp_total_masked * 100.
        eff_nom_cob_masked = (nfp_masked & nfp_nom_cob_masked).sum() / nfp_total_masked * 100.
        eff_ext_cob_masked = (nfp_masked & nfp_ext_cob_masked).sum() / nfp_total_masked * 100.
    else:
        eff_ext_flux_masked = 0
        eff_nom_cob_masked = 0
        eff_ext_cob_masked = 0

    # Secondary mask efficiency for masked targets
    fp_single_total_masked = fp_single_contaminant_masked.sum() if fp_single_contaminant_masked.size > 0 else 0
    if fp_single_total_masked > 0:
        eff_secondary_masked = secondary_mask_conditions_masked.sum() / fp_single_total_masked * 100.
        eff_sec_cob_masked = secondary_mask_cob_conditions_masked.sum() / fp_single_total_masked * 100.
    else:
        eff_secondary_masked = 0
        eff_sec_cob_masked = 0

    # Store efficiencies for masked targets
    eff_extended_ring_masked.append(eff_ext_flux_masked)
    eff_secondary_ring_masked.append(eff_secondary_masked)
    eff_nom_cob_ring_masked.append(eff_nom_cob_masked)
    eff_ext_cob_ring_masked.append(eff_ext_cob_masked)
    eff_sec_cob_ring_masked.append(eff_sec_cob_masked)

    # Compute uncertainties for masked targets using binomial error
    delta_eff_ext_flux_masked = binomial_error(int((nfp_masked & nfp_ext_mask_masked).sum()), int(nfp_total_masked))
    eff_extended_ring_masked_error.append(delta_eff_ext_flux_masked)

    delta_eff_nom_cob_masked = binomial_error(int((nfp_masked & nfp_nom_cob_masked).sum()), int(nfp_total_masked))
    eff_nom_cob_ring_masked_error.append(delta_eff_nom_cob_masked)

    delta_eff_ext_cob_masked = binomial_error(int((nfp_masked & nfp_ext_cob_masked).sum()), int(nfp_total_masked))
    eff_ext_cob_ring_masked_error.append(delta_eff_ext_cob_masked)

    delta_eff_secondary_masked = binomial_error(int(secondary_mask_conditions_masked.sum()), int(fp_single_total_masked))
    eff_secondary_ring_masked_error.append(delta_eff_secondary_masked)

    delta_eff_sec_cob_masked = binomial_error(int(secondary_mask_cob_conditions_masked.sum()), int(fp_single_total_masked))
    eff_sec_cob_ring_masked_error.append(delta_eff_sec_cob_masked)

    # Print efficiencies for the current ring with error bars
    print(f"Ring {i} Efficiencies for All Targets:")
    print(f"  Extended Flux Efficiency: {eff_ext_flux_all:.2f}% ± {delta_eff_ext_flux_all:.2f}%")
    print(f"  Secondary Mask Efficiency: {eff_secondary_all:.2f}% ± {delta_eff_secondary_all:.2f}%")
    print(f"  Nominal COB Efficiency: {eff_nom_cob_all:.2f}% ± {delta_eff_nom_cob_all:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob_all:.2f}% ± {delta_eff_ext_cob_all:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_sec_cob_all:.2f}% ± {delta_eff_sec_cob_all:.2f}%")

    print(f"Ring {i} Efficiencies for Targets with eta_nom_bt_24_cameras > {eta_nom_threshold}:")
    print(f"  Extended Flux Efficiency: {eff_ext_flux_masked:.2f}% ± {delta_eff_ext_flux_masked:.2f}%")
    print(f"  Secondary Mask Efficiency: {eff_secondary_masked:.2f}% ± {delta_eff_secondary_masked:.2f}%")
    print(f"  Nominal COB Efficiency: {eff_nom_cob_masked:.2f}% ± {delta_eff_nom_cob_masked:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob_masked:.2f}% ± {delta_eff_ext_cob_masked:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_sec_cob_masked:.2f}% ± {delta_eff_sec_cob_masked:.2f}%")
    print("-----------------------------------------------------")

# After the ring loop

# Compute median SPR_tot for all targets per ring
median_spr_tot_nom_all = [np.median(spr_tot_nom_all[ring]) for ring in range(number_of_circles)]
median_spr_tot_ext_all = [np.median(spr_tot_ext_all[ring]) for ring in range(number_of_circles)]
median_spr_tot_sec_all = [np.median(spr_tot_sec_all[ring]) for ring in range(number_of_circles)]

# Compute median SPR_tot for masked targets per ring
median_spr_tot_nom_masked = [
    np.median(spr_tot_nom_masked[ring]) if len(spr_tot_nom_masked[ring]) > 0 else np.nan 
    for ring in range(number_of_circles)
]
median_spr_tot_ext_masked = [
    np.median(spr_tot_ext_masked[ring]) if len(spr_tot_ext_masked[ring]) > 0 else np.nan 
    for ring in range(number_of_circles)
]
median_spr_tot_sec_masked = [
    np.median(spr_tot_sec_masked[ring]) if len(spr_tot_sec_masked[ring]) > 0 else np.nan 
    for ring in range(number_of_circles)
]

# Define ring numbers for x-axis (1 to number_of_circles)
ring_numbers = np.arange(1, number_of_circles + 1)

# ----- Plot 1: Median SPR_tot vs Ring Number for All Targets -----
plt.figure(figsize=(10, 6), dpi=120)

plt.plot(
    ring_numbers, 
    median_spr_tot_nom_all, 
    marker='o', 
    linestyle='-', 
    color='brown', 
    label='SPR_tot_nom (All Targets)'
)
plt.plot(
    ring_numbers, 
    median_spr_tot_ext_all, 
    marker='s', 
    linestyle='-', 
    color='blue', 
    label='SPR_tot_ext (All Targets)'
)
plt.plot(
    ring_numbers, 
    median_spr_tot_sec_all, 
    marker='^', 
    linestyle='-', 
    color='green', 
    label='SPR_tot_sec (All Targets)'
)

plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'$\rm SPR_{tot}$', fontsize=14)
plt.title("Median SPR_tot as a Function of Ring Number (All Targets)", fontsize=16)
plt.xticks(ring_numbers, [f"Ring {i}" for i in ring_numbers], fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)

plt.tight_layout()
plt.savefig("Median_SPR_tot_vs_Ring_All_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ----- Plot 2: Median SPR_tot vs Ring Number for Masked Targets -----
plt.figure(figsize=(10, 6), dpi=120)

plt.plot(
    ring_numbers, 
    median_spr_tot_nom_masked, 
    marker='o', 
    linestyle='-', 
    color='red', 
    label=f'SPR_tot_nom (eta_nom > {eta_nom_threshold})'
)
plt.plot(
    ring_numbers, 
    median_spr_tot_ext_masked, 
    marker='s', 
    linestyle='-', 
    color='purple', 
    label=f'SPR_tot_ext (eta_nom > {eta_nom_threshold})'
)
plt.plot(
    ring_numbers, 
    median_spr_tot_sec_masked, 
    marker='^', 
    linestyle='-', 
    color='orange', 
    label=f'SPR_tot_sec (eta_nom > {eta_nom_threshold})'
)

plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'$\rm SPR_{tot}$', fontsize=14)
plt.title(f"Median SPR_tot as a Function of Ring Number (eta_nom_bt_24_cameras > {eta_nom_threshold})", fontsize=16)
plt.xticks(ring_numbers, [f"Ring {i}" for i in ring_numbers], fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)

plt.tight_layout()
plt.savefig("Median_SPR_tot_vs_Ring_Masked_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ----- Plot 3: Efficiencies Across Rings for All Targets with Error Bars -----
plt.figure(figsize=(10, 6), dpi=120)

plt.errorbar(
    ring_numbers, 
    eff_nom_cob_ring, 
    yerr=eff_nom_cob_ring_error, 
    fmt='o-', 
    color='brown', 
    label="Nominal COB Efficiency (All Targets)",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_ext_cob_ring, 
    yerr=eff_ext_cob_ring_error, 
    fmt='s-', 
    color='blue', 
    label="Extended COB Efficiency (All Targets)",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_sec_cob_ring, 
    yerr=eff_sec_cob_ring_error, 
    fmt='^-', 
    color='purple', 
    label="Secondary COB Efficiency (All Targets)",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_secondary_ring, 
    yerr=eff_secondary_ring_error, 
    fmt='^-', 
    color='green', 
    label="Secondary Flux Efficiency (All Targets)",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_extended_ring, 
    yerr=eff_extended_ring_error, 
    fmt='d-', 
    color='red', 
    label="Extended Flux Efficiency (All Targets)",
    capsize=5
)

# Adjust x-axis labels
circle_labels = [f"Ring {i}" for i in ring_numbers]  # Ring 1 to Ring 7
plt.xticks(ring_numbers, circle_labels, fontsize=10)

plt.ylabel("Efficiency (%)", fontsize=14)
plt.title("Efficiencies Across Rings (All Targets)", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Efficiencies_Across_Rings_All_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ----- Plot 4: Efficiencies Across Rings for Masked Targets with Error Bars -----
plt.figure(figsize=(10, 6), dpi=120)

plt.errorbar(
    ring_numbers, 
    eff_nom_cob_ring_masked, 
    yerr=eff_nom_cob_ring_masked_error, 
    fmt='o-', 
    color='brown', 
    label=f"Nominal COB Efficiency (eta_nom > {eta_nom_threshold})",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_ext_cob_ring_masked, 
    yerr=eff_ext_cob_ring_masked_error, 
    fmt='s-', 
    color='blue', 
    label=f"Extended COB Efficiency (eta_nom > {eta_nom_threshold})",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_sec_cob_ring_masked, 
    yerr=eff_sec_cob_ring_masked_error, 
    fmt='^-', 
    color='purple', 
    label=f"Secondary COB Efficiency (eta_nom > {eta_nom_threshold})",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_secondary_ring_masked, 
    yerr=eff_secondary_ring_masked_error, 
    fmt='^-', 
    color='green', 
    label=f"Secondary Flux Efficiency (eta_nom > {eta_nom_threshold})",
    capsize=5
)
plt.errorbar(
    ring_numbers, 
    eff_extended_ring_masked, 
    yerr=eff_extended_ring_masked_error, 
    fmt='d-', 
    color='red', 
    label=f"Extended Flux Efficiency (eta_nom > {eta_nom_threshold})",
    capsize=5
)

# Adjust x-axis labels
plt.xticks(ring_numbers, circle_labels, rotation=45, fontsize=10)

plt.ylabel("Efficiency (%)", fontsize=14)
plt.title(f"Efficiencies Across Rings (eta_nom_bt_24_cameras > {eta_nom_threshold})", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Efficiencies_Across_Rings_Masked_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ----- Existing Plotting Code for Efficiencies Across Rings -----
# (If you have other plots like quadrant plots, include them here with similar error bar additions)