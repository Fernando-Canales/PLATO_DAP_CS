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

# Define ring numbers for x-axis (1 to number_of_circles)
number_of_circles = 7  # Ensure this matches your data
ring_numbers = np.arange(1, number_of_circles + 1)  # [1, 2, ..., number_of_circles]

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

# Initialize lists to store median SPR_tot for masked targets per ring
median_spr_tot_nom_masked_ring = []
median_spr_tot_ext_masked_ring = []
median_spr_tot_sec_masked_ring = []

# Initialize lists to store errors for median SPR_tot for masked targets
median_spr_tot_nom_masked_ring_error = []
median_spr_tot_ext_masked_ring_error = []
median_spr_tot_sec_masked_ring_error = []

# Initialize lists to store median SPR_tot for all targets per ring
median_spr_tot_nom_ring = []
median_spr_tot_ext_ring = []
median_spr_tot_sec_ring = []

# Initialize lists to store errors for median SPR_tot for all targets
median_spr_tot_nom_ring_error = []
median_spr_tot_ext_ring_error = []
median_spr_tot_sec_ring_error = []

# Define parameters
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

# Function to estimate median error via bootstrapping
def bootstrap_median(data, num_samples=1000, random_seed=None):
    """
    Estimates the standard error of the median using bootstrapping.
    
    Args:
        data (array-like): The data to bootstrap.
        num_samples (int): Number of bootstrap samples.
        random_seed (int, optional): Seed for reproducibility.
    
    Returns:
        float: Estimated standard error of the median.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    medians = []
    n = len(data)
    if n == 0:
        return np.nan  # Return NaN if data is empty
    
    for _ in range(num_samples):
        sample = np.random.choice(data, size=n, replace=True)
        medians.append(np.median(sample))
    
    return np.std(medians)

# Loop through each ring, load corresponding .npy files, and compute metrics
for i in range(number_of_circles):
    # Load .npy files for the current ring
    try:
        ring_nominal = np.load(f"{dataDIR}ring_{i}_nominal.npy")
        ring_secondary = np.load(f"{dataDIR}ring_{i}_secondary.npy")
        ring_extended = np.load(f"{dataDIR}ring_{i}_extended.npy")
    except FileNotFoundError as e:
        print(f"Error loading data for Ring {i}: {e}")
        # Append np.nan to median lists if data is missing
        median_spr_tot_nom_ring.append(np.nan)
        median_spr_tot_ext_ring.append(np.nan)
        median_spr_tot_sec_ring.append(np.nan)
        median_spr_tot_nom_ring_error.append(np.nan)
        median_spr_tot_ext_ring_error.append(np.nan)
        median_spr_tot_sec_ring_error.append(np.nan)
        # Append zeros to efficiency lists to maintain list length
        eff_extended_ring.append(0)
        eff_secondary_ring.append(0)
        eff_nom_cob_ring.append(0)
        eff_ext_cob_ring.append(0)
        eff_sec_cob_ring.append(0)
        
        # Append zeros to uncertainties lists
        eff_extended_ring_error.append(0)
        eff_secondary_ring_error.append(0)
        eff_nom_cob_ring_error.append(0)
        eff_ext_cob_ring_error.append(0)
        eff_sec_cob_ring_error.append(0)
        
        # Append zeros to masked efficiencies and uncertainties
        eff_extended_ring_masked.append(0)
        eff_secondary_ring_masked.append(0)
        eff_nom_cob_ring_masked.append(0)
        eff_ext_cob_ring_masked.append(0)
        eff_sec_cob_ring_masked.append(0)
        
        eff_extended_ring_masked_error.append(0)
        eff_secondary_ring_masked_error.append(0)
        eff_nom_cob_ring_masked_error.append(0)
        eff_ext_cob_ring_masked_error.append(0)
        eff_sec_cob_ring_masked_error.append(0)
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
        try:
            sprk_10first = ring_nominal[j, 17:27]   # Adjust indices if necessary
            nsr_1h_24_cameras_nominal_mask = ring_nominal[j, 7]  # Adjust index if necessary
            SPR_tot = ring_nominal[j, 11]  # Adjust index if necessary
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in SPR_tot extraction: {e}")
            sprk_10first = np.zeros(10)
            nsr_1h_24_cameras_nominal_mask = 0
            SPR_tot = 0

        # Compute eta_nom_bt_24_cameras and eta_nom_bt_6_cameras
        try:
            eta_nom_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                / (ring_nominal[j, 7] * (1 - ring_nominal[j, 11]))
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_nom_bt_24_cameras calculation: {e}")
            eta_nom_bt_24_cameras[j, :] = 0

        try:
            eta_nom_bt_6_cameras[j, :] = (
                gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
                / (ring_nominal[j, 148] * (1 - ring_nominal[j, 11]))  # Adjust index if necessary
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_nom_bt_6_cameras calculation: {e}")
            eta_nom_bt_6_cameras[j, :] = 0

        # Compute eta_ext_bt_24_cameras and eta_ext_bt_6_cameras
        try:
            eta_ext_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                / (ring_extended[j, 4] * (1 - ring_extended[j, 13]))
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_ext_bt_24_cameras calculation: {e}")
            eta_ext_bt_24_cameras[j, :] = 0

        try:
            eta_ext_bt_6_cameras[j, :] = (
                gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr)
                / (ring_extended[j, 44] * (1 - ring_extended[j, 13]))  # Adjust index if necessary
            )
        except IndexError as e:
            print(f"Index error for Ring {i}, Target {j} in eta_ext_bt_6_cameras calculation: {e}")
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
        try:
            eta_cob_sec = ring_secondary[j, 9] if ring_secondary.shape[1] > 9 else 0
        except IndexError:
            eta_cob_sec = 0

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

    # --- New Code for Median SPR_tot Calculations Starts Here ---

    # Compute median SPR_tot_nom, SPR_tot_ext, and SPR_tot_sec for the current ring
    if SPR_tot_nom.size > 0:
        median_spr_nom = np.median(SPR_tot_nom)
        error_spr_nom = bootstrap_median(SPR_tot_nom)
    else:
        median_spr_nom = np.nan
        error_spr_nom = np.nan

    if SPR_tot_ext.size > 0:
        median_spr_ext = np.median(SPR_tot_ext)
        error_spr_ext = bootstrap_median(SPR_tot_ext)
    else:
        median_spr_ext = np.nan
        error_spr_ext = np.nan

    if SPR_tot_sec.size > 0:
        median_spr_sec = np.median(SPR_tot_sec)
        error_spr_sec = bootstrap_median(SPR_tot_sec)
    else:
        median_spr_sec = np.nan
        error_spr_sec = np.nan

    # Append the median values and their errors to the respective lists for all targets
    median_spr_tot_nom_ring.append(median_spr_nom)
    median_spr_tot_nom_ring_error.append(error_spr_nom)
    median_spr_tot_ext_ring.append(median_spr_ext)
    median_spr_tot_ext_ring_error.append(error_spr_ext)
    median_spr_tot_sec_ring.append(median_spr_sec)
    median_spr_tot_sec_ring_error.append(error_spr_sec)

    # Similarly, compute medians and errors for masked targets
    if SPR_tot_nom_masked.size > 0:
        median_spr_nom_masked = np.median(SPR_tot_nom_masked)
        error_spr_nom_masked = bootstrap_median(SPR_tot_nom_masked)
    else:
        median_spr_nom_masked = np.nan
        error_spr_nom_masked = np.nan

    if SPR_tot_ext_masked.size > 0:
        median_spr_ext_masked = np.median(SPR_tot_ext_masked)
        error_spr_ext_masked = bootstrap_median(SPR_tot_ext_masked)
    else:
        median_spr_ext_masked = np.nan
        error_spr_ext_masked = np.nan

    if SPR_tot_sec_masked.size > 0:
        median_spr_sec_masked = np.median(SPR_tot_sec_masked)
        error_spr_sec_masked = bootstrap_median(SPR_tot_sec_masked)
    else:
        median_spr_sec_masked = np.nan
        error_spr_sec_masked = np.nan

    # Append the median values and their errors to the respective lists for masked targets
    median_spr_tot_nom_masked_ring.append(median_spr_nom_masked)
    median_spr_tot_nom_masked_ring_error.append(error_spr_nom_masked)
    median_spr_tot_ext_masked_ring.append(median_spr_ext_masked)
    median_spr_tot_ext_masked_ring_error.append(error_spr_ext_masked)
    median_spr_tot_sec_masked_ring.append(median_spr_sec_masked)
    median_spr_tot_sec_masked_ring_error.append(error_spr_sec_masked)

    # --- New Code for Median SPR_tot Calculations Ends Here ---
    # ----- Plot 1: Median SPR_tot_nom, SPR_tot_ext, and SPR_tot_sec vs Ring Number for All Targets with Error Bars -----
plt.figure(figsize=(10, 6), dpi=120)

plt.errorbar(
    ring_numbers,
    median_spr_tot_nom_ring,
    yerr=median_spr_tot_nom_ring_error,
    fmt='o-', 
    color='brown',
    label=r'Median $\rm SPR_{tot\_nom}$ (All Targets)',
    capsize=5,
    ecolor='black',
    alpha=0.7
)
plt.errorbar(
    ring_numbers,
    median_spr_tot_ext_ring,
    yerr=median_spr_tot_ext_ring_error,
    fmt='s-', 
    color='blue',
    label=r'Median $\rm SPR_{tot\_ext}$ (All Targets)',
    capsize=5,
    ecolor='black',
    alpha=0.7
)
plt.errorbar(
    ring_numbers,
    median_spr_tot_sec_ring,
    yerr=median_spr_tot_sec_ring_error,
    fmt='^-', 
    color='green',
    label=r'Median $\rm SPR_{tot\_sec}$ (All Targets)',
    capsize=5,
    ecolor='black',
    alpha=0.7
)

plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'Median $\rm SPR_{tot}$', fontsize=14)
plt.title("Median SPR_tot as a Function of Ring Number (All Targets)", fontsize=16)
plt.xticks(ring_numbers, [f"Ring {i}" for i in ring_numbers], fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("Median_SPR_tot_vs_Ring_All_Targets_with_Errors.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ----- Plot 2: Median SPR_tot_nom, SPR_tot_ext, and SPR_tot_sec vs Ring Number for Masked Targets with Error Bars -----
plt.figure(figsize=(10, 6), dpi=120)

plt.errorbar(
    ring_numbers,
    median_spr_tot_nom_masked,
    yerr=error_spr_nom_masked,
    fmt='o-', 
    color='brown',
    label=r'Median $\rm SPR_{tot\_nom}$ (Masked Targets)',
    capsize=5,
    ecolor='black',
    alpha=0.7
)
plt.errorbar(
    ring_numbers,
    median_spr_tot_ext_masked,
    yerr=error_spr_ext_masked,
    fmt='s-', 
    color='red',
    label=r'Median $\rm SPR_{tot\_ext}$ (Masked Targets)',
    capsize=5,
    ecolor='black',
    alpha=0.7
)


plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'Median $\rm SPR_{tot}$', fontsize=14)
plt.title("Median SPR_tot as a Function of Ring Number (Masked Targets)", fontsize=16)
plt.xticks(ring_numbers, [f"Ring {i}" for i in ring_numbers], fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("Median_SPR_tot_vs_Ring_Masked_Targets_with_Errors.pdf", format='pdf', bbox_inches='tight')
plt.show()