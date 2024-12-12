"""
Script for the efficiency metrics
for different quadrants in the FP 
Fernando 28th Nov. 2024
"""
import numpy as np  # type:ignore
import matplotlib.pyplot as plt  # type:ignore

dataDIR_quadrants = '/home/fercho/double-aperture-photometry/simulation_results/rings/quarters/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'

# Initialize lists to store efficiencies for each quadrant (all targets)
eff_extended_quadrant = []
eff_secondary_quadrant = []
eff_nom_cob_quadrant = []
eff_ext_cob_quadrant = []
eff_sec_cob_quadrant = []

# Initialize lists to store uncertainties for each quadrant (all targets)
eff_extended_quadrant_error = []
eff_secondary_quadrant_error = []
eff_nom_cob_quadrant_error = []
eff_ext_cob_quadrant_error = []
eff_sec_cob_quadrant_error = []

# Initialize lists to store efficiencies for masked targets
eff_extended_quadrant_masked = []
eff_secondary_quadrant_masked = []
eff_nom_cob_quadrant_masked = []
eff_ext_cob_quadrant_masked = []
eff_sec_cob_quadrant_masked = []

# Initialize lists to store uncertainties for masked targets
eff_extended_quadrant_masked_error = []
eff_secondary_quadrant_masked_error = []
eff_nom_cob_quadrant_masked_error = []
eff_ext_cob_quadrant_masked_error = []
eff_sec_cob_quadrant_masked_error = []

quadrant_names = ['Q1', 'Q2', 'Q3', 'Q4']  # List of quadrants
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Reference transit depth (ppm)
td_ref = 6.72 * 0.46**2  # Transit duration in hours
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other value could be 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for masking targets
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

# Loop through each quadrant, load corresponding .npy files, and compute metrics
for idx, quadrant_name in enumerate(quadrant_names):
    # Load .npy files for the current quadrant
    try:
        quadrant_nominal = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_nominal.npy")
        quadrant_secondary = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_secondary.npy")
        quadrant_extended = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_extended.npy")
    except FileNotFoundError as e:
        print(f"Error loading data for {quadrant_name}: {e}")
        continue  # Skip to the next quadrant if files are missing

    # Number of targets in this quadrant
    n_targets_quadrant = quadrant_nominal.shape[0]

    # Initialize arrays for this quadrant
    eta_ext_bt_24_cameras = np.zeros((n_targets_quadrant, 10))
    eta_nom_bt_24_cameras = np.zeros((n_targets_quadrant, 10))
    delta_obs = np.zeros((n_targets_quadrant, 10))
    delta_obs_ext = np.zeros((n_targets_quadrant, 10))
    sig_depth_24_cameras = np.zeros((n_targets_quadrant, 10))
    sig_depth_sec_nom_quad = np.zeros((n_targets_quadrant, 10))

    # Initialize boolean arrays for conditions
    nfp = np.zeros((n_targets_quadrant, 10), dtype=bool)
    nfp_ext_mask = np.zeros((n_targets_quadrant, 10), dtype=bool)
    nfp_nom_cob = np.zeros((n_targets_quadrant, 10), dtype=bool)
    nfp_ext_cob = np.zeros((n_targets_quadrant, 10), dtype=bool)
    fp_single_contaminant_24_cameras = np.zeros((n_targets_quadrant, 10), dtype=bool)
    secondary_mask_conditions_24_cameras = np.zeros((n_targets_quadrant, 10), dtype=bool)
    secondary_mask_cob_conditions_24_cameras = np.zeros((n_targets_quadrant, 10), dtype=bool)

    # Array to track masked targets
    masked_targets = np.zeros(n_targets_quadrant, dtype=bool)

    # Compute metrics for each target in the quadrant
    for j in range(n_targets_quadrant):
        dback = np.ones(10) * dback_ref  # Reference transit depth
        td = np.ones(10) * td_ref      # Reference transit duration

        # Flux significance for nominal mask
        try:
            eta_nom_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * quadrant_nominal[j, 17:27] * np.sqrt(td * ntr) /
                (quadrant_nominal[j, 7] * (1 - quadrant_nominal[j, 11]))
            )
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in eta_nom calculations: {e}")
            eta_nom_bt_24_cameras[j, :] = 0

        # Flux significance for extended mask
        try:
            eta_ext_bt_24_cameras[j, :] = (
                gamma_factor_significance * dback * quadrant_extended[j, 14:24] * np.sqrt(td * ntr) /
                (quadrant_extended[j, 4] * (1 - quadrant_extended[j, 13]))
            )
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in eta_ext calculations: {e}")
            eta_ext_bt_24_cameras[j, :] = 0

        # Observed transit depth
        try:
            delta_obs[j, :] = dback * quadrant_nominal[j, 17:27]
            delta_obs_ext[j, :] = dback * quadrant_extended[j, 14:24]
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in observed transit depth extraction: {e}")
            delta_obs[j, :] = 0
            delta_obs_ext[j, :] = 0

        # Compute significant transit depth variables
        try:
            # NSR1h for secondary mask (24 cameras)
            nsr1h_sec = quadrant_secondary[j, 4]
            # NSR1h for nominal mask (24 cameras)
            nsr1h_nom = quadrant_nominal[j, 7]
            # SPRtot for secondary mask
            SPR_tot_sec = quadrant_secondary[j, 5]
            # SPRtot for extended mask
            SPR_tot_ext = quadrant_extended[j, 13]
            # SPRtot for nominal mask
            SPR_tot_nom = quadrant_nominal[j, 11]
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in SPR_tot extraction: {e}")
            SPR_tot_sec = 0
            SPR_tot_ext = 0
            SPR_tot_nom = 0

        # Significant transit depth calculations
        try:
            sig_depth_secondary_mask_24_cameras = nsr1h_sec * (1 - SPR_tot_sec) / np.sqrt(td_ref * ntr)
            sig_depth_extended_mask_24_cameras = quadrant_extended[j, 4] * (1 - SPR_tot_ext) / np.sqrt(td_ref * ntr)
            sig_depth_nominal_mask_24_cameras = nsr1h_nom * (1 - SPR_tot_nom) / np.sqrt(td_ref * ntr)
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in sig_depth calculations: {e}")
            sig_depth_secondary_mask_24_cameras = 0
            sig_depth_extended_mask_24_cameras = 0
            sig_depth_nominal_mask_24_cameras = 0

        # Quadratic sum of the noises (Equation (38) from the paper)
        sig_depth_24_cameras[j, :] = np.sqrt(
            sig_depth_nominal_mask_24_cameras ** 2 + sig_depth_extended_mask_24_cameras ** 2
        )
        sig_depth_sec_nom_quad[j, :] = np.sqrt(
            sig_depth_secondary_mask_24_cameras ** 2 + sig_depth_nominal_mask_24_cameras ** 2
        )

        # COB-related variables
        try:
            eta_cob_nom_10first_24_cameras = quadrant_nominal[j, 46:56]
            eta_cob_ext_10first_24_cameras = quadrant_extended[j, 45:55]
            eta_cob_sec = quadrant_secondary[j, 9]
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in COB extraction: {e}")
            eta_cob_nom_10first_24_cameras = np.zeros(10)
            eta_cob_ext_10first_24_cameras = np.zeros(10)
            eta_cob_sec = 0

        # False Positive detection conditions
        nfp[j, :] = eta_nom_bt_24_cameras[j, :] > flux_thresh_nom_mask

        # Extended mask detection conditions
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
            eta_c = quadrant_secondary[j, 6]
            delta_obs_c = quadrant_secondary[j, 7]
            delta_obs_t = delta_obs[j, :]
        except IndexError as e:
            print(f"Index error for {quadrant_name}, Target {j} in secondary mask conditions: {e}")
            eta_c = 0
            delta_obs_c = 0
            delta_obs_t = np.zeros(10)

        secondary_mask_conditions_24_cameras[j, :] = (
            (eta_c > flux_thresh_sec_mask) &
            (delta_obs_c > delta_obs_t + depth_sig_scaling * sig_depth_sec_nom_quad[j, :]) &
            fp_single_contaminant_24_cameras[j, :]
        )
        secondary_mask_cob_conditions_24_cameras[j, :] = (
            (eta_cob_sec > cob_thresh) & fp_single_contaminant_24_cameras[j, :]
        )

        # Identify if the target is masked (eta_nom > 7.1 in any of the 10 cameras)
        if np.any(eta_nom_bt_24_cameras[j, :] > eta_nom_threshold):
            masked_targets[j] = True

    # Compute efficiencies for all targets
    nfp_total = nfp.sum()
    if nfp_total > 0:
        eff_ext_flux = (nfp & nfp_ext_mask).sum() / nfp_total * 100.0
        eff_nom_cob = (nfp & nfp_nom_cob).sum() / nfp_total * 100.0
        eff_ext_cob = (nfp & nfp_ext_cob).sum() / nfp_total * 100.0
    else:
        eff_ext_flux = 0
        eff_nom_cob = 0
        eff_ext_cob = 0

    # Secondary mask efficiency for all targets
    fp_single_total = fp_single_contaminant_24_cameras.sum()
    if fp_single_total > 0:
        eff_secondary = secondary_mask_conditions_24_cameras.sum() / fp_single_total * 100.0
        eff_secondary_cob = secondary_mask_cob_conditions_24_cameras.sum() / fp_single_total * 100.0
    else:
        eff_secondary = 0
        eff_secondary_cob = 0

    # Store efficiencies for each quadrant (all targets)
    eff_extended_quadrant.append(eff_ext_flux)
    eff_secondary_quadrant.append(eff_secondary)
    eff_nom_cob_quadrant.append(eff_nom_cob)
    eff_ext_cob_quadrant.append(eff_ext_cob)
    eff_sec_cob_quadrant.append(eff_secondary_cob)  # Corrected from eff_ext_cob

    # Compute uncertainties for each efficiency metric (all targets)
    # Using binomial error: sqrt(p*(1-p)/n) * 100
    # Where p = k/n, k is number of successes, n is number of trials

    # Number of successes for each efficiency metric
    k_nom_cob = (nfp & nfp_nom_cob).sum()
    k_ext_cob = (nfp & nfp_ext_cob).sum()
    k_sec_cob = eff_secondary_cob / 100.0 * fp_single_total if fp_single_total > 0 else 0
    k_secondary = eff_secondary / 100.0 * fp_single_total if fp_single_total > 0 else 0
    k_ext_flux = (nfp & nfp_ext_mask).sum()

    # Number of trials for each efficiency metric
    n_nom_cob = nfp_total
    n_ext_cob = nfp_total
    n_sec_cob = fp_single_total
    n_secondary = fp_single_total
    n_ext_flux = nfp_total

    # Calculate uncertainties
    delta_nom_cob = binomial_error(int(k_nom_cob), int(n_nom_cob))
    delta_ext_cob = binomial_error(int(k_ext_cob), int(n_ext_cob))
    delta_sec_cob = binomial_error(int(k_sec_cob), int(n_sec_cob))
    delta_secondary = binomial_error(int(k_secondary), int(n_secondary))
    delta_ext_flux = binomial_error(int(k_ext_flux), int(n_ext_flux))

    # Store uncertainties for each quadrant (all targets)
    eff_nom_cob_quadrant_error.append(delta_nom_cob)
    eff_ext_cob_quadrant_error.append(delta_ext_cob)
    eff_sec_cob_quadrant_error.append(delta_sec_cob)
    eff_secondary_quadrant_error.append(delta_secondary)
    eff_extended_quadrant_error.append(delta_ext_flux)

    # --- Compute efficiencies for masked targets ---
    if np.any(masked_targets):
        # Extract masked targets' data
        nfp_masked = nfp[masked_targets, :]
        nfp_ext_mask_masked = nfp_ext_mask[masked_targets, :]
        nfp_nom_cob_masked = nfp_nom_cob[masked_targets, :]
        nfp_ext_cob_masked = nfp_ext_cob[masked_targets, :]
        fp_single_contaminant_masked = fp_single_contaminant_24_cameras[masked_targets, :]
        secondary_mask_conditions_masked = secondary_mask_conditions_24_cameras[masked_targets, :]
        secondary_mask_cob_conditions_masked = secondary_mask_cob_conditions_24_cameras[masked_targets, :]

        # Total number of False Positives (FP) for masked targets
        nfp_total_masked = nfp_masked.sum()

        if nfp_total_masked > 0:
            eff_ext_flux_masked = (nfp_masked & nfp_ext_mask_masked).sum() / nfp_total_masked * 100.0
            eff_nom_cob_masked = (nfp_masked & nfp_nom_cob_masked).sum() / nfp_total_masked * 100.0
            eff_ext_cob_masked = (nfp_masked & nfp_ext_cob_masked).sum() / nfp_total_masked * 100.0
        else:
            eff_ext_flux_masked = 0
            eff_nom_cob_masked = 0
            eff_ext_cob_masked = 0

        # Secondary mask efficiency for masked targets
        fp_single_total_masked = fp_single_contaminant_masked.sum()
        if fp_single_total_masked > 0:
            eff_secondary_masked = secondary_mask_conditions_masked.sum() / fp_single_total_masked * 100.0
            eff_sec_cob_masked = secondary_mask_cob_conditions_masked.sum() / fp_single_total_masked * 100.0
        else:
            eff_secondary_masked = 0
            eff_sec_cob_masked = 0

        # Store efficiencies for masked targets
        eff_extended_quadrant_masked.append(eff_ext_flux_masked)
        eff_secondary_quadrant_masked.append(eff_secondary_masked)
        eff_nom_cob_quadrant_masked.append(eff_nom_cob_masked)
        eff_ext_cob_quadrant_masked.append(eff_ext_cob_masked)
        eff_sec_cob_quadrant_masked.append(eff_sec_cob_masked)

        # Compute uncertainties for masked targets
        # Number of successes for each efficiency metric
        k_nom_cob_masked = (nfp_masked & nfp_nom_cob_masked).sum()
        k_ext_cob_masked = (nfp_masked & nfp_ext_cob_masked).sum()
        k_sec_cob_masked = eff_sec_cob_masked / 100.0 * fp_single_total_masked if fp_single_total_masked > 0 else 0
        k_secondary_masked = eff_secondary_masked / 100.0 * fp_single_total_masked if fp_single_total_masked > 0 else 0
        k_ext_flux_masked = (nfp_masked & nfp_ext_mask_masked).sum()

        # Number of trials for each efficiency metric
        n_nom_cob_masked = nfp_total_masked
        n_ext_cob_masked = nfp_total_masked
        n_sec_cob_masked = fp_single_total_masked
        n_secondary_masked = fp_single_total_masked
        n_ext_flux_masked = nfp_total_masked

        # Calculate uncertainties
        delta_nom_cob_masked = binomial_error(int(k_nom_cob_masked), int(n_nom_cob_masked))
        delta_ext_cob_masked = binomial_error(int(k_ext_cob_masked), int(n_ext_cob_masked))
        delta_sec_cob_masked = binomial_error(int(k_sec_cob_masked), int(n_sec_cob_masked))
        delta_secondary_masked = binomial_error(int(k_secondary_masked), int(n_secondary_masked))
        delta_ext_flux_masked = binomial_error(int(k_ext_flux_masked), int(n_ext_flux_masked))

        # Store uncertainties for masked targets
        eff_nom_cob_quadrant_masked_error.append(delta_nom_cob_masked)
        eff_ext_cob_quadrant_masked_error.append(delta_ext_cob_masked)
        eff_sec_cob_quadrant_masked_error.append(delta_sec_cob_masked)
        eff_secondary_quadrant_masked_error.append(delta_secondary_masked)
        eff_extended_quadrant_masked_error.append(delta_ext_flux_masked)
    else:
        # If no masked targets, append zeros
        eff_extended_quadrant_masked.append(0)
        eff_secondary_quadrant_masked.append(0)
        eff_nom_cob_quadrant_masked.append(0)
        eff_ext_cob_quadrant_masked.append(0)
        eff_sec_cob_quadrant_masked.append(0)

        eff_extended_quadrant_masked_error.append(0)
        eff_secondary_quadrant_masked_error.append(0)
        eff_nom_cob_quadrant_masked_error.append(0)
        eff_ext_cob_quadrant_masked_error.append(0)
        eff_sec_cob_quadrant_masked_error.append(0)

    # Print efficiencies for the current quadrant with error bars
    print(f"Quadrant {quadrant_name} Efficiencies:")
    print(f"  Nominal COB Efficiency: {eff_nom_cob:.2f}% ± {delta_nom_cob:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob:.2f}% ± {delta_ext_cob:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_secondary_cob:.2f}% ± {delta_sec_cob:.2f}%")
    print(f"  Secondary Flux Efficiency: {eff_secondary:.2f}% ± {delta_secondary:.2f}%")
    print(f"  Extended Flux Efficiency: {eff_ext_flux:.2f}% ± {delta_ext_flux:.2f}%\n")

    # If there are masked targets, print their efficiencies as well
    if np.any(masked_targets):
        print(f"Quadrant {quadrant_name} Masked Targets Efficiencies:")
        print(f"  Nominal COB Efficiency (Masked): {eff_nom_cob_masked:.2f}% ± {delta_nom_cob_masked:.2f}%")
        print(f"  Extended COB Efficiency (Masked): {eff_ext_cob_masked:.2f}% ± {delta_ext_cob_masked:.2f}%")
        print(f"  Secondary COB Efficiency (Masked): {eff_sec_cob_masked:.2f}% ± {delta_sec_cob_masked:.2f}%")
        print(f"  Secondary Flux Efficiency (Masked): {eff_secondary_masked:.2f}% ± {delta_secondary_masked:.2f}%")
        print(f"  Extended Flux Efficiency (Masked): {eff_ext_flux_masked:.2f}% ± {delta_ext_flux_masked:.2f}%\n")

# --- Plotting Efficiencies for All Targets ---
fig_all, ax_all = plt.subplots(figsize=(12, 8))

# Define positions for each bar group
x = np.arange(len(quadrant_names))  # the label locations
width = 0.15  # the width of the bars

# Plot each efficiency metric with corresponding error bars
rects1 = ax_all.bar(x - 2*width, eff_nom_cob_quadrant, width, yerr=eff_nom_cob_quadrant_error,
                    label='Nominal COB Efficiency', color='brown', capsize=5)
rects2 = ax_all.bar(x - width, eff_ext_cob_quadrant, width, yerr=eff_ext_cob_quadrant_error,
                    label='Extended COB Efficiency', color='blue', capsize=5)
rects3 = ax_all.bar(x, eff_sec_cob_quadrant, width, yerr=eff_sec_cob_quadrant_error,
                    label='Secondary COB Efficiency', color='purple', capsize=5)
rects4 = ax_all.bar(x + width, eff_secondary_quadrant, width, yerr=eff_secondary_quadrant_error,
                    label='Secondary Flux Efficiency', color='green', capsize=5)
rects5 = ax_all.bar(x + 2*width, eff_extended_quadrant, width, yerr=eff_extended_quadrant_error,
                    label='Extended Flux Efficiency', color='red', capsize=5)

# Add labels, title, and custom x-axis tick labels
ax_all.set_xlabel('Quadrant', fontsize=14)
ax_all.set_ylabel('Efficiency (%)', fontsize=14)
ax_all.set_title('Efficiencies Across Quadrants with Error Bars (All Targets)', fontsize=16)
ax_all.set_xticks(x)
ax_all.set_xticklabels(quadrant_names, fontsize=12)
ax_all.legend(fontsize=12)
ax_all.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("Efficiencies_Across_Quadrants_All_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()

# --- Plotting Efficiencies for Masked Targets ---
fig_masked, ax_masked = plt.subplots(figsize=(12, 8))

# Plot each efficiency metric with corresponding error bars
rects1_masked = ax_masked.bar(x - 2*width, eff_nom_cob_quadrant_masked, width, yerr=eff_nom_cob_quadrant_masked_error,
                                label='Nominal COB Efficiency (Masked)', color='brown', capsize=5)
rects2_masked = ax_masked.bar(x - width, eff_ext_cob_quadrant_masked, width, yerr=eff_ext_cob_quadrant_masked_error,
                                label='Extended COB Efficiency (Masked)', color='blue', capsize=5)
rects3_masked = ax_masked.bar(x, eff_sec_cob_quadrant_masked, width, yerr=eff_sec_cob_quadrant_masked_error,
                                label='Secondary COB Efficiency (Masked)', color='purple', capsize=5)
rects4_masked = ax_masked.bar(x + width, eff_secondary_quadrant_masked, width, yerr=eff_secondary_quadrant_masked_error,
                                label='Secondary Flux Efficiency (Masked)', color='green', capsize=5)
rects5_masked = ax_masked.bar(x + 2*width, eff_extended_quadrant_masked, width, yerr=eff_extended_quadrant_masked_error,
                                label='Extended Flux Efficiency (Masked)', color='red', capsize=5)

# Add labels, title, and custom x-axis tick labels
ax_masked.set_xlabel('Quadrant', fontsize=14)
ax_masked.set_ylabel('Efficiency (%)', fontsize=14)
ax_masked.set_title('Efficiencies Across Quadrants with Error Bars (Masked Targets)', fontsize=16)
ax_masked.set_xticks(x)
ax_masked.set_xticklabels(quadrant_names, fontsize=12)
ax_masked.legend(fontsize=12)
ax_masked.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("Efficiencies_Across_Quadrants_Masked_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()