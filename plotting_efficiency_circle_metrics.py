import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

# Some file parameters
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/rings/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'
eta_ext_bt_24_cameras_ring = []
eta_ext_bt_6_cameras_ring = []
eta_nom_bt_24_cameras_ring = []
eta_nom_bt_6_cameras_ring = []
delta_obs_ring = []
delta_obs_ext_ring = []
delta_obs_ext_6_cameras_ring = []

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

number_of_circles = 7  # Adjust as needed
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Reference value for transit depth (ppm)
td_ref = 6.72 * 0.46**2  # Transit duration in hours
seed = 123434434
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other possible value: 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Threshold for eta_nom_bt_24_cameras
eta_nom_threshold = 7.1

# Loop through each ring, load corresponding .npy files, and compute metrics
for i in range(number_of_circles):
    # Load .npy files for the current ring
    ring_nominal = np.load(f"{dataDIR}ring_{i}_nominal.npy")
    ring_secondary = np.load(f"{dataDIR}ring_{i}_secondary.npy")
    ring_extended = np.load(f"{dataDIR}ring_{i}_extended.npy")

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

        # Flux significance for nominal mask
        eta_nom_bt_24_cameras[j, :] = gamma_factor_significance * dback * ring_nominal[j, 17:27] * np.sqrt(td * ntr) / \
                                      (ring_nominal[j, 7] * (1 - ring_nominal[j, 11]))
        eta_nom_bt_6_cameras[j, :] = gamma_factor_significance * dback * ring_nominal[j, 17:27] * np.sqrt(td * ntr) / \
                                     (ring_nominal[j, 148] * (1 - ring_nominal[j, 11]))

        # Flux significance for extended mask
        eta_ext_bt_24_cameras[j, :] = gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr) / \
                                      (ring_extended[j, 4] * (1 - ring_extended[j, 13]))
        eta_ext_bt_6_cameras[j, :] = gamma_factor_significance * dback * ring_extended[j, 14:24] * np.sqrt(td * ntr) / \
                                     (ring_extended[j, 44] * (1 - ring_extended[j, 13]))

        # Observed transit depth
        delta_obs[j, :] = dback * ring_nominal[j, 17:27]
        delta_obs_ext[j, :] = dback * ring_extended[j, 14:24]

        # Compute significant transit depth variables
        # NSR1h for secondary mask (24 cameras)
        nsr1h_sec = ring_secondary[j, 4]
        # NSR1h for nominal mask (24 cameras)
        nsr1h_nom = ring_nominal[j, 7]
        # SPRtot for secondary mask
        SPR_tot_sec = ring_secondary[j, 5]
        # SPRtot for extended mask
        SPR_tot_ext = ring_extended[j, 13]
        # SPRtot for nominal mask
        SPR_tot_nom = ring_nominal[j, 11]

        # Significant transit depth calculations
        sig_depth_secondary_mask_24_cameras = nsr1h_sec * (1 - SPR_tot_sec) / np.sqrt(td_ref * ntr)
        sig_depth_extended_mask_24_cameras = ring_extended[j, 4] * (1 - SPR_tot_ext) / np.sqrt(td_ref * ntr)
        sig_depth_nominal_mask_24_cameras = nsr1h_nom * (1 - SPR_tot_nom) / np.sqrt(td_ref * ntr)

        # Quadratic sum of the noises (Equation (38) from the paper)
        sig_depth_24_cameras[j, :] = np.sqrt(sig_depth_nominal_mask_24_cameras ** 2 + sig_depth_extended_mask_24_cameras ** 2)
        sig_depth_sec_nom_quad[j, :] = np.sqrt(sig_depth_secondary_mask_24_cameras ** 2 + sig_depth_nominal_mask_24_cameras ** 2)

        # COB-related variables
        eta_cob_nom_10first_24_cameras = ring_nominal[j, 46:56]
        eta_cob_ext_10first_24_cameras = ring_extended[j, 45:55]
        eta_cob_sec = ring_secondary[j, 9]

        # False Positive detection conditions
        nfp[j, :] = eta_nom_bt_24_cameras[j, :] > flux_thresh_nom_mask

        # Extended mask detection conditions
        nfp_ext_mask[j, :] = (eta_ext_bt_24_cameras[j, :] > flux_thresh_ext_mask) & \
                             (delta_obs_ext[j, :] > delta_obs[j, :] + depth_sig_scaling * sig_depth_24_cameras[j, :])

        # COB detection conditions
        nfp_nom_cob[j, :] = eta_cob_nom_10first_24_cameras > cob_thresh
        nfp_ext_cob[j, :] = eta_cob_ext_10first_24_cameras > cob_thresh

        # Single contaminant false positives
        fp_single_contaminant_24_cameras[j, :] = nfp[j, :]

        # Secondary mask conditions
        eta_c = ring_secondary[j, 6]
        delta_obs_c = ring_secondary[j, 7]
        delta_obs_t = delta_obs[j, :]
        secondary_mask_conditions_24_cameras[j, :] = (
            (eta_c > flux_thresh_sec_mask) &
            (delta_obs_c > delta_obs_t + depth_sig_scaling * sig_depth_sec_nom_quad[j, :]) &
            fp_single_contaminant_24_cameras[j, :]
        )

        # Secondary mask COB conditions
        eta_cob_sec = ring_secondary[j, 9]
        secondary_mask_cob_conditions_24_cameras[j, :] = (eta_cob_sec > cob_thresh) & fp_single_contaminant_24_cameras[j, :]

    # --- New Code for SPR_tot Accumulation Starts Here ---

    # Extract SPR_tot_nom, SPR_tot_ext, and SPR_tot_sec for all targets in this ring
    SPR_tot_nom = ring_nominal[:, 11]   # Adjust index if necessary
    SPR_tot_ext = ring_extended[:, 13]  # Adjust index if necessary
    SPR_tot_sec = ring_secondary[:, 9]  # Adjust index if necessary

    # Append to all targets lists
    spr_tot_nom_all.append(SPR_tot_nom)
    spr_tot_ext_all.append(SPR_tot_ext)
    spr_tot_sec_all.append(SPR_tot_sec)

    # Extract SPR_tot for masked targets
    eta_nom_above_threshold = np.any(eta_nom_bt_24_cameras > eta_nom_threshold, axis=1)  # Shape (n_targets_ring,)
    SPR_tot_nom_masked = SPR_tot_nom[eta_nom_above_threshold]
    SPR_tot_ext_masked = SPR_tot_ext[eta_nom_above_threshold]
    SPR_tot_sec_masked = SPR_tot_sec[eta_nom_above_threshold]

    # Append to masked targets lists
    spr_tot_nom_masked.append(SPR_tot_nom_masked)
    spr_tot_ext_masked.append(SPR_tot_ext_masked)
    spr_tot_sec_masked.append(SPR_tot_sec_masked)


    # Create a boolean mask for targets where any eta_nom_bt_24_cameras > threshold
    eta_nom_above_threshold = np.any(eta_nom_bt_24_cameras > eta_nom_threshold, axis=1)  # Shape (n_targets_ring,)

    # Compute efficiencies for all targets (as before)

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

    # Now compute efficiencies only for targets where eta_nom_bt_24_cameras > threshold

    # Create masked versions of the arrays based on eta_nom_above_threshold
    nfp_masked = nfp[eta_nom_above_threshold, :]
    nfp_ext_mask_masked = nfp_ext_mask[eta_nom_above_threshold, :]
    nfp_nom_cob_masked = nfp_nom_cob[eta_nom_above_threshold, :]
    nfp_ext_cob_masked = nfp_ext_cob[eta_nom_above_threshold, :]
    fp_single_contaminant_masked = fp_single_contaminant_24_cameras[eta_nom_above_threshold, :]
    secondary_mask_conditions_masked = secondary_mask_conditions_24_cameras[eta_nom_above_threshold, :]
    secondary_mask_cob_conditions_masked = secondary_mask_cob_conditions_24_cameras[eta_nom_above_threshold, :]

    # Total number of False Positives (FP) for masked targets
    nfp_total_masked = nfp_masked.sum()

    if nfp_total_masked > 0:
        eff_ext_flux_masked = (nfp_masked & nfp_ext_mask_masked).sum() / nfp_total_masked * 100.
        eff_nom_cob_masked = (nfp_masked & nfp_nom_cob_masked).sum() / nfp_total_masked * 100.
        eff_ext_cob_masked = (nfp_masked & nfp_ext_cob_masked).sum() / nfp_total_masked * 100.
    else:
        eff_ext_flux_masked = 0
        eff_nom_cob_masked = 0
        eff_ext_cob_masked = 0

    # Secondary mask efficiency for masked targets
    fp_single_total_masked = fp_single_contaminant_masked.sum()
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

    # Print efficiencies for the current ring
    print(f"Ring {i} Efficiencies for All Targets:")
    print(f"  Extended Flux Efficiency: {eff_ext_flux_all:.2f}%")
    print(f"  Secondary Mask Efficiency: {eff_secondary_all:.2f}%")
    print(f"  Nominal COB Efficiency: {eff_nom_cob_all:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob_all:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_sec_cob_all:.2f}%")

    print(f"Ring {i} Efficiencies for Targets with eta_nom_bt_24_cameras > {eta_nom_threshold}:")
    print(f"  Extended Flux Efficiency: {eff_ext_flux_masked:.2f}%")
    print(f"  Secondary Mask Efficiency: {eff_secondary_masked:.2f}%")
    print(f"  Nominal COB Efficiency: {eff_nom_cob_masked:.2f}%")
    print(f"  Extended COB Efficiency: {eff_ext_cob_masked:.2f}%")
    print(f"  Secondary COB Efficiency: {eff_sec_cob_masked:.2f}%")
    print("-----------------------------------------------------")

# Plotting efficiencies across rings with 'o-' format
# For All Targets
plt.figure(figsize=(10, 6))
plt.plot(range(number_of_circles), eff_nom_cob_ring, 'o-', color='brown', label="Nominal COB Efficiency (All Targets)")
plt.plot(range(number_of_circles), eff_ext_cob_ring, 's-', color='blue', label="Extended COB Efficiency (All Targets)")
plt.plot(range(number_of_circles), eff_sec_cob_ring, 's-', color='purple', label="Secondary COB Efficiency (All Targets)")
plt.plot(range(number_of_circles), eff_secondary_ring, '^-', color='green', label="Secondary Flux Efficiency (All Targets)")
plt.plot(range(number_of_circles), eff_extended_ring, 'd-', color='red', label="Extended Flux Efficiency (All Targets)")
# Adjust x-axis labels
circle_labels = [f"Ring {i}" for i in range(1, number_of_circles + 1)]  # Circle 1 to Circle 7
plt.xticks(ticks=range(number_of_circles), labels=circle_labels, fontsize=10)

plt.ylabel("Efficiency (%)", fontsize=14)
plt.title("Efficiencies Across Rings (All Targets)")
plt.legend()
plt.grid(True)
plt.show()

# Plotting efficiencies across rings with 'o-' format
# For Targets with eta_nom_bt_24_cameras > threshold
plt.figure(figsize=(10, 6))
plt.plot(range(number_of_circles), eff_nom_cob_ring_masked, 'o-', color='brown', label=f"Nominal COB Efficiency (eta_nom > {eta_nom_threshold})")
plt.plot(range(number_of_circles), eff_ext_cob_ring_masked, 's-', color='blue', label=f"Extended COB Efficiency (eta_nom > {eta_nom_threshold})")
plt.plot(range(number_of_circles), eff_sec_cob_ring_masked, 's-', color='purple', label=f"Secondary COB Efficiency (eta_nom > {eta_nom_threshold})")
plt.plot(range(number_of_circles), eff_secondary_ring_masked, '^-', color='green', label=f"Secondary Flux Efficiency (eta_nom > {eta_nom_threshold})")
plt.plot(range(number_of_circles), eff_extended_ring_masked, 'd-', color='red', label=f"Extended Flux Efficiency (eta_nom > {eta_nom_threshold})")
# Adjust x-axis labels
plt.xticks(ticks=range(number_of_circles), labels=circle_labels, fontsize=10)

plt.ylabel("Efficiency (%)", fontsize=14)
plt.title(f"Efficiencies Across Rings (eta_nom_bt_24_cameras > {eta_nom_threshold})")
plt.legend()
plt.grid(True)
plt.show()

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


plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'$\rm SPR_{tot}$', fontsize=14)
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

plt.xlabel("Ring Number", fontsize=14)
plt.ylabel(r'$\rm SPR_{tot}$', fontsize=14)
plt.xticks(ring_numbers, [f"Ring {i}" for i in ring_numbers], fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)

plt.tight_layout()
plt.savefig("Median_SPR_tot_vs_Ring_Masked_Targets.pdf", format='pdf', bbox_inches='tight')
plt.show()