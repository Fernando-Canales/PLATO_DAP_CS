"""
Script for the efficiency metrics
for different quadrants in the FP 
Fernando 28th Nov. 2024
"""
import numpy as np #type:ignore
import matplotlib.pyplot as plt # type:ignore

dataDIR_quadrants = '/home/fercho/double-aperture-photometry/simulation_results/rings/quarters/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'


# Initialize lists to store efficiencies for each quadrant
eff_extended_quadrant = []
eff_secondary_quadrant = []
eff_nom_cob_quadrant = []
eff_ext_cob_quadrant = []

quadrant_names = ['Q1', 'Q2', 'Q3', 'Q4']  # List of quadrants
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Reference transit depth (ppm)
td_ref = 6.72 * 0.46**2  # Transit duration in hours
depth_sig_scaling = 3
gamma_factor_significance = 1  # Other value could be 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

# Loop through each quadrant, load corresponding .npy files, and compute metrics
for idx, quadrant_name in enumerate(quadrant_names):
    # Load .npy files for the current quadrant
    quadrant_nominal = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_nominal.npy")
    quadrant_secondary = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_secondary.npy")
    quadrant_extended = np.load(f"{dataDIR_quadrants}quadrant_{quadrant_name}_extended.npy")

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

    # Compute metrics for each target in the quadrant
    for j in range(n_targets_quadrant):
        dback = np.ones(10) * dback_ref  # Reference transit depth
        td = np.ones(10) * td_ref  # Reference transit duration

        # Flux significance for nominal mask
        eta_nom_bt_24_cameras[j, :] = (
            gamma_factor_significance * dback * quadrant_nominal[j, 17:27] * np.sqrt(td * ntr) /
            (quadrant_nominal[j, 7] * (1 - quadrant_nominal[j, 11]))
        )

        # Flux significance for extended mask
        eta_ext_bt_24_cameras[j, :] = (
            gamma_factor_significance * dback * quadrant_extended[j, 14:24] * np.sqrt(td * ntr) /
            (quadrant_extended[j, 4] * (1 - quadrant_extended[j, 13]))
        )

        # Observed transit depth
        delta_obs[j, :] = dback * quadrant_nominal[j, 17:27]
        delta_obs_ext[j, :] = dback * quadrant_extended[j, 14:24]

        # Compute significant transit depth variables
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

        # Significant transit depth calculations
        sig_depth_secondary_mask_24_cameras = nsr1h_sec * (1 - SPR_tot_sec) / np.sqrt(td_ref * ntr)
        sig_depth_extended_mask_24_cameras = quadrant_extended[j, 4] * (1 - SPR_tot_ext) / np.sqrt(td_ref * ntr)
        sig_depth_nominal_mask_24_cameras = nsr1h_nom * (1 - SPR_tot_nom) / np.sqrt(td_ref * ntr)

        # Quadratic sum of the noises (Equation (38) from the paper)
        sig_depth_24_cameras[j, :] = np.sqrt(
            sig_depth_nominal_mask_24_cameras ** 2 + sig_depth_extended_mask_24_cameras ** 2
        )
        sig_depth_sec_nom_quad[j, :] = np.sqrt(
            sig_depth_secondary_mask_24_cameras ** 2 + sig_depth_nominal_mask_24_cameras ** 2
        )

        # COB-related variables
        eta_cob_nom_10first_24_cameras = quadrant_nominal[j, 46:56]
        eta_cob_ext_10first_24_cameras = quadrant_extended[j, 45:55]

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
        eta_c = quadrant_secondary[j, 6]
        delta_obs_c = quadrant_secondary[j, 7]
        delta_obs_t = delta_obs[j, :]
        secondary_mask_conditions_24_cameras[j, :] = (
            (eta_c > flux_thresh_sec_mask) &
            (delta_obs_c > delta_obs_t + depth_sig_scaling * sig_depth_sec_nom_quad[j, :]) &
            fp_single_contaminant_24_cameras[j, :]
        )

    # Compute efficiencies
    nfp_total = nfp.sum()
    if nfp_total > 0:
        eff_ext_flux = (nfp & nfp_ext_mask).sum() / nfp_total * 100.0
        eff_nom_cob = (nfp & nfp_nom_cob).sum() / nfp_total * 100.0
        eff_ext_cob = (nfp & nfp_ext_cob).sum() / nfp_total * 100.0
    else:
        eff_ext_flux = 0
        eff_nom_cob = 0
        eff_ext_cob = 0

    # Secondary mask efficiency
    fp_single_total = fp_single_contaminant_24_cameras.sum()
    if fp_single_total > 0:
        eff_secondary = secondary_mask_conditions_24_cameras.sum() / fp_single_total * 100.0
    else:
        eff_secondary = 0

    # Store efficiencies for each quadrant
    eff_extended_quadrant.append(eff_ext_flux)
    eff_secondary_quadrant.append(eff_secondary)
    eff_nom_cob_quadrant.append(eff_nom_cob)
    eff_ext_cob_quadrant.append(eff_ext_cob)

    # Print efficiencies for the current quadrant
    print(f"Quadrant {quadrant_name} Extended Flux Efficiency: {eff_ext_flux:.2f}%")
    print(f"Quadrant {quadrant_name} Secondary Mask Efficiency: {eff_secondary:.2f}%")
    print(f"Quadrant {quadrant_name} Nominal COB Efficiency: {eff_nom_cob:.2f}%")
    print(f"Quadrant {quadrant_name} Extended COB Efficiency: {eff_ext_cob:.2f}%\n")

# Plotting efficiencies across quadrants
plt.figure(figsize=(10, 6))
x = np.arange(len(quadrant_names))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots(figsize=(10, 6))

rects1 = ax.bar(x - 1.5*width, eff_nom_cob_quadrant, width, label='Nominal COB Efficiency', color='blue')
rects2 = ax.bar(x - 0.5*width, eff_ext_cob_quadrant, width, label='Extended COB Efficiency', color='green')
rects3 = ax.bar(x + 0.5*width, eff_secondary_quadrant, width, label='Secondary Flux Efficiency', color='red')
rects4 = ax.bar(x + 1.5*width, eff_extended_quadrant, width, label='Extended Flux Efficiency', color='purple')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Quadrant')
ax.set_ylabel('Efficiency (%)')
ax.set_title('Efficiencies across Quadrants')
ax.set_xticks(x)
ax.set_xticklabels(quadrant_names)
ax.legend()
ax.grid(True)

plt.tight_layout()
plt.show()