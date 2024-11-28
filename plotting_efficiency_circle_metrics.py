import numpy as np #type:ignore
import matplotlib.pyplot as plt # type:ignore

# Some files parameters
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


number_of_circles = 6  # Adjust as needed
ntr = 3  # Number of transits in one hour
dback_ref = 132000  # Example reference value for transit depth (ppm)
td_ref = 6.72 * 0.46**2  # Example transit duration in hours
seed = 123434434
depth_sig_scaling = 3
gamma_factor_significance = 1 ## other 0.46
flux_thresh_nom_mask, cob_thresh = 6, 3
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3

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

     # Initialize variables to count metrics satisfying conditions for efficiency
    count_nominal = 0
    count_extended = 0
    count_secondary = 0

    # Initialize counts for COB efficiencies
    count_nom_cob = 0
    count_ext_cob = 0

     # Compute metrics for each target in the ring
    for j in range(n_targets_ring):
        dback = np.ones(10) * dback_ref  # Reference transit depth
        td = np.ones(10) * td_ref  # Reference transit duration

        # Flux significance for nominal mask (assuming columns 17-27 and 7, 11 are correct for nominal data)
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

    # Compute efficiencies as per your method
    nfp_total = nfp.sum()
    if nfp_total > 0:
        eff_ext_flux = (nfp & nfp_ext_mask).sum() / nfp_total * 100.
        eff_nom_cob = (nfp & nfp_nom_cob).sum() / nfp_total * 100.
        eff_ext_cob = (nfp & nfp_ext_cob).sum() / nfp_total * 100.
    else:
        eff_ext_flux = 0
        eff_nom_cob = 0
        eff_ext_cob = 0

    # Secondary mask efficiency
    fp_single_total = fp_single_contaminant_24_cameras.sum()
    if fp_single_total > 0:
        eff_secondary = secondary_mask_conditions_24_cameras.sum() / fp_single_total * 100.
    else:
        eff_secondary = 0

    # Store efficiencies for each ring
    eff_extended_ring.append(eff_ext_flux)
    eff_secondary_ring.append(eff_secondary)
    eff_nom_cob_ring.append(eff_nom_cob)
    eff_ext_cob_ring.append(eff_ext_cob)

    # Print efficiencies for the current ring
    print(f"Ring {i} Extended Flux Efficiency: {eff_ext_flux:.2f}%")
    print(f"Ring {i} Secondary Mask Efficiency: {eff_secondary:.2f}%")
    print(f"Ring {i} Nominal COB Efficiency: {eff_nom_cob:.2f}%")
    print(f"Ring {i} Extended COB Efficiency: {eff_ext_cob:.2f}%")

# Plotting efficiencies across rings with 'o-' format
plt.figure(figsize=(10, 6))
plt.plot(range(number_of_circles), eff_nom_cob_ring, 'o-', color='blue', label="Nominal COB Efficiency")
plt.plot(range(number_of_circles), eff_ext_cob_ring, 's-', color='green', label="Extended COB Efficiency")
plt.plot(range(number_of_circles), eff_secondary_ring, '^-', color='red', label="Secondary Flux Efficiency")
plt.plot(range(number_of_circles), eff_extended_ring, 'd-', color='purple', label="Extended Flux Efficiency")
plt.xlabel("Ring Index")
plt.ylabel("Efficiency (%)")
plt.title("Efficiencies across Rings")
plt.legend()
plt.grid(True)
plt.show()                          