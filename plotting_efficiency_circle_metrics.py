import numpy as np #type:ignore
import matplotlib.pyplot as plt #type:ignore



# Some files parameters
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/rings/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'
# Initialize variables for storing ring-based results
eta_ext_bt_24_cameras_ring = []
eta_ext_bt_6_cameras_ring = []
eta_nom_bt_24_cameras_ring = []
eta_nom_bt_6_cameras_ring = []
delta_obs_ring = []
delta_obs_ext_ring = []
delta_obs_ext_6_cameras_ring = []


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
    delta_obs_ext_6_cameras = np.zeros((n_targets_ring, 10))

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
        delta_obs_ext_6_cameras[j, :] = dback * ring_extended[j, 14:24]

    # Append results for this ring to main lists
    eta_ext_bt_24_cameras_ring.append(eta_ext_bt_24_cameras)
    eta_ext_bt_6_cameras_ring.append(eta_ext_bt_6_cameras)
    eta_nom_bt_24_cameras_ring.append(eta_nom_bt_24_cameras)
    eta_nom_bt_6_cameras_ring.append(eta_nom_bt_6_cameras)
    delta_obs_ring.append(delta_obs)
    delta_obs_ext_ring.append(delta_obs_ext)
    delta_obs_ext_6_cameras_ring.append(delta_obs_ext_6_cameras)



                                       

