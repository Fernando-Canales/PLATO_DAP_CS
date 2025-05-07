import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from mpl_toolkits.axes_grid1.inset_locator import inset_axes #type: ignore
from matplotlib.ticker import FuncFormatter  # type: ignore

from imagette import ran_unique_int

dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Fixed_transit_depths_and_durations/magnitude_bins/fixed_dback_132000ppm_and_td_1_422_hr/1000_targets_per_magnitude_bin/different_PSFs/Noblesse_PSF/' 
#dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/1000_targets_per_magnitude_bin/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/'
#dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Long_Observational_Phase_Nord/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/'
# Parameters for the plots
Pmin = 10
Pmax = 13
binsize = 0.5
nP = int((Pmax - Pmin) / binsize + 1)
fsize = 14
flux_thresh_nom_mask, cob_thresh = 7.1, 3 
flux_thresh_ext_mask, flux_thresh_sec_mask = 3, 3
n_tar = 1000
depth_sig_scaling = 3
gamma_factor_significance = 1 ## other 0.46
td_ref = 6.72*0.46**2
dback_ref = 132000

#data_catalogue = np.load(cataDIR + 'LOPN1_DR3_20241011_gr0.npy') # star catalogue from GAIA
data_catalogue = np.load(cataDIR + 'SFP_DR3_20230101.npy')
magnitude_all_stars = data_catalogue[:, 2]
data = np.load(dataDIR + 'targets_P5.npy')

data_sec = np.load(dataDIR + 'targets_P5_secondary.npy')

data_ext = np.load(dataDIR + 'targets_P5_extended.npy')

data_bray = np.load(dataDIR + 'targets_P5_bray.npy')


eta_cob_ext_2_pix = data_ext[:, 244]
sigma_1_24_ext_2_pix = data_ext[:, 245]
delta_cob_ext_2_pix = data_ext[:, 246]

mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])
ntr = 3        # number of transits in one hour
n = data.shape[0]
td_from_dap = data[:, 126:136]
dback_from_dap = data[:, 136:146]
seed = 123434434

# We obtain the magnitude of all the targets and the magnitude of the most problematic contaminants'
mag = data[:, 1]
mag_bad = data[:, 3]

# We obtain the number of contaminant stars that could create a false positive for each target (i.e. N_bad)
n_bad = data[:, 8]

# We create a useful mask for getting the magnitude range for the P5 sample only (P = 10.66 - 12.66)
mask_p5 = (mag >= 10) & (mag <= 13)

# Now we apply the mask for getting the magnitude range
mag_p5 = mag[mask_p5]

# Now we apply the mask again to estimate the number of N_bad in the P5 sample magnitude range given the nominal mask
n_bad_p5 = n_bad[mask_p5]

# Now we plot a percentage histogram like the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

# Now we get all the etas and delta_obs
eta_t = data[:, 12]
eta_t_24_cameras = data[:, 37]
eta_t_6_cameras_earth_like = data[:, 38]
eta_t_24_cameras_superearth = data[:, 39]
eta_t_6_cameras_superearth = data[:, 40]
eta_t_6_cameras = data[:, 45]
delta_obs_t = data[:, 13]
eta_c = data_sec[:, 6]
eta_c_6_cameras = data_sec[:, 11]
delta_obs_c = data_sec[:, 7]
eta_ext_for_the_most_significant_contaminant = data_ext[:, 7]
delta_obs_ext_single_contaminant = data_ext[:, 8]
nsr1h = data[:, 7]
nsr1h_sec = data_sec[:, 4]
nsr1h_sec_6_cameras = data_sec[:, 18]
spr_crit = data[:, 9]
nsr1h_ext = data_ext[:, 4]
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]

# We also obtain the shape and size of every mask
key_nom = data[:, 5]
key_sec = data_sec[:, 2]
key_ext = data_ext[:, 2]
size_nom = data[:, 6]
size_sec = data_sec[:, 3]
size_e = data_ext[:, 3]

# Now all the metrics related to the COB taking into account only one contaminant
delta_cob = data[:, 14]
eta_cob = data[:, 15]
sigma_cob = data[:, 16]
eta_cob_6_cameras = data[:, 41]
delta_cob_6_cameras = data[:, 42]
sigma_cob_6_cameras = data[:, 43]
delta_cob_sec = data_sec[:, 8]
eta_cob_sec_24_cameras = data_sec[:, 9]
sigma_cob_sec_24_cameras = data_sec[:, 10]
eta_cob_sec_6_cameras = data_sec[:, 12]
delta_cob_sec_6_cameras = data_sec[:, 13]
sigma_cob_sec_6_cameras = data_sec[:, 14]

# We get now the 10 first values for each param. of the nominal cob shift
SPRK10_first = data[:, 17:27]
eta_cob_nom_10first_24_cameras = data[:, 46:56]
sigma_cob_10first_24_cameras = data[:, 56:66]
abs_cob_shift_10first = data[:, 66:76]
eta_cob_nom_10first_6_cameras = data[:, 76:86]
sigma_cob_10first_6_cameras = data[:, 86:96]
abs_cob_shift_10first = data[:, 96:106]
delta_cob_ext = data_ext[:, 10]
eta_cob_ext = data_ext[:, 11]
sigma_cob_ext = data_ext[:, 12]
gamma_cob_nom_10first_24_cameras = data[:, 106:116]


# We get now the 10 first values for each param. of the extended cob shift
SPRK10_first_ext = data_ext[:,14:24]
eta_cob_ext_10first_24_cameras = data_ext[:, 45:55]
sigma_cob_ext_10first_24_cameras = data_ext[:, 55:65]
delta_cob_ext_10first_24_cameras = data_ext[:, 65:75]
eta_cob_ext_10first_6_cameras = data_ext[:, 75:85]
sigma_cob_ext_10first_6_cameras = data_ext[:, 85:95]
delta_cob_ext_10first_6_cameras = data_ext[:, 95:105]
gamma_cob_ext_10first_24_cameras = data_ext[:, 105:115]
mag_target_10first_contaminants = data[:, 169:179]
dist_target_to_10first_contaminants = data[:, 179:189]

# We get now the new varaibles for the significant transit depth
sig_depth_secondary_mask_24_cameras = nsr1h_sec*(1 - data_sec[:, 5])/np.sqrt(td_ref*ntr) # NSR*(1-SPRtot)
sig_depth_secondary_mask_6_cameras = nsr1h_sec_6_cameras*(1 - data_sec[:, 5])/np.sqrt(td_ref*ntr) # NSR*(1-SPRtot)
sig_depth_extended_mask_24_cameras = data_ext[:, 4]*(1 - data_ext[:, 13])/np.sqrt(td_ref*ntr) # NSR_1h^ext*(1-SPR_tot^ext) / sqrt(td ntr) [24 cameras]
sig_depth_extended_mask_6_cameras = data_ext[:, 44]*(1 - data_ext[:, 13])/np.sqrt(td_ref*ntr) # NSR_1h^ext*(1-SPR_tot^ext) / sqrt(td ntr) [6 cameras]
sig_depth_nominal_mask_24_cameras = data[:, 7]*(1 - data[:, 11])/np.sqrt(td_ref*ntr) # NSR_1h^nom*(1-SPR_tot^nom) / sqrt(td ntr) [24 cameras]
sig_depth_nominal_mask_6_cameras = data[:, 148]*(1 - data[:, 11])/np.sqrt(td_ref*ntr) # NSR_1h^nom*(1-SPR_tot^nom) / sqrt(td ntr) [6 cameras]

# We obtainn the value of the quadratic sum of the noises (Eq. (38) of the paper)
sig_depth_24_cameras = np.sqrt(sig_depth_nominal_mask_24_cameras**2 + sig_depth_extended_mask_24_cameras**2)
sig_depth_6_cameras = np.sqrt(sig_depth_nominal_mask_6_cameras**2 + sig_depth_extended_mask_6_cameras**2)

# First, the expressions for the efficiency false positives detections for the different masks
fp_single_contaminant_24_cameras = (eta_t > flux_thresh_nom_mask)  # false positive
sdr_flux = fp_single_contaminant_24_cameras & (eta_c > flux_thresh_nom_mask) & (delta_obs_c > delta_obs_t)  # secondary mask false positive detection rate
fp_single_contaminant_6_cameras = (eta_t_6_cameras > flux_thresh_nom_mask)
secondary_mask_conditions_24_cameras = (eta_c > flux_thresh_sec_mask) & (delta_obs_c > delta_obs_t  + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_24_cameras**2 + sig_depth_nominal_mask_24_cameras**2)) & fp_single_contaminant_24_cameras  # secondary mask efficiency condition for 24 cameras
secondary_mask_conditions_6_cameras = (eta_c_6_cameras > flux_thresh_sec_mask) & (delta_obs_c > delta_obs_t + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_6_cameras**2 + sig_depth_nominal_mask_6_cameras**2)) & fp_single_contaminant_6_cameras # secondary mask efficiency condition for 6 cameras
eficiency_extended_mask_highest_spr_contaminant = fp_single_contaminant_24_cameras & (eta_ext_for_the_most_significant_contaminant > flux_thresh_ext_mask) & (delta_obs_ext_single_contaminant > delta_obs_t + depth_sig_scaling*sig_depth_24_cameras)  # extended mask false positive detection rate
cd = fp_single_contaminant_24_cameras & (eta_cob > cob_thresh)  # nominal mask false positive detection rate via cob shift
cd_6_cameras = fp_single_contaminant_6_cameras & (eta_cob_6_cameras > cob_thresh)
secondary_mask_conditions_cob_24_cameras = (eta_cob_sec_24_cameras > cob_thresh) & fp_single_contaminant_24_cameras  # secondary mask efficiency condition for 24 cameras and cob shift
secondary_mask_conditions_cob_6_cameras = (eta_cob_sec_6_cameras > cob_thresh) & fp_single_contaminant_6_cameras # secondary mask efficiency condition for 6 cameras and cob shift
ecd = fp_single_contaminant_24_cameras & (eta_cob_ext > cob_thresh)  # extended mask false positive detection rate via cob shift

# We reshape some of them for 24 cameras
sig_depth_extended_mask_24_cameras_10first = np.reshape(sig_depth_extended_mask_24_cameras, (n, 1))
sig_depth_extended_mask_24_cameras_10first = np.repeat(sig_depth_extended_mask_24_cameras_10first, 10, axis=1)
sig_depth_nominal_mask_24_cameras_10first = np.reshape(sig_depth_nominal_mask_24_cameras, (n, 1))
sig_depth_nominal_mask_24_cameras_10first = np.repeat(sig_depth_nominal_mask_24_cameras_10first, 10, axis=1)


# We reshape some of them for 6 cameras
sig_depth_extended_mask_6_cameras_10first = np.reshape(sig_depth_extended_mask_6_cameras, (n, 1))
sig_depth_extended_mask_6_cameras_10first = np.repeat(sig_depth_extended_mask_6_cameras_10first, 10, axis=1)
sig_depth_nominal_mask_6_cameras_10first = np.reshape(sig_depth_nominal_mask_6_cameras, (n, 1))
sig_depth_nominal_mask_6_cameras_10first = np.repeat(sig_depth_nominal_mask_6_cameras_10first, 10, axis=1)

sig_depth_24_cameras_10first = np.sqrt(sig_depth_nominal_mask_24_cameras_10first**2 + sig_depth_extended_mask_24_cameras_10first**2)
sig_depth_6_cameras_10first = np.sqrt(sig_depth_nominal_mask_6_cameras_10first**2 + sig_depth_extended_mask_6_cameras_10first**2)


nbad_sp = np.zeros(n) # small planet 2<R_Ee -> 4*84ppm
eta_ext_bt_24_cameras= np.zeros((n, 10))
eta_ext_bt_6_cameras = np.zeros((n,10))
eta_nom_bt_24_cameras= np.zeros((n, 10))
eta_nom_bt_6_cameras = np.zeros((n,10))
delta_obs = np.zeros((n,10))
delta_obs_ext = np.zeros((n,10))
delta_obs_ext_6_cameras = np.zeros((n, 10))
delta_obs_ext_2_pix_24_cameras = np.zeros((n, 10))
delta_obs_ext_2_pix_6_cameras = np.zeros((n, 10))
eta_ext_2_pix_bt_24_cameras = np.zeros((n, 10))
eta_ext_2_pix_bt_6_cameras = np.zeros((n, 10))

mag_2d = np.repeat(mag[:,np.newaxis], 10, axis=1)

for i in range(n):
    dback= dback_from_dap[i, ::] #np.ones(10)*dback_ref
    td = td_from_dap[i, ::] #np.ones(10)*td_ref
    eta_nom_bt_24_cameras[i, :] = gamma_factor_significance*dback*data[i, 17:27]*np.sqrt(td*ntr)/(data[i, 7]*(1 - data[i, 11])) # Eq.(12) from the paper
    eta_nom_bt_6_cameras[i, :] = gamma_factor_significance*dback*data[i, 17:27]*np.sqrt(td*ntr)/(data[i, 148]*(1 - data[i, 11]))
    eta_ext_bt_24_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 4] * (1 - data_ext[i, 13])) # Eq.(18) from the paper
    eta_ext_bt_6_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 44] * (1 - data_ext[i, 13]))
    eta_ext_2_pix_bt_24_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 126:136]*np.sqrt(td*ntr)/(data_ext[i, 238] * (1 - data_ext[i, 248]))
    eta_ext_2_pix_bt_6_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 126:136]*np.sqrt(td*ntr)/(data_ext[i, 156] * (1 - data_ext[i, 248]))
    delta_obs[i,:] = dback*SPRK10_first[i,:] # observed transit depth

    delta_obs[i,:] = dback*SPRK10_first[i,:] # observed transit depth
    #delta_int = delta_obs[i,:]/(1. -data_nommask[i,9] ) # inferred intrinsic transit depth
    delta_obs_ext[i,:] = dback*data_ext[i,14:24] # observed transit depth
    delta_obs_ext_6_cameras[i, :] = dback*data_ext[i,14:24] # observed transit depth with 6 cameras
    delta_obs_ext_2_pix_24_cameras[i, :] = dback*data_ext[i, 126:136]
    delta_obs_ext_2_pix_6_cameras[i, :] = dback*data_ext[i, 126:136]

    delta_int = delta_obs_t[i]/ (1 - data[i, 11])
    nbad_sp[i] = np.sum( (eta_nom_bt_24_cameras[i,:]>7.1) & (delta_int<4*84. ))

# Save the eta_nom_bt_24_camerasarray
np.save(dataDIR+'eta_bt_24_cameras.npy', eta_nom_bt_24_cameras)
np.save(dataDIR+'eta_nom_bt_6_cameras', eta_nom_bt_6_cameras)
np.save(dataDIR+'eta_ext_bt_24_cameras.npy', eta_ext_bt_24_cameras)
np.save(dataDIR+'eta_ext_bt_6_cameras.npy', eta_ext_bt_6_cameras)

n_bad_bray_p5 = n_bad_bray[mask_p5]

# Custom function to format y-axis ticks
def percentage_formatter(x, _):
    return f"{int(x * 100)}%"


"""
Now we plot the cumulative or total number of nominal mask shapes to address the total number of target stars 
"""
cumulative_total = []
P_values = []
for i in range(nP):
    Pi = Pmin + i * binsize
    # Apply the filter to select target stars for the current magnitude bin
    m = (mag <= Pi + binsize/2.) & fp_single_contaminant_24_cameras
    
    # Get unique mask shapes for each mask type
    unique_nominal = np.unique(np.array(key_nom[m], dtype=np.int64))
    unique_sec = np.unique(np.array(key_sec[m], dtype=np.int64))
    unique_ext = np.unique(np.array(key_ext[m], dtype=np.int64))
    
    # Combine the mask arrays
    combined_masks = np.concatenate((unique_nominal, unique_sec, unique_ext))
    
    # Get the unique mask shapes across all types
    total_masks = len(np.unique(combined_masks))
    
    # Store the cumulative total
    cumulative_total.append(total_masks)
    P_values.append(Pi)
    

"""
Now we plot the comparison between extended mask and the correct version of it as a function of the target P magnitude
"""
# To store whether a label has been added for each type
labels_added = {
    'sec_24': False,
    'sec_6': False,
    'ext_24': False,
    'ext_6': False,
    'ext_2_pix_24': False,
    'ext_2_pix_6': False
}

for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    
    # Data calculations
    fp_ext_overall_24_cameras = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)[m, :].sum()
    fp_ext_overall_6_cameras = (eta_nom_bt_6_cameras > flux_thresh_nom_mask)[m, :].sum()
    
    eff_ext_overall_24_cameras = ((eta_ext_bt_24_cameras > flux_thresh_ext_mask) & (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask))[m, :].sum() / fp_ext_overall_24_cameras * 100
    eff_ext_overall_6_cameras = ((eta_ext_bt_6_cameras > flux_thresh_ext_mask) & (delta_obs_ext_6_cameras > delta_obs + depth_sig_scaling * sig_depth_6_cameras_10first) & (eta_nom_bt_6_cameras > flux_thresh_nom_mask))[m, :].sum() / fp_ext_overall_6_cameras * 100
    eff_sec_6_cameras = (secondary_mask_conditions_6_cameras[m].sum() / fp_single_contaminant_6_cameras[m].sum()) * 100
    eff_sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100

    # Errors
    error = np.sqrt(eff_ext_overall_24_cameras * (100 - eff_ext_overall_24_cameras) / m.sum())
    error_6_cameras = np.sqrt(eff_ext_overall_6_cameras * (100 - eff_ext_overall_6_cameras) / m.sum())
    error_sec_6_cameras = np.sqrt(eff_sec_6_cameras * (100 - eff_sec_6_cameras) / m.sum())
    error_sec = np.sqrt(eff_sec * (100 - eff_sec) / m.sum())

"""
Now we plot the efficiency of the COB shift (all contaminants)
    
"""
# Initialize lists to store data points for the inset
Pi_values = []
eff_ext_cob_overall_values = []
eff_ext_cob_overall_6_cameras_values = []
eff_cob_values = []
eff_cob_6_cameras_values = []

for i in range(nP):
    Pi = Pmin + i * binsize
    Pi_values.append(Pi)  # Store Pi for later

    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    s_24_cameras = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)[m,:].sum()
    s_6_cameras = (eta_nom_bt_6_cameras > flux_thresh_nom_mask)[m,:].sum()
    eff_ext_cob_overall = ((eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask))[m,:].sum() / s_24_cameras * 100.
    eff_ext_cob_overall_6_cameras = ((eta_cob_ext_10first_6_cameras > cob_thresh) &  (eta_nom_bt_6_cameras > flux_thresh_nom_mask))[m,:].sum() / s_6_cameras * 100.
    eff_cob = ((eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask))[m,:].sum() / s_24_cameras * 100.
    eff_cob_6_cameras = ((eta_cob_nom_10first_6_cameras > cob_thresh) & (eta_nom_bt_6_cameras > flux_thresh_nom_mask))[m,:].sum() / s_6_cameras * 100.
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_sec_6_cameras = (secondary_mask_conditions_cob_6_cameras[m].sum() / fp_single_contaminant_6_cameras[m].sum()) * 100

    # Store values for the inset plot
    eff_ext_cob_overall_values.append(eff_ext_cob_overall)
    eff_ext_cob_overall_6_cameras_values.append(eff_ext_cob_overall_6_cameras)
    eff_cob_values.append(eff_cob)
    eff_cob_6_cameras_values.append(eff_cob_6_cameras)
        
    #Computing the errors:
    error_ext_cob = np.sqrt(eff_ext_cob_overall * (100 - eff_ext_cob_overall) / m.sum())
    error_ext_cob_6_cameras = np.sqrt(eff_ext_cob_overall_6_cameras * (100 - eff_ext_cob_overall_6_cameras) / m.sum())
    error_cob = np.sqrt(eff_cob * (100 - eff_cob) / m.sum())
    error_cob_6_cameras = np.sqrt(eff_cob_6_cameras * (100 - eff_cob_6_cameras) / m.sum())
    error_cob_sec = np.sqrt(eff_cob_sec * (100 - eff_cob_sec) / m.sum())
    error_cob_sec_6_cameras = np.sqrt(eff_cob_sec_6_cameras * (100 - eff_cob_sec_6_cameras) / m.sum())
    

nfp = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_mask = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) & (delta_obs_ext > delta_obs + depth_sig_scaling*sig_depth_24_cameras_10first) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_flux_without_significant_transit_depth_condition = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) &  (delta_obs_ext > delta_obs) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_mask_single_contaminant = (eficiency_extended_mask_highest_spr_contaminant)
nfp_nom_cob = (eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_cob = (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)

nfp_sec_mask = (secondary_mask_conditions_24_cameras)
nfp_sec_cob = (secondary_mask_conditions_cob_24_cameras)
nfp_highest_contaminant = (eta_nom_bt_24_cameras[:,0]>flux_thresh_nom_mask)
nfp_ext_mask_highest_contaminant = (eta_ext_bt_24_cameras[:,0]>flux_thresh_ext_mask) & (delta_obs_ext[:,0]>delta_obs[:,0]+depth_sig_scaling*sig_depth_24_cameras_10first[:,0]) & nfp_highest_contaminant

magnitudes = mag  # Using the variable with all the magnitude values of the targets

#################### CONSIDERING THE BIAS FROM THE STAR COUNT OF EACH MAGNITUDE BIN #########################

star_counts = np.zeros(nP, dtype=int)
total_stars = 0
# Loop over each magnitude bin
for i in range(nP):
    Pi = Pmin + i * binsize  # Define bin edge
    mask_bin = (magnitude_all_stars >= Pi - binsize / 2.) & (magnitude_all_stars <= Pi + binsize / 2.)

    # Count the number of stars in the current bin
    star_counts[i] = mask_bin.sum()

    # Count the total number of stars
    total_stars = total_stars + star_counts[i]

# Compute the weights for the current bin
weights = star_counts / float(total_stars)

# Initialize variables for weighted sums and counts
weighted_eff_ext_flux = 0
weighted_eff_sec_flux = 0
weighted_eff_ext_flux_single_contaminant = 0
weighted_eff_ext_flux_without_significant_transit_depth_condition = 0
weighted_eff_nom_cob = 0
weighted_eff_ext_cob = 0
weighted_eff_sec_cob = 0
weighted_fraction_fp_ext_cob_no_ext_flux = 0
weighted_fraction_fp_nom_cob_no_ext_flux = 0
weighted_fraction_fp_ext_flux_no_nom_cob = 0
weighted_fraction_fp_ext_cob_no_nom_cob = 0
weighted_fraction_fp_ext_flux_no_ext_cob = 0
weighted_fraction_fp_sec_flux_no_ext_flux = 0
weighted_fraction_fp_ext_flux_no_sec_flux = 0
weighted_fraction_fp_sec_flux_no_nom_cob = 0
weighted_fraction_fp_nom_cob_no_sec_flux = 0
weighted_fraction_fp_sec_flux_no_sec_cob = 0
weighted_fraction_fp_sec_cob_no_sec_flux = 0
weighted_fraction_fp_sec_flux_no_ext_cob = 0
weighted_fraction_fp_ext_cob_no_sec_flux = 0
weighted_fraction_fp_sec_flux_no_ext_cob = 0
weighted_fraction_fp_ext_cob_no_sec_flux = 0

weighted_fraction_fp_sec_cob_no_ext_flux = 0
weighted_fraction_fp_sec_cob_no_nom_cob = 0
weighted_fraction_fp_nom_cob_no_sec_cob = 0
weighted_fraction_fp_sec_cob_no_ext_cob = 0
weighted_fraction_fp_ext_cob_no_sec_cob = 0


# Initialize errors for fractions
weighted_variance_fraction_fp_ext_cob_no_ext_flux = 0
weighted_variance_fraction_fp_nom_cob_no_ext_flux = 0
weighted_variance_fraction_fp_ext_flux_no_nom_cob = 0
weighted_variance_fraction_fp_ext_cob_no_nom_cob = 0
weighted_variance_fraction_fp_ext_flux_no_ext_cob = 0
weighted_variance_fraction_fp_sec_flux_no_ext_flux = 0
weighted_variance_fraction_fp_ext_flux_no_sec_flux = 0
weighted_variance_fraction_fp_sec_flux_no_nom_cob = 0
weighted_variance_fraction_fp_nom_cob_no_sec_flux = 0
weighted_variance_fraction_fp_sec_flux_no_sec_cob = 0
weighted_variance_fraction_fp_sec_cob_no_sec_flux = 0
weighted_variance_fraction_fp_sec_flux_no_ext_cob = 0
weighted_variance_fraction_fp_ext_cob_no_sec_flux = 0
weighted_variance_fraction_fp_sec_cob_no_ext_flux = 0
weighted_variance_fraction_fp_sec_cob_no_ext_flux = 0
weighted_variance_fraction_fp_sec_cob_no_nom_cob = 0
weighted_variance_fraction_fp_nom_cob_no_sec_cob = 0
weighted_variance_fraction_fp_sec_cob_no_ext_cob = 0
weighted_variance_fraction_fp_ext_cob_no_sec_cob = 0
weighted_variance_fraction_fp_sec_flux_no_ext_cob = 0
weighted_variance_fraction_fp_ext_cob_no_sec_flux = 0

# Initialize variance for efficiencies
weighted_variance_eff_ext_flux = 0
weighted_variance_eff_sec_flux = 0
weighted_variance_eff_ext_flux_single_contaminant = 0
weighted_variance_eff_ext_flux_without_significant_transit_depth_condition = 0
weighted_variance_eff_nom_cob = 0
weighted_variance_eff_ext_cob = 0
weighted_variance_eff_sec_cob = 0

n_star = np.zeros(nP)
# Loop over each magnitude bin
for i in range(nP):
    Pi = Pmin + i * binsize  # Define bin edge
    m = (mag >= Pi - binsize / 2.) & (mag <= Pi + binsize / 2.)
    n_star[i] = m.sum()

    # Compute efficiencies and fractions for the current bin
    eff_ext_flux = nfp_ext_mask[m].sum() / nfp[m].sum() * 100.
    eff_sec_flux = (nfp_sec_mask[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100
    eff_ext_flux_single_contaminant = (nfp_ext_mask_single_contaminant[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100.
    eff_ext_flux_without_significant_transit_depth_condition = (nfp[m] & nfp_ext_flux_without_significant_transit_depth_condition[m]).sum() / nfp[m].sum() * 100.
    eff_nom_cob = nfp_nom_cob[m].sum() / nfp[m].sum() * 100.
    eff_ext_cob = nfp_ext_cob[m].sum() / nfp[m].sum() * 100. 
    eff_sec_cob = (secondary_mask_conditions_cob_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100.
    fraction_fp_ext_cob_no_ext_flux = ((nfp_ext_mask[m]==False) & nfp_ext_cob[m] & nfp[m]).sum()/nfp[m].sum()
    fraction_fp_nom_cob_no_ext_flux = ((nfp_ext_mask[m]==False) & nfp_nom_cob[m] & nfp[m]).sum()/nfp[m].sum()
    fraction_fp_ext_flux_no_nom_cob = ((nfp_nom_cob[m]==False) & nfp_ext_mask[m] & nfp[m]).sum()/nfp[m].sum()
    fraction_fp_ext_flux_no_ext_cob = ((nfp_ext_cob[m]==False) & nfp_ext_mask[m] & nfp[m]).sum()/nfp[m].sum()
    fraction_fp_ext_cob_no_nom_cob = ((nfp_nom_cob[m]==False) & nfp_ext_cob[m] & nfp[m]).sum()/nfp[m].sum()
    # Calculate the fraction of FPs detected by secondary flux but not by extended flux (first column only)
    fraction_fp_sec_flux_no_ext_flux = ((nfp_ext_mask[m][:, 0] == False) & nfp_sec_mask[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    # Calculate the reverse comparison: FPs detected by extended flux but not by secondary flux
    fraction_fp_ext_flux_no_sec_flux = ((nfp_sec_mask[m] == False) & nfp_ext_mask[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_flux_no_nom_cob = ((nfp_nom_cob[m][:, 0] == False) & nfp_sec_mask[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_nom_cob_no_sec_flux = ((nfp_sec_mask[m] == False) & nfp_nom_cob[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_flux_no_sec_cob = ((nfp_sec_cob[m] == False) & nfp_sec_mask[m] & fp_single_contaminant_24_cameras[m]).sum()  / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_cob_no_sec_flux = ((nfp_sec_mask[m] == False) & nfp_sec_cob[m] & fp_single_contaminant_24_cameras[m]).sum()  / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_flux_no_ext_cob = ((nfp_ext_cob[m][:, 0] == False) & nfp_sec_mask[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_ext_cob_no_sec_flux = ((nfp_sec_mask[m] == False) & nfp_ext_cob[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_cob_no_ext_flux = ((nfp_ext_mask[m][:, 0] == False) & secondary_mask_conditions_cob_24_cameras[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_cob_no_nom_cob = ((nfp_nom_cob[m][:, 0] == False) & secondary_mask_conditions_cob_24_cameras[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_nom_cob_no_sec_cob = ((secondary_mask_conditions_cob_24_cameras[m] == False) & nfp_nom_cob[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_cob_no_ext_cob = ((nfp_ext_cob[m][:, 0] == False) & secondary_mask_conditions_cob_24_cameras[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_ext_cob_no_sec_cob = ((secondary_mask_conditions_cob_24_cameras[m] == False) & nfp_ext_cob[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_sec_flux_no_ext_cob = ((nfp_ext_cob[m][:, 0] == False) & secondary_mask_conditions_24_cameras[m] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()
    fraction_fp_ext_cob_no_sec_flux = ((secondary_mask_conditions_24_cameras[m] == False) & nfp_ext_cob[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()


    # Accumulate weighted sums
    weighted_eff_ext_flux +=  weights[i] * eff_ext_flux
    weighted_eff_sec_flux +=  weights[i] * eff_sec_flux
    weighted_eff_ext_flux_single_contaminant += weights[i] * eff_ext_flux_single_contaminant
    weighted_eff_ext_flux_without_significant_transit_depth_condition += weights[i] * eff_ext_flux_without_significant_transit_depth_condition
    weighted_eff_nom_cob += weights[i] * eff_nom_cob
    weighted_eff_ext_cob += weights[i] * eff_ext_cob
    weighted_eff_sec_cob += weights[i] * eff_sec_cob
    weighted_fraction_fp_ext_cob_no_ext_flux += weights[i] * fraction_fp_ext_cob_no_ext_flux
    weighted_fraction_fp_nom_cob_no_ext_flux += weights[i] * fraction_fp_nom_cob_no_ext_flux
    weighted_fraction_fp_ext_flux_no_nom_cob += weights[i] * fraction_fp_ext_flux_no_nom_cob
    weighted_fraction_fp_ext_cob_no_nom_cob += weights[i] * fraction_fp_ext_cob_no_nom_cob
    weighted_fraction_fp_ext_flux_no_ext_cob += weights[i] * fraction_fp_ext_flux_no_ext_cob
    weighted_fraction_fp_sec_flux_no_ext_flux += weights[i] * fraction_fp_sec_flux_no_ext_flux
    weighted_fraction_fp_ext_flux_no_sec_flux += weights[i] * fraction_fp_ext_flux_no_sec_flux
    weighted_fraction_fp_sec_flux_no_nom_cob += weights[i] * fraction_fp_sec_flux_no_nom_cob
    weighted_fraction_fp_nom_cob_no_sec_flux += weights[i] * fraction_fp_nom_cob_no_sec_flux
    weighted_fraction_fp_sec_flux_no_sec_cob += weights[i] * fraction_fp_sec_flux_no_sec_cob
    weighted_fraction_fp_sec_cob_no_sec_flux += weights[i] * fraction_fp_sec_cob_no_sec_flux 
    weighted_fraction_fp_sec_cob_no_ext_flux += weights[i] * fraction_fp_sec_cob_no_ext_flux
    weighted_fraction_fp_sec_cob_no_nom_cob += weights[i] * fraction_fp_sec_cob_no_nom_cob
    weighted_fraction_fp_nom_cob_no_sec_cob += weights[i] * fraction_fp_nom_cob_no_sec_cob
    weighted_fraction_fp_sec_cob_no_ext_cob += weights[i] * fraction_fp_sec_cob_no_ext_cob
    weighted_fraction_fp_ext_cob_no_sec_cob += weights[i] * fraction_fp_ext_cob_no_sec_cob
    weighted_fraction_fp_sec_flux_no_ext_cob += weights[i] * fraction_fp_sec_flux_no_ext_cob
    weighted_fraction_fp_ext_cob_no_sec_flux += weights[i] * fraction_fp_ext_cob_no_sec_flux


    # Calculate variance for the efficiencies
    variance_eff_ext_flux = (eff_ext_flux * (100 - eff_ext_flux)) / n_star[i]
    variance_eff_sec_flux = (eff_sec_flux * (100 - eff_sec_flux)) / n_star[i]
    variance_eff_ext_flux_single_contaminant = (eff_ext_flux_single_contaminant * (100 - eff_ext_flux_single_contaminant)) / n_star[i]
    variance_eff_ext_flux_without_significant_transit_depth_condition = (eff_ext_flux_without_significant_transit_depth_condition * (100 - eff_ext_flux_without_significant_transit_depth_condition)) / n_star[i]
    variance_eff_nom_cob = (eff_nom_cob * (100 - eff_nom_cob)) / n_star[i]
    variance_eff_ext_cob = (eff_ext_cob * (100 - eff_ext_cob)) / n_star[i]
    variance_eff_sec_cob = (eff_sec_cob * (100 - eff_sec_cob)) / n_star[i]

    # Accumulate weighted variance
    weighted_variance_eff_ext_flux += weights[i] * variance_eff_ext_flux
    weighted_variance_eff_sec_flux += weights[i] * variance_eff_sec_flux
    weighted_variance_eff_ext_flux_single_contaminant += weights[i] * variance_eff_ext_flux_single_contaminant
    weighted_variance_eff_ext_flux_without_significant_transit_depth_condition += weights[i] * variance_eff_ext_flux_without_significant_transit_depth_condition 
    weighted_variance_eff_nom_cob += weights[i] * variance_eff_nom_cob
    weighted_variance_eff_ext_cob += weights[i] * variance_eff_ext_cob
    weighted_variance_eff_sec_cob += weights[i] * variance_eff_sec_cob

    # Calculate variance for the fractions
    variance_fraction_fp_ext_cob_no_ext_flux = (fraction_fp_ext_cob_no_ext_flux * (1 - fraction_fp_ext_cob_no_ext_flux)) / n_star[i]
    variance_fraction_fp_nom_cob_no_ext_flux = (fraction_fp_nom_cob_no_ext_flux * (1 - fraction_fp_nom_cob_no_ext_flux)) / n_star[i]
    variance_fraction_fp_ext_flux_no_nom_cob = (fraction_fp_ext_flux_no_nom_cob * (1 - fraction_fp_ext_flux_no_nom_cob)) / n_star[i]
    variance_fraction_fp_ext_cob_no_nom_cob = (fraction_fp_ext_cob_no_nom_cob *(1 - fraction_fp_ext_cob_no_nom_cob)) / n_star[i]
    variance_fraction_fp_ext_flux_no_ext_cob = (fraction_fp_ext_flux_no_ext_cob * (1 - fraction_fp_ext_flux_no_ext_cob)) / n_star[i]
    # Calculate variances
    variance_fraction_fp_sec_flux_no_ext_flux = (fraction_fp_sec_flux_no_ext_flux * (1 - fraction_fp_sec_flux_no_ext_flux)) / n_star[i]
    variance_fraction_fp_ext_flux_no_sec_flux = (fraction_fp_ext_flux_no_sec_flux * (1 - fraction_fp_ext_flux_no_sec_flux)) / n_star[i]
    variance_fraction_fp_sec_flux_no_nom_cob = (fraction_fp_sec_flux_no_nom_cob * (1 - fraction_fp_sec_flux_no_nom_cob)) / n_star[i]
    variance_fraction_fp_nom_cob_no_sec_flux = (fraction_fp_nom_cob_no_sec_flux * (1 - fraction_fp_nom_cob_no_sec_flux)) / n_star[i]
    variance_fraction_fp_sec_flux_no_sec_cob = (fraction_fp_sec_flux_no_sec_cob * (1 - fraction_fp_sec_flux_no_sec_cob)) / n_star[i]
    variance_fraction_fp_sec_cob_no_sec_flux = (fraction_fp_sec_cob_no_sec_flux * (1 - fraction_fp_sec_cob_no_sec_flux)) / n_star[i]
    variance_fraction_fp_sec_cob_no_ext_flux = (fraction_fp_sec_cob_no_ext_flux * (1 - fraction_fp_sec_cob_no_ext_flux)) / n_star[i]
    variance_fraction_fp_sec_cob_no_nom_cob = (fraction_fp_sec_cob_no_nom_cob * (1 - fraction_fp_sec_cob_no_nom_cob)) / n_star[i]
    variance_fraction_fp_nom_cob_no_sec_cob = (fraction_fp_nom_cob_no_sec_cob * (1 - fraction_fp_nom_cob_no_sec_cob)) / n_star[i]
    variance_fraction_fp_sec_cob_no_ext_cob = (fraction_fp_sec_cob_no_ext_cob * (1 - fraction_fp_sec_cob_no_ext_cob)) / n_star[i]
    variance_fraction_fp_ext_cob_no_sec_cob = (fraction_fp_ext_cob_no_sec_cob * (1 - fraction_fp_ext_cob_no_sec_cob)) / n_star[i]
    variance_fraction_fp_sec_flux_no_ext_cob = (fraction_fp_sec_flux_no_ext_cob * (1 - fraction_fp_sec_flux_no_ext_cob)) / n_star[i]
    variance_fraction_fp_ext_cob_no_sec_flux = (fraction_fp_ext_cob_no_sec_flux * (1 - fraction_fp_ext_cob_no_sec_flux)) / n_star[i]

     

    # Accumulate weighted variance for the fractions
    weighted_variance_fraction_fp_ext_cob_no_ext_flux += weights[i] * variance_fraction_fp_ext_cob_no_ext_flux
    weighted_variance_fraction_fp_nom_cob_no_ext_flux += weights[i] * variance_fraction_fp_nom_cob_no_ext_flux
    weighted_variance_fraction_fp_ext_flux_no_nom_cob += weights[i] * variance_fraction_fp_ext_flux_no_nom_cob
    weighted_variance_fraction_fp_ext_cob_no_nom_cob += weights[i] * variance_fraction_fp_ext_cob_no_nom_cob
    weighted_variance_fraction_fp_ext_flux_no_ext_cob += weights[i] * variance_fraction_fp_ext_flux_no_ext_cob
    weighted_variance_fraction_fp_sec_flux_no_ext_flux += weights[i] * variance_fraction_fp_sec_flux_no_ext_flux
    weighted_variance_fraction_fp_ext_flux_no_sec_flux += weights[i] * variance_fraction_fp_ext_flux_no_sec_flux
    weighted_variance_fraction_fp_sec_flux_no_nom_cob += weights[i] * variance_fraction_fp_sec_flux_no_nom_cob
    weighted_variance_fraction_fp_nom_cob_no_sec_flux += weights[i] * variance_fraction_fp_nom_cob_no_sec_flux
    weighted_variance_fraction_fp_sec_flux_no_sec_cob += weights[i] * variance_fraction_fp_sec_flux_no_sec_cob
    weighted_variance_fraction_fp_sec_cob_no_sec_flux += weights[i] * variance_fraction_fp_sec_cob_no_sec_flux
    weighted_variance_fraction_fp_sec_cob_no_ext_flux += weights[i] * variance_fraction_fp_sec_cob_no_ext_flux
    weighted_variance_fraction_fp_sec_cob_no_nom_cob += weights[i] * variance_fraction_fp_sec_cob_no_nom_cob
    weighted_variance_fraction_fp_nom_cob_no_sec_cob += weights[i] * variance_fraction_fp_nom_cob_no_sec_cob
    weighted_variance_fraction_fp_sec_cob_no_ext_cob += weights[i] * variance_fraction_fp_sec_cob_no_ext_cob
    weighted_variance_fraction_fp_ext_cob_no_sec_cob += weights[i] * variance_fraction_fp_ext_cob_no_sec_cob
    weighted_variance_fraction_fp_sec_flux_no_ext_cob += weights[i] * variance_fraction_fp_sec_flux_no_ext_cob
    weighted_variance_fraction_fp_ext_cob_no_sec_flux += weights[i] * variance_fraction_fp_ext_cob_no_sec_flux




# Compute the final errors as the square root of the weighted variances
weighted_error_eff_ext_flux = np.sqrt(weighted_variance_eff_ext_flux)
weighted_error_eff_sec_flux = np.sqrt(weighted_variance_eff_sec_flux)
weighted_error_eff_ext_flux_single_contaminant = np.sqrt(weighted_variance_eff_ext_flux_single_contaminant)
weighted_error_eff_ext_flux_without_significant_transit_depth_condition = np.sqrt(weighted_variance_eff_ext_flux_without_significant_transit_depth_condition)
weighted_error_eff_nom_cob = np.sqrt(weighted_variance_eff_nom_cob)
weighted_error_eff_ext_cob = np.sqrt(weighted_variance_eff_ext_cob)
weighted_error_eff_sec_cob = np.sqrt(weighted_variance_eff_sec_cob)

# Computing the errors (standard deviaiton -> square root of the variance previously obtained)
weighted_error_fraction_fp_ext_cob_no_ext_flux = np.sqrt(weighted_variance_fraction_fp_ext_cob_no_ext_flux)
weighted_error_fraction_fp_nom_cob_no_ext_flux = np.sqrt(weighted_variance_fraction_fp_nom_cob_no_ext_flux)
weighted_error_fraction_fp_ext_flux_no_nom_cob = np.sqrt(weighted_variance_fraction_fp_ext_flux_no_nom_cob)
weighted_error_fraction_fp_ext_cob_no_nom_cob = np.sqrt(weighted_variance_fraction_fp_ext_cob_no_nom_cob)
weighted_error_fraction_fp_ext_flux_no_ext_cob = np.sqrt(weighted_variance_fraction_fp_ext_flux_no_ext_cob)
weighted_error_fraction_fp_sec_flux_no_ext_flux = np.sqrt(weighted_variance_fraction_fp_sec_flux_no_ext_flux)
weighted_error_fraction_fp_ext_flux_no_sec_flux = np.sqrt(weighted_variance_fraction_fp_ext_flux_no_sec_flux)
weighted_error_fraction_fp_sec_flux_no_nom_cob = np.sqrt(weighted_variance_fraction_fp_sec_flux_no_nom_cob)
weighted_error_fraction_fp_nom_cob_no_sec_flux = np.sqrt(weighted_variance_fraction_fp_nom_cob_no_sec_flux)
weighted_error_fraction_fp_sec_flux_no_sec_cob = np.sqrt(weighted_variance_fraction_fp_sec_flux_no_sec_cob)
weighted_error_fraction_fp_sec_cob_no_sec_flux = np.sqrt(weighted_variance_fraction_fp_sec_cob_no_sec_flux)
weighted_error_fraction_fp_sec_cob_no_ext_flux = np.sqrt(weighted_variance_fraction_fp_sec_cob_no_ext_flux)
weighted_error_fraction_fp_sec_cob_no_nom_cob = np.sqrt(weighted_variance_fraction_fp_sec_cob_no_nom_cob)
weighted_error_fraction_fp_nom_cob_no_sec_cob = np.sqrt(weighted_variance_fraction_fp_nom_cob_no_sec_cob)
weighted_error_fraction_fp_sec_cob_no_ext_cob = np.sqrt(weighted_variance_fraction_fp_sec_cob_no_ext_cob)
weighted_error_fraction_fp_ext_cob_no_sec_cob = np.sqrt(weighted_variance_fraction_fp_ext_cob_no_sec_cob)
weighted_error_fraction_fp_sec_flux_no_ext_cob = np.sqrt(weighted_variance_fraction_fp_sec_flux_no_ext_cob)
weighted_error_fraction_fp_ext_cob_no_sec_flux = np.sqrt(weighted_variance_fraction_fp_ext_cob_no_sec_flux)

# Print weighted efficiency results with errors
print(f'Weighted extended flux efficiency: {weighted_eff_ext_flux:.2f}% ± {weighted_error_eff_ext_flux:.2f}%')
print(f'Weighted secondary flux efficiency: {weighted_eff_sec_flux:.2f}% ± {weighted_error_eff_sec_flux:.2f}%')
print(f'Weighted extended flux for a single contaminant efficiency: {weighted_eff_ext_flux_single_contaminant:.2f}% ± {weighted_error_eff_ext_flux_single_contaminant:.2f}%')
print(f'Weighted extended flux withouth the significant transit condition: {weighted_eff_ext_flux_without_significant_transit_depth_condition:.2f}% ± {weighted_error_eff_ext_flux_without_significant_transit_depth_condition:.2f}%')
print(f'Weighted nominal COB efficiency: {weighted_eff_nom_cob:.2f}% ± {weighted_error_eff_nom_cob:.2f}%')
print(f'Weighted extended COB efficiency: {weighted_eff_ext_cob:.2f}% ± {weighted_error_eff_ext_cob:.2f}%')
print(f'Weighted secondary mask COB efficiency: {weighted_eff_sec_cob:.2f}% ± {weighted_error_eff_sec_cob:.2f}%')


def calculate_effective_efficiency(data, data_sec, data_ext, 
                                  eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
                                  eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
                                  fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
                                  flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
                                  delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first,
                                  depth_sig_scaling):
    """
    Calculate the effective efficiency of the PLATO false positive detection strategy.
    
    Parameters:
    -----------
    data, data_sec, data_ext : numpy arrays
        Arrays containing target and contaminant data
    eta_nom_bt_24_cameras, eta_ext_bt_24_cameras : numpy arrays
        Arrays containing eta values for nominal and extended flux
    eta_c : numpy array
        Array containing eta values for secondary flux
    eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras : numpy arrays
        Arrays containing eta values for COB shifts
    fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras : numpy arrays
        Boolean masks for false positives
    flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh : float
        Threshold values for detection methods
    delta_obs, delta_obs_ext, delta_obs_t : numpy arrays
        Observed transit depths
    sig_depth_24_cameras_10first : numpy array
        Signal depth significance
    depth_sig_scaling : float
        Scaling factor for depth significance
        
    Returns:
    --------
    tuple
        Tuple containing (effective_efficiency, metric_counts, resource_usage, resource_limits, 
                         n_counts, method_efficiencies, total_fps, detected_fps, n_targets, within_limits)
    """
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("+                                                                                 +")
    print("+ Calculating effective efficiency of the PLATO false positive detection strategy +")
    print("+                                                                                 +")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
    # Extract necessary data
    mag = data[:, 1]  # Target magnitudes
    n_bad = data[:, 8]  # Number of potential FPs per target
    print(len(n_bad))
    n_targets = len(mag)
    
    
    # Arrays to store assigned metrics and detected FPs
    assigned_metrics = np.empty(n_targets, dtype='U4')  # EFX, SFX, NCOB, ECOB
    fps_detected_by_assigned = np.zeros(n_targets)
    total_fps_per_target = np.zeros(n_targets)
    
    # variables to count metric assignments
    count_efx = 0
    count_sfx = 0
    count_ncob = 0
    count_ecob = 0
    count_efx_zero = 0
    
    # variables to track FPs
    total_fps = 0 # Total count of FPs across all targets
    detected_fps = 0 # Total count of detected FPs across all targets 
    
    # Define detection capability for each method (for each target)
    #efx_detection = np.zeros((n_targets, 10), dtype=bool)
    #sfx_detection = np.zeros(n_targets, dtype=bool)
    #ncob_detection = np.zeros((n_targets, 10), dtype=bool)
    #ecob_detection = np.zeros((n_targets, 10), dtype=bool)
    
    # Calculate detection capability for each method
    efx_detection = (eta_ext_bt_24_cameras > flux_thresh_ext_mask) & \
                    (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    
    sfx_detection = secondary_mask_conditions_24_cameras
    
    ncob_detection = (eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    
    ecob_detection = (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    print(np.shape(eta_nom_bt_24_cameras))
    # For each target, follow the decision flow to assign the best metric
    for i in range(n_targets):
        # Count the number of potential FPs for this target
        N = n_bad[i]
        total_fps_per_target[i] = N
        total_fps += N
        
        # Implement the decision logic from the flow diagram
        if N == 0:
            # Assign Extended Flux (Low Priority)
            assigned_metrics[i] = "EFX"
            count_efx += 1
            # No FPs to detect
            fps_detected_by_assigned[i] = 0
            
        elif N == 1:
            # Check if SFX can detect the highest SPR contaminant
            if sfx_detection[i]:
                assigned_metrics[i] = "SFX"
                count_sfx += 1
                fps_detected_by_assigned[i] = 1
                #fps_detected_by_assigned[i] = min(fps_detected_by_assigned[i], total_fps_per_target[i])
                detected_fps += fps_detected_by_assigned[i]
            else:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                # Check if NCOB detects the FP
                #if ncob_detection[i, 0]:  # Check only highest SPR contaminant
                fps_detected_by_assigned[i] = 1
                    #fps_detected_by_assigned[i] = min(fps_detected_by_assigned[i], total_fps_per_target[i])
                detected_fps += fps_detected_by_assigned[i]
                #else:
                    #fps_detected_by_assigned[i] = 0
                    
        else:  # N ≥ 2
            # Count how many FPs each method can detect
            efx_count = np.sum(efx_detection[i])
            ncob_count = np.sum(ncob_detection[i])
            ecob_count = np.sum(ecob_detection[i])
            
            # Find the method that detects the most
            if efx_count >= ncob_count and efx_count >= ecob_count:
                assigned_metrics[i] = "EFX"
                count_efx += 1
                fps_detected_by_assigned[i] = efx_count
                if efx_count > N:
                    print('efx_count:', efx_count)
                    print('N is:', N)
                #fps_detected_by_assigned[i] = min(fps_detected_by_assigned[i], total_fps_per_target[i])
                detected_fps += fps_detected_by_assigned[i]
            elif ncob_count >= ecob_count:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                fps_detected_by_assigned[i] = ncob_count
                #fps_detected_by_assigned[i] = min(fps_detected_by_assigned[i], total_fps_per_target[i])
                detected_fps += fps_detected_by_assigned[i]
            else:
                assigned_metrics[i] = "ECOB"
                print(efx_count, ncob_count, ecob_count, N)
                count_ecob += 1
                fps_detected_by_assigned[i] = ecob_count
                #fps_detected_by_assigned[i] = min(fps_detected_by_assigned[i], total_fps_per_target[i])
                detected_fps += fps_detected_by_assigned[i]
    
    # Calculate effective efficiency
    effective_efficiency = (detected_fps / total_fps) * 100
    
    # Metric distribution percentages
    metric_counts = np.zeros(4)  # [EFX, SFX, NCOB, ECOB]
    metric_counts[0] = (count_efx / n_targets) * 100
    metric_counts[1] = (count_sfx / n_targets) * 100
    metric_counts[2] = (count_ncob / n_targets) * 100
    metric_counts[3] = (count_ecob / n_targets) * 100
    
    
    # Distribution of targets by N value
    n_counts = np.zeros(3)  # [N=0, N=1, N≥2]
    n_counts[0] = np.sum(n_bad == 0) / n_targets * 100
    n_counts[1] = np.sum(n_bad == 1) / n_targets * 100
    n_counts[2] = np.sum(n_bad >= 2) / n_targets * 100
    
    # Calculate detection efficiency by detection method
    method_efficiencies = np.zeros(5)  # [EFX, SFX, NCOB, ECOB, Overall]
    
    if np.sum(assigned_metrics == "EFX") > 0:
        method_efficiencies[0] = np.sum(fps_detected_by_assigned[assigned_metrics == "EFX"]) / np.sum(total_fps_per_target[assigned_metrics == "EFX"]) * 100
        print(np.where(fps_detected_by_assigned[assigned_metrics == "EFX"] > total_fps_per_target[assigned_metrics == "EFX"]))
    
    if np.sum(assigned_metrics == "SFX") > 0:
        method_efficiencies[1] = np.sum(fps_detected_by_assigned[assigned_metrics == "SFX"]) / np.sum(total_fps_per_target[assigned_metrics == "SFX"]) * 100
    
    if np.sum(assigned_metrics == "NCOB") > 0:
        method_efficiencies[2] = np.sum(fps_detected_by_assigned[assigned_metrics == "NCOB"]) / np.sum(total_fps_per_target[assigned_metrics == "NCOB"]) * 100
    
    if np.sum(assigned_metrics == "ECOB") > 0:
        method_efficiencies[3] = np.sum(fps_detected_by_assigned[assigned_metrics == "ECOB"]) / np.sum(total_fps_per_target[assigned_metrics == "ECOB"]) * 100
    
    method_efficiencies[4] = effective_efficiency  # Overall efficiency
    
    # Print results
    print(f"Effective efficiency: {effective_efficiency:.2f}%")
    print("\nMetric distribution:")
    print(f"  EFX: {metric_counts[0]:.2f}%")
    print(f"  SFX: {metric_counts[1]:.2f}%")
    print(f"  NCOB: {metric_counts[2]:.2f}%")
    print(f"  ECOB: {metric_counts[3]:.2f}%")
    
    print("\nTarget distribution by FP count:")
    print(f"  N=0: {n_counts[0]:.2f}%")
    print(f"  N=1: {n_counts[1]:.2f}%")
    print(f"  N≥2: {n_counts[2]:.2f}%")
    
    print("\nMethod efficiency for assigned targets:")
    print(f"  EFX: {method_efficiencies[0]:.2f}%")
    print(f"  SFX: {method_efficiencies[1]:.2f}%")
    print(f"  NCOB: {method_efficiencies[2]:.2f}%")
    print(f"  ECOB: {method_efficiencies[3]:.2f}%")
    print(f"  Overall: {method_efficiencies[4]:.2f}%")
    
    return (effective_efficiency, metric_counts, n_counts, method_efficiencies, total_fps, detected_fps, n_targets)

# Example call

results = calculate_effective_efficiency(
    data, data_sec, data_ext,
    eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
    eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
    fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
    flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
    delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first, 
    depth_sig_scaling
)