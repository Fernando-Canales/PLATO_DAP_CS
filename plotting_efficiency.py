import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from mpl_toolkits.axes_grid1.inset_locator import inset_axes #type: ignore
from matplotlib.ticker import FuncFormatter  # type: ignore

from imagette import ran_unique_int

dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/1000_targets_per_magnitude_bin/' 
#dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Fixed_transit_depths_and_durations/magnitude_bins/fixed_dback_132000ppm_and_td_1_422_hr/1000_targets_per_magnitude_bin/standard_results/'
#dataDIR = '/home/fercho//double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/EBs_rate/1000_targets_per_magnitude_bin/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/'
#dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/Long_Observational_Phase_Nord/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/'
#DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/EBs_rate/Distribution_transit_depth_and_durations/'
DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/all_contaminants_are_EBs/Distribution_transit_depth_and_durations/'
#DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/all_contaminants_are_EBs/Fixed_transit_depth_and_durations/'
#DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/EBs_rate/Distribution_transit_depth_and_durations/'
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
td_10first = data[:, 126:136]
dback_10first = data[:, 136:146]


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
eta_sec_24_cameras= data_sec[:, 6]
eta_sec_6_cameras = data_sec[:, 11]
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
delta_cob_ext = data_ext[:, 10]
eta_cob_ext = data_ext[:, 11]
sigma_cob_ext = data_ext[:, 12]

# We get now the 10 first values for each param. of the nominal cob shift
SPRK10_first = data[:, 17:27]
eta_cob_nom_10first_24_cameras = data[:, 46:56]
sigma_cob_10first_24_cameras = data[:, 56:66]
delta_cob_10first_24_cameras = data[:, 66:76]
eta_cob_nom_10first_6_cameras = data[:, 76:86]
sigma_cob_10first_6_cameras = data[:, 86:96]
delta_cob_10first_6_cameras = data[:, 96:106]
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

# Check for extreme values
extreme_mask = td_10first > 100  # More than 4 days
if extreme_mask.any():
    print(f"  Warning: {extreme_mask.sum()} values exceed 100 hours (max: {td_10first.max():.1f} hrs)")

# Create safe versions for calculations (avoiding division by zero)
# Use small value instead of zero to prevent NaN/Inf
MIN_TD = 0.1  # Minimum transit duration in hours (6 minutes)
td_10first_safe = np.where(td_10first > 0, td_10first, MIN_TD)

# We get now the new varaibles for the significant transit depth as Réza suggested
# Corrected significance depth calculations with safety:
sig_depth_secondary_mask_24_cameras = nsr1h_sec * (1 - data_sec[:, 5]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR*(1-SPRtot)
sig_depth_secondary_mask_6_cameras = nsr1h_sec_6_cameras * (1 - data_sec[:, 5]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR*(1-SPRtot)
sig_depth_extended_mask_24_cameras = data_ext[:, 4] * (1 - data_ext[:, 13]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR_1h^ext*(1-SPR_tot^ext) / sqrt(td ntr) [24 cameras]
sig_depth_extended_mask_6_cameras = data_ext[:, 44] * (1 - data_ext[:, 13]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR_1h^ext*(1-SPR_tot^ext) / sqrt(td ntr) [6 cameras]
sig_depth_nominal_mask_24_cameras = data[:, 7] * (1 - data[:, 11]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR_1h^nom*(1-SPR_tot^nom) / sqrt(td ntr) [24 cameras]
sig_depth_nominal_mask_6_cameras = data[:, 148] * (1 - data[:, 11]) / np.sqrt(td_10first_safe[:, 0] * ntr) # NSR_1h^nom*(1-SPR_tot^nom) / sqrt(td ntr) [6 cameras]

# We obtainn the value of the quadratic sum of the noises (Eq. (38) of the paper)
sig_depth_24_cameras = np.sqrt(sig_depth_nominal_mask_24_cameras**2 + sig_depth_extended_mask_24_cameras**2)
sig_depth_6_cameras = np.sqrt(sig_depth_nominal_mask_6_cameras**2 + sig_depth_extended_mask_6_cameras**2)

# First, the expressions for the efficiency false positives detections for the different masks
fp_single_contaminant_24_cameras = (eta_t > flux_thresh_nom_mask)  # false positive
sdr_flux = fp_single_contaminant_24_cameras & (eta_sec_24_cameras> flux_thresh_nom_mask) & (delta_obs_c > delta_obs_t)  # secondary mask false positive detection rate
fp_single_contaminant_6_cameras = (eta_t_6_cameras > flux_thresh_nom_mask)
secondary_mask_conditions_24_cameras = (eta_sec_24_cameras> flux_thresh_sec_mask) & (delta_obs_c > delta_obs_t  + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_24_cameras**2 + sig_depth_nominal_mask_24_cameras**2)) & fp_single_contaminant_24_cameras  # secondary mask efficiency condition for 24 cameras
secondary_mask_conditions_6_cameras = (eta_sec_6_cameras > flux_thresh_sec_mask) & (delta_obs_c > delta_obs_t + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_6_cameras**2 + sig_depth_nominal_mask_6_cameras**2)) & fp_single_contaminant_6_cameras # secondary mask efficiency condition for 6 cameras
eficiency_extended_mask_highest_spr_contaminant = fp_single_contaminant_24_cameras & (eta_ext_for_the_most_significant_contaminant > flux_thresh_ext_mask) & (delta_obs_ext_single_contaminant > delta_obs_t + depth_sig_scaling*sig_depth_24_cameras)  # extended mask false positive detection rate
secondary_mask_conditions_cob_24_cameras = (eta_cob_sec_24_cameras > cob_thresh) & fp_single_contaminant_24_cameras &  (delta_cob_sec > 10*sigma_cob_sec_24_cameras/np.sqrt(td_10first_safe[:, 0] * ntr)) # secondary mask efficiency condition for 24 cameras and cob shift
secondary_mask_conditions_cob_6_cameras = (eta_cob_sec_6_cameras > cob_thresh)  & fp_single_contaminant_6_cameras  &  (delta_cob_sec_6_cameras > 10*sigma_cob_sec_6_cameras/np.sqrt(td_10first_safe[:, 0] * ntr)) # secondary mask efficiency condition for 6 cameras and cob shift

# We reshape some of them for 24 cameras
# For 10-first arrays with broadcasting (using safe version):
sig_depth_extended_mask_24_cameras_10first = data_ext[:, 4][:, None] * (1 - data_ext[:, 13][:, None]) / np.sqrt(td_10first_safe * ntr)
sig_depth_extended_mask_6_cameras_10first = data_ext[:, 44][:, None] * (1 - data_ext[:, 13][:, None]) / np.sqrt(td_10first_safe * ntr)
sig_depth_nominal_mask_24_cameras_10first = data[:, 7][:, None] * (1 - data[:, 11][:, None]) / np.sqrt(td_10first_safe * ntr)
sig_depth_nominal_mask_6_cameras_10first = data[:, 148][:, None] * (1 - data[:, 11][:, None]) / np.sqrt(td_10first_safe * ntr)


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
    dback= dback_10first[i, ::] #np.ones(10)*dback_ref
    #td = td_from_dap[i, ::] #np.ones(10)*td_ref
    td = td_10first_safe[i, :]
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

# Set consistent figure size for both plots
fig_width = 8
fig_height = 6
plt.figure(4, figsize=(fig_width, fig_height))
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

    delta_eff_ext_flux = eff_ext_overall_24_cameras - eff_ext_overall_6_cameras
    combined_error_ext_flux = np.sqrt(error**2 + error_6_cameras**2)
    test_significance_ext_flux = delta_eff_ext_flux/combined_error_ext_flux

    delta_eff_sec_flux = eff_sec - eff_sec_6_cameras
    combined_error_sec_flux = np.sqrt(error_sec**2 + error_sec_6_cameras**2)
    test_significance_sec_flux = delta_eff_sec_flux/combined_error_sec_flux


    #print("Diagnostic extended flux:", test_significance_ext_flux)
    #print("Diagnostic secondary flux:", delta_eff_sec_flux)

    # Plotting errorbars with labels only once
    plt.errorbar(Pi, eff_sec, yerr=error_sec, fmt='o', color='purple', ecolor='purple', capsize=5, label='Sec. Mask (24 cameras)' if not labels_added['sec_24'] else "", markersize=4)
    plt.errorbar(Pi, eff_sec_6_cameras, yerr=error_sec_6_cameras, fmt='o', color='green', ecolor='green', capsize=5, label='Sec. Mask (6 cameras)' if not labels_added['sec_6'] else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_24_cameras, yerr=error, fmt='s', color='blue', ecolor='blue', capsize=5, label='Ext. Mask (24 cameras)' if not labels_added['ext_24'] else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_6_cameras, yerr=error_6_cameras, fmt='s', color='red', ecolor='red', capsize=5, label='Ext. Mask (6 cameras)' if not labels_added['ext_6'] else "", markersize=4)
    #plt.fill_between([9, 11.7], [47, 47], [100, 100], color='aqua', alpha=0.1) # type: ignore 
    #plt.fill_between([11, 13.4], [47, 47], [100, 100], color='plum', alpha=0.1) # type: ignore
    plt.fill_between([9, 11.7], [40, 40], [100, 100], color='aqua', alpha=0.1) # type: ignore this is for the variabel transit parameters case
    plt.fill_between([10.7, 13.4], [40, 40], [100, 100], color='plum', alpha=0.1) # type: ignore this is for the variabel transit parameters case

    # Update label tracking
    labels_added['sec_24'] = True
    labels_added['sec_6'] = True
    labels_added['ext_24'] = True
    labels_added['ext_6'] = True
    labels_added['ext_2_pix_24'] = True
    labels_added['ext_2_pix_6'] = True

    # Connecting lines for the previous data points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_sec, eff_sec], color='purple', linestyle='-', markersize=0) # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_sec_6_cameras, eff_sec_6_cameras], color='green', linestyle='-', markersize=0) # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall, eff_ext_overall_24_cameras], color='blue', linestyle='-', markersize=0) # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall_6_cameras, eff_ext_overall_6_cameras], color='red', linestyle='-', markersize=0) # type: ignore
    
    # Update previous values
    prev_Pi, prev_eff_sec, prev_eff_sec_6_cameras, prev_eff_ext_overall, prev_eff_ext_overall_6_cameras = Pi, eff_sec, eff_sec_6_cameras, eff_ext_overall_24_cameras, eff_ext_overall_6_cameras

# Additional plot settings
#plt.vlines(11.7, ymin=47, ymax=100, linestyles='dashed', colors='green') # type: ignore
#plt.vlines(11, ymin=47, ymax=100, linestyles='dashdot', colors='red') # type: ignore
plt.vlines(11.7, ymin=40, ymax=100, linestyles='dashed', colors='green') # type: ignore this is for the variabel transit parameters case
plt.vlines(10.7, ymin=40, ymax=100, linestyles='dashdot', colors='red') # type: ignore this is for the variabel transit parameters case
#plt.ylim(55, 100)
plt.ylim(40, 100) # this is for the variabel transit parameters case
plt.xlim(9.9, 13.1)
#plt.text(10, 78.1, 'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold')
plt.text(10, 75.1, 'Earth-like planet detection\nregion (24 cameras)',  color='green', weight='bold') # This line is fo the variable transit parameters case
plt.text(11.2, 63, 'On-board light curve processing region', color='red', weight='bold')

# Display legend below the plot
plt.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', borderaxespad=0., ncol=2)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)

# Adjust layout to accommodate the legend below
plt.tight_layout(rect=[0, 0, 1, 0.88])
#plt.savefig(DIRout + "DAP_Flux_efficiency_standard.pdf", format='pdf', bbox_inches='tight')
plt.savefig(DIRout + "DAP_Flux_efficiency_variable_standard_results.pdf", format='pdf', bbox_inches='tight') # this is for the variable transit parameters case
plt.show()


def analyze_centroid_issues(delta_cob_nom, delta_cob_ext):
    # 1. Check for zeros/NaNs
    print(f"Nominal shifts: {np.sum(delta_cob_nom == 0)} zeros, {np.sum(np.isnan(delta_cob_nom))} NaNs")
    print(f"Extended shifts: {np.sum(delta_cob_ext == 0)} zeros, {np.sum(np.isnan(delta_cob_ext))} NaNs")

# 2. Collect data for contaminants producing zero centroid shifts
    zero_nom_distances = []
    zero_nom_cont_mags = []
    zero_nom_target_mags = []
    zero_nom_mag_diffs = []
    
    zero_ext_distances = []
    zero_ext_cont_mags = []
    zero_ext_target_mags = []
    zero_ext_mag_diffs = []
    
    for i in range(len(dist_target_to_10first_contaminants)):
        for j in range(10):
            dist = dist_target_to_10first_contaminants[i,j]
            shift_nom = delta_cob_10first_24_cameras[i,j]
            shift_ext = delta_cob_ext_10first_24_cameras[i,j]
            target_mag = mag[i]
            cont_mag = mag_target_10first_contaminants[i,j]
            mag_diff = cont_mag - target_mag
            
            if shift_nom == 0:
                zero_nom_distances.append(dist)
                zero_nom_cont_mags.append(cont_mag)
                zero_nom_target_mags.append(target_mag)
                zero_nom_mag_diffs.append(mag_diff)
            
            if shift_ext == 0:
                zero_ext_distances.append(dist)
                zero_ext_cont_mags.append(cont_mag)
                zero_ext_target_mags.append(target_mag)
                zero_ext_mag_diffs.append(mag_diff)
    
    # Convert to numpy arrays
    zero_nom_distances = np.array(zero_nom_distances)
    zero_nom_cont_mags = np.array(zero_nom_cont_mags)
    zero_nom_mag_diffs = np.array(zero_nom_mag_diffs)
    
    zero_ext_distances = np.array(zero_ext_distances)
    zero_ext_cont_mags = np.array(zero_ext_cont_mags)
    zero_ext_mag_diffs = np.array(zero_ext_mag_diffs)
    
    print(f"\nNominal zero shifts:\n - Distance stats: mean={np.mean(zero_nom_distances):.2f}, std={np.std(zero_nom_distances):.2f}")
    print(f"Nominal zero shifts:\n - Mag diff stats: mean={np.mean(zero_nom_mag_diffs):.2f}, std={np.std(zero_nom_mag_diffs):.2f}")
    print(f"Extended zero shifts\n - Distance stats: mean={np.mean(zero_ext_distances):.2f}, std={np.std(zero_ext_distances):.2f}")
    print(f"Extended zero shifts\n - Mag diff stats: mean={np.mean(zero_ext_mag_diffs):.2f}, std={np.std(zero_ext_mag_diffs):.2f}")
    
    return zero_nom_distances, zero_nom_mag_diffs, zero_ext_distances, zero_ext_mag_diffs

# Unpack the returned tuple
nom_distances, nom_mag_diffs, ext_distances, ext_mag_diffs = analyze_centroid_issues(delta_cob_10first_24_cameras, delta_cob_ext_10first_24_cameras) # type: ignore

# Now use the unpacked variables directly
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
plt.hist(nom_distances, bins=30, alpha=0.7, label='Nominal')
plt.hist(ext_distances, bins=30, alpha=0.7, label='Extended')
plt.xlabel('Distance (px)')
plt.ylabel('Count')
plt.title('Distance Distribution for Zero Shifts')
plt.legend()

plt.subplot(2, 2, 2)
plt.hist(nom_mag_diffs, bins=30, alpha=0.7, label='Nominal')
plt.hist(ext_mag_diffs, bins=30, alpha=0.7, label='Extended')
plt.xlabel('Magnitude Difference')
plt.ylabel('Count')
plt.title('Mag Difference Distribution for Zero Shifts')
plt.legend()

plt.tight_layout()
plt.show()

"""
Now we plot the efficiency of the COB shift (all contaminants)
    
"""
# Initialize lists to store data points for the inset
Pi_values = []
eff_ext_cob_overall_values = []
eff_ext_cob_overall_6_cameras_values = []
eff_cob_values = []
eff_cob_6_cameras_values = []

plt.figure(5, figsize=(fig_width, fig_height))
for i in range(nP):
    Pi = Pmin + i * binsize
    Pi_values.append(Pi)  # Store Pi for later

    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    s_24_cameras = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)[m,:].sum()
    s_6_cameras = (eta_nom_bt_6_cameras > flux_thresh_nom_mask)[m,:].sum()
    eff_ext_cob_overall = ((eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (delta_cob_ext_10first_24_cameras > 10*sigma_cob_ext_10first_24_cameras/ np.sqrt(td_10first_safe * ntr)))[m,:].sum() / s_24_cameras * 100.
    eff_ext_cob_overall_6_cameras = ((eta_cob_ext_10first_6_cameras > cob_thresh) &  (eta_nom_bt_6_cameras > flux_thresh_nom_mask) & (delta_cob_ext_10first_6_cameras > 10*sigma_cob_ext_10first_6_cameras/ np.sqrt(td_10first_safe * ntr)))[m,:].sum() / s_6_cameras * 100.
    eff_cob = ((eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (delta_cob_10first_24_cameras > 10*sigma_cob_10first_24_cameras/ np.sqrt(td_10first_safe * ntr)))[m,:].sum() / s_24_cameras * 100.
    eff_cob_6_cameras = ((eta_cob_nom_10first_6_cameras > cob_thresh) & (eta_nom_bt_6_cameras > flux_thresh_nom_mask) & (delta_cob_10first_6_cameras > 10*sigma_cob_10first_6_cameras/ np.sqrt(td_10first_safe * ntr)))[m,:].sum() / s_6_cameras * 100.
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

    delta_eff_sec_cob = eff_cob_sec - eff_cob_sec_6_cameras
    combined_error_sec_cob = np.sqrt(error_cob_sec**2 + error_cob_sec_6_cameras**2)
    test_significance_sec_cob = delta_eff_sec_cob/combined_error_sec_cob
    
    plt.errorbar(Pi, eff_ext_cob_overall, fmt='s', yerr=error_ext_cob, label='Ext. Mask (24 cameras)' if i == 0 else "", color='blue', ecolor='blue', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_ext_cob_overall_6_cameras, fmt='s', yerr=error_ext_cob_6_cameras, label='Ext. Mask (6 cameras)' if i == 0 else "", color='red', ecolor='red', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob, fmt='*', yerr=error_cob, label='Nom. Mask (24 cameras)' if i == 0 else "", color='orange', ecolor='orange', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_6_cameras, fmt='*', yerr=error_cob_6_cameras, label='Nom. Mask (6 cameras)' if i == 0 else "", color='olive', ecolor='olive', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec, fmt='o', yerr=error_cob_sec, label='Sec. Mask (24 cameras)' if i == 0 else "", color='purple', ecolor='purple', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec_6_cameras, fmt='o', yerr=error_cob_sec_6_cameras, label='Sec. Mask (6 cameras)' if i == 0 else "", color='green', ecolor='green', capsize=5, markersize=4)
    #plt.fill_between([9, 11.7], [47, 47], [100, 100], color='aqua', alpha=0.1) # type: ignore
    #plt.fill_between([11, 13.4], [47,47], [100, 100], color='plum', alpha=0.1) # type: ignore
    plt.fill_between([9, 11.7], [40, 40], [100, 100], color='aqua', alpha=0.1) # type: ignore this is for the variable transit parameters case
    plt.fill_between([10.7, 13.4], [40,40], [100, 100], color='plum', alpha=0.1) # type: ignore this is for the variable transit parameters case
    
    # Plot lines connecting the points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall, eff_ext_cob_overall], color='blue', linestyle='-') # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_cob, eff_cob], color='orange', linestyle='-') # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec, eff_cob_sec], color='purple', linestyle='-') # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall_6_cameras, eff_ext_cob_overall_6_cameras], color='red', linestyle='-') # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_cob_6_cameras, eff_cob_6_cameras], color='olive', linestyle='-') # type: ignore
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec_6_cameras, eff_cob_sec_6_cameras], color='green', linestyle='-')    # type: ignore
    # Update previous values
    prev_Pi, prev_eff_ext_cob_overall, prev_eff_cob, prev_eff_cob_sec, prev_eff_ext_cob_overall_6_cameras, prev_eff_cob_6_cameras, prev_eff_cob_sec_6_cameras = Pi, eff_ext_cob_overall, eff_cob, eff_cob_sec, eff_ext_cob_overall_6_cameras, eff_cob_6_cameras, eff_cob_sec_6_cameras


    #plt.vlines(11.7, ymin=47, ymax = 100, linestyles='dashed', colors='green') # type: ignore
    #plt.vlines(11, ymin=47, ymax=100, linestyles='dashdot', colors='red') # type: ignore
    plt.vlines(11.7, ymin=40, ymax = 100, linestyles='dashed', colors='green') # type: ignore this is for the variable transit parameters case
    plt.vlines(10.7, ymin=40, ymax=100, linestyles='dashdot', colors='red') # type: ignore this is for the variable transit parameters case
    #plt.ylim(47, 100)
    plt.ylim(40, 100) # this is for the variable transit parameters case
    plt.xlim(9.9, 13.1)
    # Move text down a bit to make room for inset if inside plot
    plt.text(10., 61.5,'Earth-like planet detection \nregion (24 cameras)', color='green', weight='bold')
    plt.text(10.9, 50.5, 'On-board light curve process. region', color='red', weight='bold')
           
# Finalize the main plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.21), borderaxespad=0., fancybox=True, ncol=3, columnspacing=0.6)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)
plt.tight_layout(rect=[0, 0, 1, 0.88])


#plt.savefig(DIRout + "DAP_CS_efficiency_standard.pdf", format='pdf', bbox_inches='tight')
plt.savefig(DIRout + "DAP_CS_efficiency_variable_standard_results.pdf", format='pdf', bbox_inches='tight') # this is for the variable 
plt.show()

# Flatten centroid arrays and repeat corresponding arrays to match
eta_cob_flat = eta_cob_nom_10first_24_cameras.flatten()
eta_cob_ext_flat = eta_cob_ext_10first_24_cameras.flatten()

# Repeat magnitude and secondary flux arrays to match flattened centroid arrays
mag_repeated = np.repeat(mag, 10)
eta_sec_24_cameras_repeated = np.repeat(eta_sec_24_cameras, 10)

# Create the comparison plot
plt.figure(6, figsize=(12, 8))

# Remove invalid values
# Compute valid mask (same as before)
valid_mask = (eta_sec_24_cameras_repeated > 0) & (eta_cob_flat > 0) & (eta_cob_ext_flat > 0)

# Compute log ratios
log_ratio_nom = np.log10(eta_sec_24_cameras_repeated[valid_mask] / eta_cob_flat[valid_mask])
log_ratio_ext = np.log10(eta_sec_24_cameras_repeated[valid_mask] / eta_cob_ext_flat[valid_mask])

plt.subplot(1, 2, 1)
plt.hist(log_ratio_nom, bins=50, color='red', alpha=0.7)
plt.axvline(0, color='black', linestyle='--', linewidth=1)  # reference line
plt.xlabel(r'$\log_{10}(\eta_{\mathrm{k_{max}}}^{\mathrm{sec}} / \eta_{k}^{\mathrm{nom}, \Delta C})$', fontsize=fsize)
plt.ylabel('Counts', fontsize=fsize)


plt.subplot(1, 2, 2)
plt.hist(log_ratio_ext, bins=50, color='green', alpha=0.7)
plt.axvline(0, color='black', linestyle='--', linewidth=1)  # reference line
plt.xlabel(r'$\log_{10}(\eta_{\mathrm{k_{max}}}^{\mathrm{sec}} / \eta_{k}^{\mathrm{ext}, \Delta C})$', fontsize=fsize)
plt.ylabel('Counts', fontsize=fsize)

plt.tight_layout()
plt.savefig(DIRout + "Hist_log_ratios_eta_sec_vs_nom_and_ext_cob.pdf", format='pdf', bbox_inches='tight')
plt.show()




nfp = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_mask = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) & (delta_obs_ext > delta_obs + depth_sig_scaling*sig_depth_24_cameras_10first) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_flux_without_significant_transit_depth_condition = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) &  (delta_obs_ext > delta_obs) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_mask_single_contaminant = (eficiency_extended_mask_highest_spr_contaminant)
nfp_nom_cob = (eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (delta_cob_10first_24_cameras > 10*sigma_cob_10first_24_cameras/np.sqrt(td_10first_safe * ntr))
nfp_ext_cob = (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (delta_cob_ext_10first_24_cameras > 10*sigma_cob_ext_10first_24_cameras/np.sqrt(td_10first_safe * ntr))


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
weighted_fraction_fp_nom_cob_no_ext_cob = 0
weighted_fraction_fp_sec_cob_no_ext_cob = 0
weighted_fraction_fp_ext_cob_no_sec_cob = 0
weighted_fraction_fp_ext_flux_no_sec_cob = 0


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
weighted_variance_fraction_fp_ext_flux_no_sec_cob = 0
weighted_variance_fraction_fp_nom_cob_no_ext_cob = 0


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
    n_star[i] = m.sum() # type: ignore

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
    fraction_fp_nom_cob_no_ext_cob = ((nfp_ext_cob[m]==False) & nfp_nom_cob[m] & nfp[m]).sum()/nfp[m].sum()
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
    fraction_fp_ext_flux_no_sec_cob = ((secondary_mask_conditions_cob_24_cameras[m] == False) & nfp_ext_mask[m][:, 0] & fp_single_contaminant_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum()




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
    weighted_fraction_fp_ext_flux_no_sec_cob += weights[i] * fraction_fp_ext_flux_no_sec_cob
    weighted_fraction_fp_nom_cob_no_ext_cob += weights[i] * fraction_fp_nom_cob_no_ext_cob


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
    variance_fraction_fp_ext_flux_no_sec_cob = (fraction_fp_ext_flux_no_sec_cob * (1 - fraction_fp_ext_flux_no_sec_cob)) / n_star[i]
    variance_fraction_fp_nom_cob_no_ext_cob = (fraction_fp_nom_cob_no_ext_cob * (1 - fraction_fp_nom_cob_no_ext_cob)) / n_star[i]

     

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
    weighted_variance_fraction_fp_ext_flux_no_sec_cob += weights[i] * variance_fraction_fp_ext_flux_no_sec_cob
    weighted_variance_fraction_fp_nom_cob_no_ext_cob += weights[i] * variance_fraction_fp_nom_cob_no_ext_cob




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
weighted_error_fraction_fp_ext_flux_no_sec_cob = np.sqrt(weighted_variance_fraction_fp_ext_flux_no_sec_cob)
weighted_error_fraction_fp_nom_cob_no_ext_cob = np.sqrt(weighted_variance_fraction_fp_nom_cob_no_ext_cob)

# Print weighted efficiency results with errors
print(f'Weighted extended flux efficiency: {weighted_eff_ext_flux:.2f}% ± {weighted_error_eff_ext_flux:.2f}%')
print(f'Weighted secondary flux efficiency: {weighted_eff_sec_flux:.2f}% ± {weighted_error_eff_sec_flux:.2f}%')
print(f'Weighted extended flux for a single contaminant efficiency: {weighted_eff_ext_flux_single_contaminant:.2f}% ± {weighted_error_eff_ext_flux_single_contaminant:.2f}%')
print(f'Weighted extended flux withouth the significant transit condition: {weighted_eff_ext_flux_without_significant_transit_depth_condition:.2f}% ± {weighted_error_eff_ext_flux_without_significant_transit_depth_condition:.2f}%')
print(f'Weighted nominal COB efficiency: {weighted_eff_nom_cob:.2f}% ± {weighted_error_eff_nom_cob:.2f}%')
print(f'Weighted extended COB efficiency: {weighted_eff_ext_cob:.2f}% ± {weighted_error_eff_ext_cob:.2f}%')
print(f'Weighted secondary mask COB efficiency: {weighted_eff_sec_cob:.2f}% ± {weighted_error_eff_sec_cob:.2f}%')

print(f'Weighted fraction of FPs detected by NCOB but not by SCOB: {weighted_fraction_fp_nom_cob_no_sec_cob * 100:.2f}% ± {weighted_error_fraction_fp_nom_cob_no_sec_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SCOB but not by ECOB: {weighted_fraction_fp_sec_cob_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_cob_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by SCOB: {weighted_fraction_fp_ext_cob_no_sec_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_sec_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SFX but not by ECOB: {weighted_fraction_fp_sec_flux_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_flux_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by SFX: {weighted_fraction_fp_ext_cob_no_sec_flux * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_sec_flux * 100:.2f}%')

print(f'Weighted fraction of FPs detected by ECOB but not by EFX: {weighted_fraction_fp_ext_cob_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by NCOB but not by EFX: {weighted_fraction_fp_nom_cob_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_nom_cob_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by NCOB: {weighted_fraction_fp_ext_cob_no_nom_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_nom_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by EFX but not by ECOB: {weighted_fraction_fp_ext_flux_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_flux_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SFX but not by EFX: {weighted_fraction_fp_sec_flux_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_sec_flux_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by EFX but not by SFX: {weighted_fraction_fp_ext_flux_no_sec_flux * 100:.2f}% ± {weighted_error_fraction_fp_ext_flux_no_sec_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SFX but not by NCOB: {weighted_fraction_fp_sec_flux_no_nom_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_flux_no_nom_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by NCOB but not by SFX: {weighted_fraction_fp_nom_cob_no_sec_flux * 100:.2f}% ± {weighted_error_fraction_fp_nom_cob_no_sec_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SFX but not by SCOB: {weighted_fraction_fp_sec_flux_no_sec_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_flux_no_sec_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SCOB but not by SFX: {weighted_fraction_fp_sec_cob_no_sec_flux * 100:.2f}% ± {weighted_error_fraction_fp_sec_cob_no_sec_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SCOB but not by EFX: {weighted_fraction_fp_sec_cob_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_sec_cob_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SCOB but not by NCOB: {weighted_fraction_fp_sec_cob_no_nom_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_cob_no_nom_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by NCOB but not by ECOB: {weighted_fraction_fp_nom_cob_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_nom_cob_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SCOB but not by ECOB: {weighted_fraction_fp_sec_cob_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_cob_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by SCOB: {weighted_fraction_fp_ext_cob_no_sec_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_sec_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by SFX but not by ECOB: {weighted_fraction_fp_sec_flux_no_ext_cob * 100:.2f}% ± {weighted_error_fraction_fp_sec_flux_no_ext_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by SFX: {weighted_fraction_fp_ext_cob_no_sec_flux * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_sec_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by EFX but not by SCOB: {weighted_fraction_fp_ext_flux_no_sec_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_flux_no_sec_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by EFX but not by NCOB: {weighted_fraction_fp_ext_flux_no_nom_cob * 100:.1f}% ± {weighted_error_fraction_fp_ext_flux_no_nom_cob * 100:.1f}%')
 