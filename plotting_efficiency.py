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
secondary_mask_conditions_cob_24_cameras = (eta_cob_sec_24_cameras > cob_thresh) & (delta_cob_sec > 10*sigma_cob_sec_24_cameras) & fp_single_contaminant_24_cameras  # secondary mask efficiency condition for 24 cameras and cob shift
secondary_mask_conditions_cob_6_cameras = (eta_cob_sec_6_cameras > cob_thresh) &  (delta_cob_sec_6_cameras > 10*sigma_cob_sec_6_cameras) & fp_single_contaminant_6_cameras # secondary mask efficiency condition for 6 cameras and cob shift
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
        # Plotting errorbars with labels only once
    plt.errorbar(Pi, eff_sec, yerr=error_sec, fmt='o', color='purple', ecolor='purple', capsize=5, label='Sec. Mask (24 cameras)' if not labels_added['sec_24'] else "", markersize=4)
    plt.errorbar(Pi, eff_sec_6_cameras, yerr=error_sec_6_cameras, fmt='o', color='green', ecolor='green', capsize=5, label='Sec. Mask (6 cameras)' if not labels_added['sec_6'] else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_24_cameras, yerr=error, fmt='s', color='blue', ecolor='blue', capsize=5, label='Ext. Mask (24 cameras)' if not labels_added['ext_24'] else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_6_cameras, yerr=error_6_cameras, fmt='s', color='red', ecolor='red', capsize=5, label='Ext. Mask (6 cameras)' if not labels_added['ext_6'] else "", markersize=4)
    #plt.fill_between([9, 11.7], [47, 47], [100, 100], color='aqua', alpha=0.1) # type: ignore 
    #plt.fill_between([11, 13.4], [47, 47], [100, 100], color='plum', alpha=0.1) # type: ignore
    plt.fill_between([9, 11.7], [40, 40], [100, 100], color='aqua', alpha=0.1) # type: ignore this is for the variabel transit parameters case
    plt.fill_between([11, 13.4], [40, 40], [100, 100], color='plum', alpha=0.1) # type: ignore this is for the variabel transit parameters case

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
plt.vlines(11, ymin=40, ymax=100, linestyles='dashdot', colors='red') # type: ignore this is for the variabel transit parameters case
#plt.ylim(55, 100)
plt.ylim(40, 100) # this is for the variabel transit parameters case
plt.xlim(9.9, 13.1)
#plt.text(10, 78.1, 'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold')
plt.text(10, 72.1, 'Earth-like planet detection\nregion (24 cameras)',  color='green', weight='bold') # This line is fo the variable transit parameters case
plt.text(11.2, 50, 'On-board light curve processing region', color='red', weight='bold')

# Display legend below the plot
plt.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', borderaxespad=0., ncol=2)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)

# Adjust layout to accommodate the legend below
plt.tight_layout(rect=[0, 0, 1, 0.88])
#plt.savefig(DIRout + "DAP_Flux_efficiency_standard.pdf", format='pdf', bbox_inches='tight')
plt.savefig(DIRout + "DAP_Flux_efficiency_variable.pdf", format='pdf', bbox_inches='tight') # this is for the variable transit parameters case
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

plt.figure(5)
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
    
    plt.errorbar(Pi, eff_ext_cob_overall, fmt='s', yerr=error_ext_cob, label='Ext. Mask (24 cameras)' if i == 0 else "", color='blue', ecolor='blue', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_ext_cob_overall_6_cameras, fmt='s', yerr=error_ext_cob_6_cameras, label='Ext. Mask (6 cameras)' if i == 0 else "", color='red', ecolor='red', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob, fmt='*', yerr=error_cob, label='Nom. Mask (24 cameras)' if i == 0 else "", color='orange', ecolor='orange', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_6_cameras, fmt='*', yerr=error_cob_6_cameras, label='Nom. Mask (6 cameras)' if i == 0 else "", color='olive', ecolor='olive', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec, fmt='o', yerr=error_cob_sec, label='Sec. Mask (24 cameras)' if i == 0 else "", color='purple', ecolor='purple', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec_6_cameras, fmt='o', yerr=error_cob_sec_6_cameras, label='Sec. Mask (6 cameras)' if i == 0 else "", color='green', ecolor='green', capsize=5, markersize=4)
    #plt.fill_between([9, 11.7], [47, 47], [100, 100], color='aqua', alpha=0.1) # type: ignore
    #plt.fill_between([11, 13.4], [47,47], [100, 100], color='plum', alpha=0.1) # type: ignore
    plt.fill_between([9, 11.7], [40, 40], [100, 100], color='aqua', alpha=0.1) # type: ignore this is for the variable transit parameters case
    plt.fill_between([11, 13.4], [40,40], [100, 100], color='plum', alpha=0.1) # type: ignore this is for the variable transit parameters case
    
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
    plt.vlines(11, ymin=40, ymax=100, linestyles='dashdot', colors='red') # type: ignore this is for the variable transit parameters case
    #plt.ylim(47, 100)
    plt.ylim(40, 100) # this is for the variable transit parameters case
    plt.xlim(9.9, 13.1)
    # Move text down a bit to make room for inset if inside plot
    plt.text(10., 71.5,'Earth-like planet detection \nregion (24 cameras)', color='green', weight='bold')
    plt.text(11.1, 60.5, 'On-board light curve process. region', color='red', weight='bold')
           
# Finalize the main plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.21), borderaxespad=0., fancybox=True, ncol=3, columnspacing=0.6)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)
plt.tight_layout(rect=[0, 0, 1, 0.88])
# After loop: plot inset with accumulated data
#ax_inset = inset_axes(plt.gca(), width="35%", height="35%", loc="center")
ax_inset = inset_axes(
    plt.gca(), 
    width=1.9,    # Width in inches
    height=0.5,   # Height in inches
    loc="center", 
    bbox_to_anchor=(0.48, 0.8),  # Center it horizontally (0.5) and position lower vertically (0.3)
    bbox_transform=plt.gca().transAxes
)
ax_inset.errorbar(Pi_values, eff_ext_cob_overall_values, fmt='s', yerr=error_ext_cob, color='blue', linewidth=1.5, label='Ext. Mask (24 cameras)') # type: ignore
ax_inset.errorbar(Pi_values, eff_ext_cob_overall_6_cameras_values, fmt='s', yerr=error_ext_cob_6_cameras, color='red', linewidth=1.5, label='Ext. Mask (6 cameras)') # type: ignore
ax_inset.errorbar(Pi_values, eff_cob_values, color='orange', fmt='*', yerr=error_cob, linewidth=1.5, label='Nom. Mask (24 cameras)') # type: ignore
ax_inset.errorbar(Pi_values, eff_cob_6_cameras_values, color='olive', fmt='*', yerr=error_cob_6_cameras,  linewidth=1.5, label='Nom. Mask (6 cameras)') # type: ignore

# Connect the dots
ax_inset.plot(Pi_values, eff_ext_cob_overall_values, color='blue', linestyle='-')
ax_inset.plot(Pi_values, eff_cob_values, color='orange', linestyle='-')
ax_inset.plot(Pi_values, eff_ext_cob_overall_6_cameras_values, color='red', linestyle='-')
ax_inset.plot(Pi_values, eff_cob_6_cameras_values, color='olive', linestyle='-')

# Add shaded areas and vertical lines to the inset plot
ax_inset.fill_between([9, 11.7], [95, 95], [100, 100], color='aqua', alpha=0.5)
ax_inset.fill_between([11, 13.4], [95, 95], [100, 100], color='plum', alpha=0.5)
ax_inset.vlines(11.7, ymin=95, ymax=100, linestyles='dashed', colors='green')
ax_inset.vlines(11, ymin=95, ymax=100, linestyles='dashdot', colors='red')

# Set inset plot limits and appearance
ax_inset.set_ylim(95, 100)
ax_inset.set_xlim(9.9, 13.1)
ax_inset.set_yticks([95, 97, 100])
ax_inset.grid(True, linestyle='--', linewidth=0.5, color='grey', alpha=0.5)
plt.gca().indicate_inset_zoom(ax_inset, edgecolor="grey")

#plt.savefig(DIRout + "DAP_CS_efficiency_standard.pdf", format='pdf', bbox_inches='tight')
plt.savefig(DIRout + "DAP_CS_efficiency_variable.pdf", format='pdf', bbox_inches='tight') # this is for the variable 
plt.show()

    

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




def calculate_effective_efficiency(data, data_sec, data_ext, 
                                  eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
                                  eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
                                  fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
                                  flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
                                  delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first,
                                  depth_sig_scaling, Nmax_threshold):
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
    n_targets = len(mag)
    
    
    # Arrays to store assigned metrics and detected FPs
    assigned_metrics = np.empty(n_targets, dtype='U4')  # EFX, SFX, NCOB, ECOB
    fps_detected_by_assigned = np.zeros(n_targets)
    total_fps_per_target = np.zeros(n_targets)
    
# variables to count metric assignments
    count_efx_zero = 0  # Count of EFX assigned for N=0 cases
    count_efx_multi = 0  # Count of EFX assigned for N≥2 cases
    count_sfx = 0
    count_ncob = 0
    count_ncob_n1 = 0  # NCOB assigned to N=1 targets
    count_ncob_n2plus = 0  # NCOB assigned to N≥2 targets
    count_ecob = 0
    count_none = 0  # For cases with N > threshold
    
    # variables to track FPs
    total_fps = 0 # Total count of FPs across all targets
    detected_fps = 0 # Total count of detected FPs across all targets
    fps_in_pathological = 0 
    
    # Define detection capability for each method (for each target)
    
    # Calculate detection capability for each method
    efx_detection = (eta_ext_bt_24_cameras > flux_thresh_ext_mask) & \
                    (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    
    sfx_detection = secondary_mask_conditions_24_cameras
    
    ncob_detection = (eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    
    ecob_detection = (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
    # For each target, follow the decision flow to assign the best metric
    for i in range(n_targets):
        # Count the number of potential FPs for this target
        N = int(n_bad[i])
        total_fps_per_target[i] = N
        total_fps += N
        if N > Nmax_threshold:
            assigned_metrics[i] = "NONE"
            count_none += 1
            fps_detected_by_assigned[i] = 0
            fps_in_pathological += N

        # Implement the decision logic from the flow diagram
        elif N == 0:
            # Assign Extended Flux (Low Priority)
            assigned_metrics[i] = "EFX"
            count_efx_zero += 1
            # No FPs to detect
            fps_detected_by_assigned[i] = 0
            
        elif N == 1:
            # Check if SFX can detect the highest SPR contaminant
            if sfx_detection[i]:
                assigned_metrics[i] = "SFX"
                count_sfx += 1
                fps_detected_by_assigned[i] = 1
                detected_fps += fps_detected_by_assigned[i]
            else:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                count_ncob_n1 += 1  # Count NCOB for N=1 case
                fps_detected_by_assigned[i] = 1
                detected_fps += fps_detected_by_assigned[i]
                    
        else:  # 1 < N ≤ Nmax_threshold
            # Count how many FPs each method can detect
            efx_count = np.sum(efx_detection[i])
            ncob_count = np.sum(ncob_detection[i])
            ecob_count = np.sum(ecob_detection[i])
            
            # Find the method that detects the most
            if efx_count >= ncob_count and efx_count >= ecob_count:
                assigned_metrics[i] = "EFX"
                count_efx_multi += 1
                fps_detected_by_assigned[i] = efx_count
                detected_fps += fps_detected_by_assigned[i]
            elif ncob_count >= ecob_count:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                count_ncob_n2plus += 1  # Count NCOB for N≥2 case
                fps_detected_by_assigned[i] = ncob_count
                detected_fps += fps_detected_by_assigned[i]
            else:
                assigned_metrics[i] = "ECOB"
                count_ecob += 1
                fps_detected_by_assigned[i] = ecob_count
                detected_fps += fps_detected_by_assigned[i]
    
    # Calculate effective efficiency
    effective_efficiency = (detected_fps / total_fps) * 100
    
    # Metric distribution percentages
    metric_counts = np.zeros(6)  # [EFX_zero, EFX_multi, SFX, NCOB, ECOB]
    metric_counts[0] = (count_efx_zero / n_targets) * 100
    metric_counts[1] = (count_efx_multi / n_targets) * 100  
    metric_counts[2] = (count_sfx / n_targets) * 100
    metric_counts[3] = (count_ncob_n1 / n_targets) * 100
    metric_counts[4] = (count_ncob_n2plus / n_targets) * 100
    metric_counts[5] = (count_ecob / n_targets) * 100
    
    # Resource usage percentages (for on-board constraints)
    # PLATO P5 targets per N-CAM: 80,000
    total_p5_targets_per_cam = 80000
    
    # Resource limits as absolute numbers and percentages
    efx_sfx_limit = 45000  # Maximum for flux measurements, either extended or secondary (no on-board distinction)
    efx_sfx_limit_pct = (efx_sfx_limit / total_p5_targets_per_cam) * 100  # 56.25%
    
    ncob_limit = 7400  # Maximum for nominal COB
    ncob_limit_pct = (ncob_limit / total_p5_targets_per_cam) * 100  # 9.25%
    
    ecob_limit = 7400  # Maximum for extended COB
    ecob_limit_pct = (ecob_limit / total_p5_targets_per_cam) * 100  # 9.25%
    
    # Calculate assignment percentages
    efx_total = count_efx_zero + count_efx_multi
    efx_sfx_used_pct = ((efx_total + count_sfx) / n_targets) * 100
    ncob_used_pct = (count_ncob / n_targets) * 100
    ecob_used_pct = (count_ecob / n_targets) * 100
    
    # Check if usage percentages would exceed limits
    within_limits = (efx_sfx_used_pct <= efx_sfx_limit_pct) and \
                    (ncob_used_pct <= ncob_limit_pct) and \
                    (ecob_used_pct <= ecob_limit_pct)
    
    # Distribution of targets by N value
    n_counts = np.zeros(3)  # [N=0, N=1, N≥2]
    n_counts[0] = np.sum(n_bad == 0) / n_targets * 100
    n_counts[1] = np.sum(n_bad == 1) / n_targets * 100
    n_counts[2] = np.sum(n_bad >= 2) / n_targets * 100
    
    # Calculate detection efficiency by detection method
    method_efficiencies = np.zeros(5)  # [EFX, SFX, NCOB, ECOB, Overall]
    
    if np.sum(assigned_metrics == "EFX") > 0:
        method_efficiencies[0] = np.sum(fps_detected_by_assigned[assigned_metrics == "EFX"]) / np.sum(total_fps_per_target[assigned_metrics == "EFX"]) * 100
    
    if np.sum(assigned_metrics == "SFX") > 0:
        method_efficiencies[1] = np.sum(fps_detected_by_assigned[assigned_metrics == "SFX"]) / np.sum(total_fps_per_target[assigned_metrics == "SFX"]) * 100
    
    if np.sum(assigned_metrics == "NCOB") > 0:
        method_efficiencies[2] = np.sum(fps_detected_by_assigned[assigned_metrics == "NCOB"]) / np.sum(total_fps_per_target[assigned_metrics == "NCOB"]) * 100
    
    if np.sum(assigned_metrics == "ECOB") > 0:
        method_efficiencies[3] = np.sum(fps_detected_by_assigned[assigned_metrics == "ECOB"]) / np.sum(total_fps_per_target[assigned_metrics == "ECOB"]) * 100
    
    method_efficiencies[4] = effective_efficiency  # Overall efficiency
    
    # Calculate NCOB breakdown percentages
    ncob_n1_pct = (count_ncob_n1 / n_targets) * 100
    ncob_n2plus_pct = (count_ncob_n2plus / n_targets) * 100

    # Calculate EFX breakdown percentages
    efx_zero_pct = (count_efx_zero / n_targets) * 100
    efx_multi_pct = (count_efx_multi / n_targets) * 100
    
    # Calculate what percentage of N=1 targets get NCOB
    n1_targets = np.sum(n_bad == 1)
    if n1_targets > 0:
        ncob_percentage_of_n1 = (count_ncob_n1 / n1_targets) * 100
    else:
        ncob_percentage_of_n1 = 0
    
    # Print results
    print(f"Effective efficiency: {effective_efficiency:.2f}%")
    print("\nMetric distribution:")
    print(f"  EFX (N=0): {metric_counts[0]:.2f}%")
    print(f"  EFX (N>=2): {metric_counts[1]:.2f}%")
    print(f"  EFX: {metric_counts[0] + metric_counts[1]:.2f}%")
    print(f"  SFX: {metric_counts[2]:.2f}%")
    print(f"  NCOB (N=1): {metric_counts[3]:.2f}%")
    print(f"  NCOB (N>=2): {metric_counts[4]:.2f}%")
    print(f"  NCOB: {metric_counts[3] + metric_counts[4]:.2f}%")
    print(f"  ECOB: {metric_counts[5]:.2f}%")
    
    print("\nTarget distribution by FP count:")
    print(f"  N=0: {n_counts[0]:.2f}%")
    print(f"  N=1: {n_counts[1]:.2f}%")
    print(f"  N≥2: {n_counts[2]:.2f}%")
    
    print("\nDetailed NCOB and EFX distribution:")
    print(f"  NCOB (N=1): {count_ncob_n1} targets ({ncob_n1_pct:.2f}% of all targets)")
    print(f"  NCOB (N≥2): {count_ncob_n2plus} targets ({ncob_n2plus_pct:.2f}% of all targets)")
    print(f"  NCOB (Total): {count_ncob} targets ({(count_ncob/n_targets)*100:.2f}% of all targets)")
    print(f"  {ncob_percentage_of_n1:.2f}% of N=1 targets are assigned NCOB")
    print(f"  EFX (N = 0): {count_efx_zero} targets ({efx_zero_pct:.2f}% of all targets)")
    print(f"  EFX (N >= 2): {count_efx_multi} targets ({efx_multi_pct:.2f}% of all targets)")
    
    print("\nResource usage (as % of all targets):")
    print(f"  EFX+SFX: {efx_sfx_used_pct:.2f}% (limit: {efx_sfx_limit_pct:.2f}%)")
    print(f"  NCOB: {ncob_used_pct:.2f}% (limit: {ncob_limit_pct:.2f}%)")
    print(f"  ECOB: {ecob_used_pct:.2f}% (limit: {ecob_limit_pct:.2f}%)")
    print(f"  Within limits: {within_limits}")
    
    print("\nMethod efficiency for assigned targets:")
    print(f"  EFX: {method_efficiencies[0]:.2f}%")
    print(f"  SFX: {method_efficiencies[1]:.2f}%")
    print(f"  NCOB: {method_efficiencies[2]:.2f}%")
    print(f"  ECOB: {method_efficiencies[3]:.2f}%")
    print(f"  Overall: {method_efficiencies[4]:.2f}%")
    
    # Prepare extended return tuple including detailed NCOB statistics
    return (effective_efficiency, metric_counts, n_counts, method_efficiencies, total_fps, detected_fps, n_targets,
            count_ncob_n1, count_ncob_n2plus, ncob_n1_pct, ncob_n2plus_pct, ncob_percentage_of_n1,
            efx_sfx_used_pct, ncob_used_pct, ecob_used_pct, within_limits)

# Example call

results = calculate_effective_efficiency(
    data, data_sec, data_ext,
    eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
    eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
    fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
    flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
    delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first, 
    depth_sig_scaling, 100
)