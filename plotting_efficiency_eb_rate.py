import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from mpl_toolkits.axes_grid1.inset_locator import inset_axes #type: ignore
from matplotlib.ticker import FuncFormatter  # type: ignore

from imagette import ran_unique_int

# ============================================================================
# Configuration 
# ============================================================================
# Data directories - point to EB rate results
dataDIR = '/home/fercho//double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/EBs_rate/1000_targets_per_magnitude_bin/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/'
DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/EBs_rate/Distribution_transit_depth_and_durations/'

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
gamma_factor_significance = 1
td_ref = 6.72*0.46**2
dback_ref = 132000
ntr = 3

# ============================================================================
# Helper Functions for EB Occurrence Rate
# ============================================================================
def get_actual_eb_count(sprk_array):
    """
    Determine the actual number of EBs for each target based on non-zero SPR values.
    
    Parameters:
    -----------
    sprk_array : numpy array of shape (n_targets, 10)
        Array containing SPR values for up to 10 contaminants per target
        
    Returns:
    --------
    numpy array of shape (n_targets,)
        Actual number of EBs per target
    """
    return np.sum(sprk_array > 0, axis=1)

def create_valid_eb_mask(sprk_array):
    """
    Create a boolean mask indicating valid EB positions.
    
    Parameters:
    -----------
    sprk_array : numpy array of shape (n_targets, 10)
        Array containing SPR values for up to 10 contaminants per target
        
    Returns:
    --------
    numpy array of shape (n_targets, 10)
        Boolean mask where True indicates a valid EB
    """
    return sprk_array > 0

def print_eb_statistics(actual_eb_counts, mag):
    """Print statistics about EB distribution."""
    print("\n" + "="*60)
    print("EB OCCURRENCE RATE STATISTICS")
    print("="*60)
    print(f"Total targets analyzed: {len(actual_eb_counts)}")
    print(f"Maximum EBs in any target: {int(np.max(actual_eb_counts))}")
    print(f"Average EBs per target: {np.mean(actual_eb_counts):.2f}")
    print(f"Median EBs per target: {np.median(actual_eb_counts):.1f}")
    print(f"\nEB Distribution:")
    print(f"  Targets with 0 EBs: {np.sum(actual_eb_counts == 0)} ({np.sum(actual_eb_counts == 0)/len(actual_eb_counts)*100:.1f}%)")
    print(f"  Targets with 1 EB: {np.sum(actual_eb_counts == 1)} ({np.sum(actual_eb_counts == 1)/len(actual_eb_counts)*100:.1f}%)")
    print(f"  Targets with 2 EBs: {np.sum(actual_eb_counts == 2)} ({np.sum(actual_eb_counts == 2)/len(actual_eb_counts)*100:.1f}%)")
    print(f"  Targets with 3+ EBs: {np.sum(actual_eb_counts >= 3)} ({np.sum(actual_eb_counts >= 3)/len(actual_eb_counts)*100:.1f}%)")
    print(f"  Targets with 5+ EBs: {np.sum(actual_eb_counts >= 5)} ({np.sum(actual_eb_counts >= 5)/len(actual_eb_counts)*100:.1f}%)")
    print(f"  Targets with 10+ EBs: {np.sum(actual_eb_counts >= 10)} ({np.sum(actual_eb_counts == 10)/len(actual_eb_counts)*100:.1f}%)")
    
    # Statistics by magnitude bin
    print(f"\nEB Distribution by Magnitude:")
    for i in range(nP):
        Pi = Pmin + i * binsize
        m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
        avg_ebs = np.mean(actual_eb_counts[m])
        print(f"  P = {Pi:.1f}: Average {avg_ebs:.2f} EBs per target")
    print("="*60 + "\n")

# ============================================================================
# Load Data
# ============================================================================
print("Loading data files...")
data_catalogue = np.load(cataDIR + 'SFP_DR3_20230101.npy')
magnitude_all_stars = data_catalogue[:, 2]
data = np.load(dataDIR + 'targets_P5.npy')
data_sec = np.load(dataDIR + 'targets_P5_secondary.npy')
data_ext = np.load(dataDIR + 'targets_P5_extended.npy')
data_bray = np.load(dataDIR + 'targets_P5_bray.npy')

# ============================================================================
# Extract Data Arrays
# ============================================================================
n = data.shape[0]
mag = data[:, 1]
n_bad = data[:, 8]
mask_p5 = (mag >= 10) & (mag <= 13)

# Get all the metrics
eta_t = data[:, 12]
eta_t_6_cameras = data[:, 45]
delta_obs_t = data[:, 13]
eta_c = data_sec[:, 6]
eta_c_6_cameras = data_sec[:, 11]
delta_obs_c = data_sec[:, 7]
nsr1h = data[:, 7]
nsr1h_6_cameras = data[:, 148]
nsr1h_sec = data_sec[:, 4]
nsr1h_sec_6_cameras = data_sec[:, 18]

# Get SPR values for first 10 contaminants
SPRK10_first = data[:, 17:27]
td_from_dap = data[:, 126:136]
dback_from_dap = data[:, 136:146]
eta_cob_ext_10first_6_cameras = data_ext[:, 75:85]
sigma_cob_ext_10first_6_cameras = data_ext[:, 85:95]
delta_cob_ext_10first_6_cameras = data_ext[:, 95:105]
gamma_cob_ext_10first_24_cameras = data_ext[:, 105:115]
mag_target_10first_contaminants = data[:, 169:179]
dist_target_to_10first_contaminants = data[:, 179:189]
eta_cob_nom_10first_6_cameras = data[:, 76:86]
sigma_cob_10first_6_cameras = data[:, 86:96]


# Get actual EB counts and create valid EB mask
actual_eb_counts = get_actual_eb_count(SPRK10_first)
valid_eb_mask = create_valid_eb_mask(SPRK10_first)

# Print EB statistics
print_eb_statistics(actual_eb_counts, mag)

# ============================================================================
# Compute Eta Arrays for Actual EBs Only
# ============================================================================
print("Computing eta arrays for actual EBs...")

# Initialize arrays
eta_nom_bt_24_cameras = np.zeros((n, 10))
eta_nom_bt_6_cameras = np.zeros((n, 10))
eta_ext_bt_24_cameras = np.zeros((n, 10))
eta_ext_bt_6_cameras = np.zeros((n, 10))
delta_obs = np.zeros((n, 10))
delta_obs_ext = np.zeros((n, 10))

# Compute values only for actual EBs
for i in range(n):
    n_ebs = actual_eb_counts[i]
    if n_ebs > 0:
        # Get values for actual EBs only
        dback = dback_from_dap[i, :n_ebs]
        td = td_from_dap[i, :n_ebs]
        
        # Compute eta values
        eta_nom_bt_24_cameras[i, :n_ebs] = (gamma_factor_significance * dback * 
                                            data[i, 17:17+n_ebs] * np.sqrt(td * ntr) / 
                                            (data[i, 7] * (1 - data[i, 11])))
        eta_nom_bt_6_cameras[i, :n_ebs] = (gamma_factor_significance * dback * 
                                           data[i, 17:17+n_ebs] * np.sqrt(td * ntr) / 
                                           (data[i, 148] * (1 - data[i, 11])))
        eta_ext_bt_24_cameras[i, :n_ebs] = (gamma_factor_significance * dback * 
                                            data_ext[i, 14:14+n_ebs] * np.sqrt(td * ntr) / 
                                            (data_ext[i, 4] * (1 - data_ext[i, 13])))
        eta_ext_bt_6_cameras[i, :n_ebs] = (gamma_factor_significance * dback * 
                                           data_ext[i, 14:14+n_ebs] * np.sqrt(td * ntr) / 
                                           (data_ext[i, 44] * (1 - data_ext[i, 13])))
        
        delta_obs[i, :n_ebs] = dback * SPRK10_first[i, :n_ebs]
        delta_obs_ext[i, :n_ebs] = dback * data_ext[i, 14:14+n_ebs]

# Save eta arrays
np.save(dataDIR + 'eta_bt_24_cameras_eb_rate.npy', eta_nom_bt_24_cameras)
np.save(dataDIR + 'eta_nom_bt_6_cameras_eb_rate.npy', eta_nom_bt_6_cameras)
np.save(dataDIR + 'eta_ext_bt_24_cameras_eb_rate.npy', eta_ext_bt_24_cameras)
np.save(dataDIR + 'eta_ext_bt_6_cameras_eb_rate.npy', eta_ext_bt_6_cameras)

# ============================================================================
# Prepare for Efficiency Calculations
# ============================================================================
# Get other needed arrays
eta_cob = data[:, 15]
eta_cob_6_cameras = data[:, 41]
eta_cob_sec_24_cameras = data_sec[:, 9]
eta_cob_sec_6_cameras = data_sec[:, 12]
eta_cob_nom_10first_24_cameras = data[:, 46:56]
eta_cob_ext_10first_24_cameras = data_ext[:, 45:55]

# COB metrics
delta_cob = data[:, 14]
delta_cob_6_cameras = data[:, 42]
delta_cob_sec = data_sec[:, 8]
delta_cob_sec_6_cameras = data_sec[:, 13]
sigma_cob_sec_24_cameras = data_sec[:, 10]
sigma_cob_sec_6_cameras = data_sec[:, 14]

# Signal depth calculations
sig_depth_secondary_mask_24_cameras = nsr1h_sec*(1 - data_sec[:, 5])/np.sqrt(td_ref*ntr) # type: ignore
sig_depth_secondary_mask_6_cameras = nsr1h_sec_6_cameras*(1 - data_sec[:, 5])/np.sqrt(td_ref*ntr) # type: ignore
sig_depth_extended_mask_24_cameras = data_ext[:, 4]*(1 - data_ext[:, 13])/np.sqrt(td_ref*ntr)
sig_depth_extended_mask_6_cameras = data_ext[:, 44]*(1 - data_ext[:, 13])/np.sqrt(td_ref*ntr)
sig_depth_nominal_mask_24_cameras = data[:, 7]*(1 - data[:, 11])/np.sqrt(td_ref*ntr)
sig_depth_nominal_mask_6_cameras = data[:, 148]*(1 - data[:, 11])/np.sqrt(td_ref*ntr)

# Reshape for 10 contaminants
sig_depth_24_cameras_10first = np.repeat(sig_depth_nominal_mask_24_cameras[:, np.newaxis], 10, axis=1)
sig_depth_24_cameras_10first = np.sqrt(sig_depth_24_cameras_10first**2 + 
                                       np.repeat(sig_depth_extended_mask_24_cameras[:, np.newaxis], 10, axis=1)**2)

# Detection conditions (for single contaminant)
fp_single_contaminant_24_cameras = (eta_t > flux_thresh_nom_mask)
fp_single_contaminant_6_cameras = (eta_t_6_cameras > flux_thresh_nom_mask)
secondary_mask_conditions_24_cameras = ((eta_c > flux_thresh_sec_mask) & 
                                       (delta_obs_c > delta_obs_t + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_24_cameras**2 + sig_depth_nominal_mask_24_cameras**2)) & 
                                       fp_single_contaminant_24_cameras)
secondary_mask_conditions_6_cameras = ((eta_c_6_cameras > flux_thresh_sec_mask) & 
                                      (delta_obs_c > delta_obs_t + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_6_cameras**2 + sig_depth_nominal_mask_6_cameras**2)) & 
                                      fp_single_contaminant_6_cameras)
secondary_mask_conditions_cob_24_cameras = ((eta_cob_sec_24_cameras > cob_thresh) & 
                                           (delta_cob_sec > 10*sigma_cob_sec_24_cameras) & 
                                           fp_single_contaminant_24_cameras)
secondary_mask_conditions_cob_6_cameras = ((eta_cob_sec_6_cameras > cob_thresh) & 
                                          (delta_cob_sec_6_cameras > 10*sigma_cob_sec_6_cameras) & 
                                          fp_single_contaminant_6_cameras)

# ============================================================================
# Plot 1: Flux Efficiency vs Magnitude
# ============================================================================
print("\nCreating flux efficiency plot...")

plt.figure(figsize=(10, 8))

# Store values for plotting
Pi_values = []
eff_sec_values = []
eff_sec_6_cameras_values = []
eff_ext_overall_24_cameras_values = []
eff_ext_overall_6_cameras_values = []

for i in range(nP):
    Pi = Pmin + i * binsize
    Pi_values.append(Pi)
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    
    # Count FPs with valid EBs
    fp_ext_overall_24_cameras = ((eta_nom_bt_24_cameras > flux_thresh_nom_mask) & valid_eb_mask)[m, :].sum()
    fp_ext_overall_6_cameras = ((eta_nom_bt_6_cameras > flux_thresh_nom_mask) & valid_eb_mask)[m, :].sum()
    
    # Calculate efficiencies
    if fp_ext_overall_24_cameras > 0:
        eff_ext_overall_24_cameras = (((eta_ext_bt_24_cameras > flux_thresh_ext_mask) & 
                                      (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & 
                                      (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                                      valid_eb_mask)[m, :].sum() / fp_ext_overall_24_cameras * 100)
    else:
        eff_ext_overall_24_cameras = 0
        
    if fp_ext_overall_6_cameras > 0:
        eff_ext_overall_6_cameras = (((eta_ext_bt_6_cameras > flux_thresh_ext_mask) & 
                                     (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & 
                                     (eta_nom_bt_6_cameras > flux_thresh_nom_mask) & 
                                     valid_eb_mask)[m, :].sum() / fp_ext_overall_6_cameras * 100)
    else:
        eff_ext_overall_6_cameras = 0
    
    eff_sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100 if fp_single_contaminant_24_cameras[m].sum() > 0 else 0
    eff_sec_6_cameras = (secondary_mask_conditions_6_cameras[m].sum() / fp_single_contaminant_6_cameras[m].sum()) * 100 if fp_single_contaminant_6_cameras[m].sum() > 0 else 0
    
    # Store values
    eff_sec_values.append(eff_sec)
    eff_sec_6_cameras_values.append(eff_sec_6_cameras)
    eff_ext_overall_24_cameras_values.append(eff_ext_overall_24_cameras)
    eff_ext_overall_6_cameras_values.append(eff_ext_overall_6_cameras)
    
    # Calculate errors
    n_targets_bin = m.sum()
    error_sec = np.sqrt(eff_sec * (100 - eff_sec) / n_targets_bin) if n_targets_bin > 0 else 0
    error_sec_6_cameras = np.sqrt(eff_sec_6_cameras * (100 - eff_sec_6_cameras) / n_targets_bin) if n_targets_bin > 0 else 0
    error_ext = np.sqrt(eff_ext_overall_24_cameras * (100 - eff_ext_overall_24_cameras) / n_targets_bin) if n_targets_bin > 0 else 0
    error_ext_6_cameras = np.sqrt(eff_ext_overall_6_cameras * (100 - eff_ext_overall_6_cameras) / n_targets_bin) if n_targets_bin > 0 else 0
    
    # Plot with error bars
    plt.errorbar(Pi, eff_sec, yerr=error_sec, fmt='o', color='purple', ecolor='purple', 
                capsize=5, label='Sec. Mask (24 cameras)' if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_sec_6_cameras, yerr=error_sec_6_cameras, fmt='o', color='green', 
                ecolor='green', capsize=5, label='Sec. Mask (6 cameras)' if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_24_cameras, yerr=error_ext, fmt='s', color='blue', 
                ecolor='blue', capsize=5, label='Ext. Mask (24 cameras)' if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_6_cameras, yerr=error_ext_6_cameras, fmt='s', color='red', 
                ecolor='red', capsize=5, label='Ext. Mask (6 cameras)' if i == 0 else "", markersize=4)

# Connect points with lines
plt.plot(Pi_values, eff_sec_values, color='purple', linestyle='-', linewidth=1)
plt.plot(Pi_values, eff_sec_6_cameras_values, color='green', linestyle='-', linewidth=1)
plt.plot(Pi_values, eff_ext_overall_24_cameras_values, color='blue', linestyle='-', linewidth=1)
plt.plot(Pi_values, eff_ext_overall_6_cameras_values, color='red', linestyle='-', linewidth=1)

# Add shaded regions and lines
plt.fill_between([9, 11.7], [47, 47], [100, 100], color='aqua', alpha=0.1) # type: ignore
plt.fill_between([11, 13.4], [47, 47], [100, 100], color='plum', alpha=0.1) # type: ignore
plt.vlines(11.7, ymin=47, ymax=100, linestyles='dashed', colors='green') # type: ignore
plt.vlines(11, ymin=47, ymax=100, linestyles='dashdot', colors='red') # type: ignore

# Labels and formatting
plt.ylim(47, 100)
plt.xlim(9.9, 13.1)
plt.text(10, 78.1, 'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold')
plt.text(11.2, 75, 'On-board light curve processing region', color='red', weight='bold')
plt.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', borderaxespad=0., ncol=2)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)
plt.title('DAP Flux Efficiency (EB Occurrence Rate = 1%)', fontsize=fsize+2)
plt.tight_layout(rect=[0, 0, 1, 0.88])
plt.savefig(DIRout + "DAP_Flux_efficiency_EB_rate.pdf", format='pdf', bbox_inches='tight')
plt.show()

# ============================================================================
# Plot 2: COB Efficiency vs Magnitude - Fixed Implementation
# ============================================================================
print("\nCreating COB efficiency plot...")

plt.figure(figsize=(10, 8))

# Initialize lists to store data points
Pi_values_cob = []
eff_ext_cob_overall_values = []
eff_ext_cob_overall_6_cameras_values = []
eff_cob_values = []
eff_cob_6_cameras_values = []
eff_cob_sec_values = []
eff_cob_sec_6_cameras_values = []

# Arrays to store values for lines (only non-zero points)
plot_data = {
    'ext_24': {'pi': [], 'eff': [], 'err': []},
    'ext_6': {'pi': [], 'eff': [], 'err': []},
    'nom_24': {'pi': [], 'eff': [], 'err': []},
    'nom_6': {'pi': [], 'eff': [], 'err': []},
    'sec_24': {'pi': [], 'eff': [], 'err': []},
    'sec_6': {'pi': [], 'eff': [], 'err': []}
}

# Process each magnitude bin
for i in range(nP):
    Pi = Pmin + i * binsize
    Pi_values_cob.append(Pi)
    
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    
    # Count targets with actual EBs for this magnitude bin
    targets_with_ebs = (actual_eb_counts[m] > 0).sum()
    n_targets_bin = m.sum()
    
    print(f"P={Pi:.1f}: {n_targets_bin} targets, {targets_with_ebs} with EBs")
    
    # Calculate denominators
    s_24_cameras = ((eta_nom_bt_24_cameras > flux_thresh_nom_mask) & valid_eb_mask)[m,:].sum()
    s_6_cameras = ((eta_nom_bt_6_cameras > flux_thresh_nom_mask) & valid_eb_mask)[m,:].sum()
    fp_single = fp_single_contaminant_24_cameras[m].sum()
    fp_single_6 = fp_single_contaminant_6_cameras[m].sum()
    
    # Calculate COB efficiencies with zero checks
    # Extended COB 24 cameras
    if s_24_cameras > 0:
        eff_ext_cob_overall = (((eta_cob_ext_10first_24_cameras > cob_thresh) & 
                               (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                               valid_eb_mask)[m,:].sum() / s_24_cameras * 100.)
        error_ext_cob = np.sqrt(eff_ext_cob_overall * (100 - eff_ext_cob_overall) / max(targets_with_ebs, 1))
        if eff_ext_cob_overall > 0:
            plot_data['ext_24']['pi'].append(Pi)
            plot_data['ext_24']['eff'].append(eff_ext_cob_overall)
            plot_data['ext_24']['err'].append(error_ext_cob)
    else:
        eff_ext_cob_overall = 0
        error_ext_cob = 0
    
    # Nominal COB 24 cameras    
    if s_24_cameras > 0:
        eff_cob = (((eta_cob_nom_10first_24_cameras > cob_thresh) & 
                   (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                   valid_eb_mask)[m,:].sum() / s_24_cameras * 100.)
        error_cob = np.sqrt(eff_cob * (100 - eff_cob) / max(targets_with_ebs, 1))
        if eff_cob > 0:
            plot_data['nom_24']['pi'].append(Pi)
            plot_data['nom_24']['eff'].append(eff_cob)
            plot_data['nom_24']['err'].append(error_cob)
    else:
        eff_cob = 0
        error_cob = 0
    
    # Extended COB 6 cameras
    if s_6_cameras > 0:
        eff_ext_cob_overall_6_cameras = (((eta_cob_ext_10first_6_cameras > cob_thresh) & 
                                         (eta_nom_bt_6_cameras > flux_thresh_nom_mask) & 
                                         valid_eb_mask)[m,:].sum() / s_6_cameras * 100.)
        error_ext_cob_6_cameras = np.sqrt(eff_ext_cob_overall_6_cameras * (100 - eff_ext_cob_overall_6_cameras) / max(targets_with_ebs, 1))
        if eff_ext_cob_overall_6_cameras > 0:
            plot_data['ext_6']['pi'].append(Pi)
            plot_data['ext_6']['eff'].append(eff_ext_cob_overall_6_cameras)
            plot_data['ext_6']['err'].append(error_ext_cob_6_cameras)
    else:
        eff_ext_cob_overall_6_cameras = 0
        error_ext_cob_6_cameras = 0
    
    # Nominal COB 6 cameras
    if s_6_cameras > 0:
        eff_cob_6_cameras = (((eta_cob_nom_10first_6_cameras > cob_thresh) & 
                            (eta_nom_bt_6_cameras > flux_thresh_nom_mask) & 
                            valid_eb_mask)[m,:].sum() / s_6_cameras * 100.)
        error_cob_6_cameras = np.sqrt(eff_cob_6_cameras * (100 - eff_cob_6_cameras) / max(targets_with_ebs, 1))
        if eff_cob_6_cameras > 0:
            plot_data['nom_6']['pi'].append(Pi)
            plot_data['nom_6']['eff'].append(eff_cob_6_cameras)
            plot_data['nom_6']['err'].append(error_cob_6_cameras)
    else:
        eff_cob_6_cameras = 0
        error_cob_6_cameras = 0
    
    # Secondary COB 24 cameras
    if fp_single > 0:
        eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_single) * 100
        error_cob_sec = np.sqrt(eff_cob_sec * (100 - eff_cob_sec) / max(targets_with_ebs, 1))
        if eff_cob_sec > 0:
            plot_data['sec_24']['pi'].append(Pi)
            plot_data['sec_24']['eff'].append(eff_cob_sec)
            plot_data['sec_24']['err'].append(error_cob_sec)
    else:
        eff_cob_sec = 0
        error_cob_sec = 0
    
    # Secondary COB 6 cameras    
    if fp_single_6 > 0:
        eff_cob_sec_6_cameras = (secondary_mask_conditions_cob_6_cameras[m].sum() / fp_single_6) * 100
        error_cob_sec_6_cameras = np.sqrt(eff_cob_sec_6_cameras * (100 - eff_cob_sec_6_cameras) / max(targets_with_ebs, 1))
        if eff_cob_sec_6_cameras > 0:
            plot_data['sec_6']['pi'].append(Pi)
            plot_data['sec_6']['eff'].append(eff_cob_sec_6_cameras)
            plot_data['sec_6']['err'].append(error_cob_sec_6_cameras)
    else:
        eff_cob_sec_6_cameras = 0
        error_cob_sec_6_cameras = 0

    # Store all values (including zeros) for completeness
    eff_ext_cob_overall_values.append(eff_ext_cob_overall)
    eff_ext_cob_overall_6_cameras_values.append(eff_ext_cob_overall_6_cameras)
    eff_cob_values.append(eff_cob)
    eff_cob_6_cameras_values.append(eff_cob_6_cameras)
    eff_cob_sec_values.append(eff_cob_sec)
    eff_cob_sec_6_cameras_values.append(eff_cob_sec_6_cameras)

# Now plot the data - only plot if there's data to show
# Extended COB 24 cameras
if len(plot_data['ext_24']['pi']) > 0:
    plt.errorbar(plot_data['ext_24']['pi'], plot_data['ext_24']['eff'], 
                yerr=plot_data['ext_24']['err'], fmt='s', color='blue', 
                ecolor='blue', capsize=5, markersize=4, label='Ext. Mask (24 cameras)')
    if len(plot_data['ext_24']['pi']) > 1:
        plt.plot(plot_data['ext_24']['pi'], plot_data['ext_24']['eff'], 
                color='blue', linestyle='-', linewidth=1, alpha=0.5)

# Extended COB 6 cameras
if len(plot_data['ext_6']['pi']) > 0:
    plt.errorbar(plot_data['ext_6']['pi'], plot_data['ext_6']['eff'], 
                yerr=plot_data['ext_6']['err'], fmt='s', color='red', 
                ecolor='red', capsize=5, markersize=4, label='Ext. Mask (6 cameras)')
    if len(plot_data['ext_6']['pi']) > 1:
        plt.plot(plot_data['ext_6']['pi'], plot_data['ext_6']['eff'], 
                color='red', linestyle='-', linewidth=1, alpha=0.5)

# Nominal COB 24 cameras
if len(plot_data['nom_24']['pi']) > 0:
    plt.errorbar(plot_data['nom_24']['pi'], plot_data['nom_24']['eff'], 
                yerr=plot_data['nom_24']['err'], fmt='*', color='orange', 
                ecolor='orange', capsize=5, markersize=8, label='Nom. Mask (24 cameras)')
    if len(plot_data['nom_24']['pi']) > 1:
        plt.plot(plot_data['nom_24']['pi'], plot_data['nom_24']['eff'], 
                color='orange', linestyle='-', linewidth=1, alpha=0.5)

# Nominal COB 6 cameras
if len(plot_data['nom_6']['pi']) > 0:
    plt.errorbar(plot_data['nom_6']['pi'], plot_data['nom_6']['eff'], 
                yerr=plot_data['nom_6']['err'], fmt='*', color='olive', 
                ecolor='olive', capsize=5, markersize=8, label='Nom. Mask (6 cameras)')
    if len(plot_data['nom_6']['pi']) > 1:
        plt.plot(plot_data['nom_6']['pi'], plot_data['nom_6']['eff'], 
                color='olive', linestyle='-', linewidth=1, alpha=0.5)

# Secondary COB 24 cameras
if len(plot_data['sec_24']['pi']) > 0:
    plt.errorbar(plot_data['sec_24']['pi'], plot_data['sec_24']['eff'], 
                yerr=plot_data['sec_24']['err'], fmt='o', color='purple', 
                ecolor='purple', capsize=5, markersize=4, label='Sec. Mask (24 cameras)')
    if len(plot_data['sec_24']['pi']) > 1:
        plt.plot(plot_data['sec_24']['pi'], plot_data['sec_24']['eff'], 
                color='purple', linestyle='-', linewidth=1, alpha=0.5)

# Secondary COB 6 cameras
if len(plot_data['sec_6']['pi']) > 0:
    plt.errorbar(plot_data['sec_6']['pi'], plot_data['sec_6']['eff'], 
                yerr=plot_data['sec_6']['err'], fmt='o', color='green', 
                ecolor='green', capsize=5, markersize=4, label='Sec. Mask (6 cameras)')
    if len(plot_data['sec_6']['pi']) > 1:
        plt.plot(plot_data['sec_6']['pi'], plot_data['sec_6']['eff'], 
                color='green', linestyle='-', linewidth=1, alpha=0.5)

# Add shaded regions and lines
plt.fill_between([9, 11.7], [0, 0], [100, 100], color='aqua', alpha=0.1) # type: ignore
plt.fill_between([11, 13.4], [0, 0], [100, 100], color='plum', alpha=0.1) # type: ignore
plt.vlines(11.7, ymin=0, ymax=100, linestyles='dashed', colors='green') # type: ignore
plt.vlines(11, ymin=0, ymax=100, linestyles='dashdot', colors='red') # type: ignore

# Labels and formatting
plt.ylim(20, 100)
plt.xlim(9.9, 13.1)
plt.text(10., 55, 'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold', fontsize=12)
plt.text(11.1, 55, 'On-board light curve\nprocessing region', color='red', weight='bold', fontsize=12)

# Only show legend if there are labeled artists - place it outside the plot
if plt.gca().get_legend_handles_labels()[0]:
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), borderaxespad=0., fancybox=True, ncol=3, columnspacing=0.6)

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)
plt.title('DAP COB Efficiency (EB Occurrence Rate = 1%)', fontsize=fsize+2)
plt.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

# Save figure
plt.tight_layout(rect=[0, 0, 1, 0.88])  # Leave space for legend below
plt.savefig(DIRout + "DAP_CS_efficiency_EB_rate.pdf", format='pdf', bbox_inches='tight')
plt.show()

# Print summary statistics
print("\nCOB Efficiency Summary:")
print("-" * 50)
for key, label in [('ext_24', 'Extended 24-cam'), ('ext_6', 'Extended 6-cam'),
                   ('nom_24', 'Nominal 24-cam'), ('nom_6', 'Nominal 6-cam'),
                   ('sec_24', 'Secondary 24-cam'), ('sec_6', 'Secondary 6-cam')]:
    if len(plot_data[key]['eff']) > 0:
        avg_eff = np.mean(plot_data[key]['eff'])
        print(f"{label}: Average {avg_eff:.1f}% (from {len(plot_data[key]['eff'])} magnitude bins with data)")
    else:
        print(f"{label}: No data points")

# ============================================================================
# Weighted Efficiency Calculations
# ============================================================================
print("\nCalculating weighted efficiencies across all magnitude bins...")

# Load star catalogue for weighting
mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])

# Calculate weights based on star distribution
star_counts = np.zeros(nP, dtype=int)
total_stars = 0
for i in range(nP):
    Pi = Pmin + i * binsize
    mask_bin = (magnitude_all_stars >= Pi - binsize / 2.) & (magnitude_all_stars <= Pi + binsize / 2.)
    star_counts[i] = mask_bin.sum()
    total_stars += star_counts[i]

weights = star_counts / float(total_stars)

# Define detection conditions for all contaminants (with EB mask)
nfp = (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & valid_eb_mask
nfp_ext_mask = ((eta_ext_bt_24_cameras > flux_thresh_ext_mask) & 
                (delta_obs_ext > delta_obs + depth_sig_scaling*sig_depth_24_cameras_10first) & 
                (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                valid_eb_mask)
nfp_nom_cob = ((eta_cob_nom_10first_24_cameras > cob_thresh) & 
               (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
               valid_eb_mask)
nfp_ext_cob = ((eta_cob_ext_10first_24_cameras > cob_thresh) & 
               (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
               valid_eb_mask)
nfp_sec_mask = secondary_mask_conditions_24_cameras
nfp_sec_cob = secondary_mask_conditions_cob_24_cameras

# Initialize weighted sums
weighted_eff_ext_flux = 0
weighted_eff_sec_flux = 0
weighted_eff_nom_cob = 0
weighted_eff_ext_cob = 0
weighted_eff_sec_cob = 0
weighted_variance_eff_ext_flux = 0
weighted_variance_eff_sec_flux = 0
weighted_variance_eff_nom_cob = 0
weighted_variance_eff_ext_cob = 0
weighted_variance_eff_sec_cob = 0

# Calculate weighted efficiencies
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize / 2.) & (mag <= Pi + binsize / 2.)
    n_targets_bin = m.sum()
    
    # Count FPs and detections
    fp_count = nfp[m].sum()
    fp_single = fp_single_contaminant_24_cameras[m].sum()
    
    # Calculate efficiencies
    if fp_count > 0:
        eff_ext_flux = nfp_ext_mask[m].sum() / fp_count * 100.
        eff_nom_cob = nfp_nom_cob[m].sum() / fp_count * 100.
        eff_ext_cob = nfp_ext_cob[m].sum() / fp_count * 100.
    else:
        eff_ext_flux = 0
        eff_nom_cob = 0
        eff_ext_cob = 0
        
    if fp_single > 0:
        eff_sec_flux = nfp_sec_mask[m].sum() / fp_single * 100
        eff_sec_cob = nfp_sec_cob[m].sum() / fp_single * 100
    else:
        eff_sec_flux = 0
        eff_sec_cob = 0
    
    # Accumulate weighted sums
    weighted_eff_ext_flux += weights[i] * eff_ext_flux
    weighted_eff_sec_flux += weights[i] * eff_sec_flux
    weighted_eff_nom_cob += weights[i] * eff_nom_cob
    weighted_eff_ext_cob += weights[i] * eff_ext_cob
    weighted_eff_sec_cob += weights[i] * eff_sec_cob
    
    # Calculate and accumulate variances
    if n_targets_bin > 0:
        weighted_variance_eff_ext_flux += weights[i] * (eff_ext_flux * (100 - eff_ext_flux)) / n_targets_bin
        weighted_variance_eff_sec_flux += weights[i] * (eff_sec_flux * (100 - eff_sec_flux)) / n_targets_bin
        weighted_variance_eff_nom_cob += weights[i] * (eff_nom_cob * (100 - eff_nom_cob)) / n_targets_bin
        weighted_variance_eff_ext_cob += weights[i] * (eff_ext_cob * (100 - eff_ext_cob)) / n_targets_bin
        weighted_variance_eff_sec_cob += weights[i] * (eff_sec_cob * (100 - eff_sec_cob)) / n_targets_bin

# Compute errors
weighted_error_eff_ext_flux = np.sqrt(weighted_variance_eff_ext_flux)
weighted_error_eff_sec_flux = np.sqrt(weighted_variance_eff_sec_flux)
weighted_error_eff_nom_cob = np.sqrt(weighted_variance_eff_nom_cob)
weighted_error_eff_ext_cob = np.sqrt(weighted_variance_eff_ext_cob)
weighted_error_eff_sec_cob = np.sqrt(weighted_variance_eff_sec_cob)

# Print weighted efficiency results
print("\n" + "="*60)
print("WEIGHTED EFFICIENCY RESULTS (EB Occurrence Rate = 1%)")
print("="*60)
print(f"Extended flux efficiency: {weighted_eff_ext_flux:.2f}% ± {weighted_error_eff_ext_flux:.2f}%")
print(f"Secondary flux efficiency: {weighted_eff_sec_flux:.2f}% ± {weighted_error_eff_sec_flux:.2f}%")
print(f"Nominal COB efficiency: {weighted_eff_nom_cob:.2f}% ± {weighted_error_eff_nom_cob:.2f}%")
print(f"Extended COB efficiency: {weighted_eff_ext_cob:.2f}% ± {weighted_error_eff_ext_cob:.2f}%")
print(f"Secondary COB efficiency: {weighted_eff_sec_cob:.2f}% ± {weighted_error_eff_sec_cob:.2f}%")
print("="*60)

# ============================================================================
# Effective Efficiency Calculation (Adapted for EB Rate)
# ============================================================================
print("\nCalculating effective efficiency with optimal metric assignment...")

def calculate_effective_efficiency_eb_rate(data, data_sec, data_ext, 
                                          eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
                                          eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
                                          fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
                                          flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
                                          delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first,
                                          depth_sig_scaling, Nmax_threshold, valid_eb_mask):
    """
    Calculate the effective efficiency for EB occurrence rate scenario.
    Adapted to handle variable numbers of EBs per target.
    """
    mag = data[:, 1]
    n_bad = data[:, 8]
    n_targets = len(mag)
    
    # Arrays to store assigned metrics and detected FPs
    assigned_metrics = np.empty(n_targets, dtype='U4') # EFX, SFX, NCOB, ECOB, SCOB
    fps_detected_by_assigned = np.zeros(n_targets)
    total_fps_per_target = np.zeros(n_targets)
    
    # Counters
    count_efx_zero = 0 # Count of EFX assigned to targets with N = 0
    count_efx_multi = 0 # Count of EFX assigned to targets with N ≥ 2 
    count_sfx = 0 # Count of assigned SFX
    count_ncob = 0 # Counf of assigned NCOB
    count_ncob_single = 0 # Count of NCOB assigned to targets with N = 1
    count_ncob_multi = 0 # Count of NCOB assignted to targets with N ≥ 2
    count_ecob = 0 # Counf of assigned ECOB
    count_none = 0 # For cases with N > threshold
    
    total_fps = 0
    detected_fps = 0
    
    # Detection capabilities with EB mask
    efx_detection = ((eta_ext_bt_24_cameras > flux_thresh_ext_mask) & 
                     (delta_obs_ext > delta_obs + depth_sig_scaling * sig_depth_24_cameras_10first) & 
                     (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                     valid_eb_mask)
    
    sfx_detection = secondary_mask_conditions_24_cameras
    
    ncob_detection = ((eta_cob_nom_10first_24_cameras > cob_thresh) & 
                      (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & 
                      valid_eb_mask)
    
    ecob_detection = ((eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & valid_eb_mask)
    
    # Process each target
    for i in range(n_targets):
        N = int(n_bad[i])
        total_fps_per_target[i] = N
        total_fps += N
        
        if N > Nmax_threshold:
            assigned_metrics[i] = "NONE"
            count_none += 1
            fps_detected_by_assigned[i] = 0
            
        elif N == 0:
            assigned_metrics[i] = "EFX"
            count_efx_zero += 1
            fps_detected_by_assigned[i] = 0
            
        elif N == 1:
            if sfx_detection[i]:
                assigned_metrics[i] = "SFX"
                count_sfx += 1
                fps_detected_by_assigned[i] = 1
            else:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                count_ncob_single += 1
                fps_detected_by_assigned[i] = 1
            detected_fps += fps_detected_by_assigned[i]
                    
        else:  # 1 < N ≤ Nmax_threshold
            # Count detections for each method
            efx_count = np.sum(efx_detection[i])
            ncob_count = np.sum(ncob_detection[i])
            ecob_count = np.sum(ecob_detection[i])
            
            # Assign best method
            if efx_count >= ncob_count and efx_count >= ecob_count:
                assigned_metrics[i] = "EFX"
                count_efx_multi += 1
                fps_detected_by_assigned[i] = efx_count
            elif ncob_count >= ecob_count:
                assigned_metrics[i] = "NCOB"
                count_ncob += 1
                fps_detected_by_assigned[i] = ncob_count
            else:
                assigned_metrics[i] = "ECOB"
                count_ecob += 1
                fps_detected_by_assigned[i] = ecob_count
            detected_fps += fps_detected_by_assigned[i]
    
    # Calculate effective efficiency
    effective_efficiency = (detected_fps / total_fps) * 100 if total_fps > 0 else 0

    # Metric distribution percentages
    metric_counts = np.zeros(6) # [EFX_zero, EFX_multi, SFX, NCOB_one, NCOB_multi, ECOB]
    metric_counts[0] = (count_efx_zero / n_targets) * 100
    metric_counts[1] = (count_efx_multi / n_targets) * 100
    metric_counts[2] = (count_sfx / n_targets) * 100
    metric_counts[3] = (count_ncob_single / n_targets) * 100
    metric_counts[4] = (count_ncob_multi / n_targets) * 100
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
    ncob_single_pct = (count_ncob_single / n_targets) * 100
    ncob_multi_pct = (count_ncob_multi / n_targets) * 100

    # Calculate EFX breakdown percentages
    efx_zero_pct = (count_efx_zero / n_targets) * 100
    efx_multi_pct = (count_efx_multi / n_targets) * 100

    # Calculate what percentage of N=1 targets get NCOB
    n1_targets = np.sum(n_bad == 1)
    if n1_targets > 0:
        ncob_percentage_of_n1 = (count_ncob_single / n1_targets) * 100
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
    print(f"  NCOB (N=1): {count_ncob_single} targets ({ncob_single_pct:.2f}% of all targets)")
    print(f"  NCOB (N≥2): {count_ncob_multi} targets ({ncob_multi_pct:.2f}% of all targets)")
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

    # After the assignment loop, check:
    print("\nEFX assignment breakdown:")
    print(f"  Targets with N=0 assigned EFX: {np.sum((n_bad == 0) & (assigned_metrics == 'EFX'))}")
    print(f"  Targets with N=1 assigned EFX: {np.sum((n_bad == 1) & (assigned_metrics == 'EFX'))}")  
    print(f"  Targets with N≥2 assigned EFX: {np.sum((n_bad >= 2) & (assigned_metrics == 'EFX'))}")

    # Also check what N≥2 targets are getting:
    n2_plus_mask = n_bad >= 2
    if n2_plus_mask.sum() > 0:
        print(f"\nFor {n2_plus_mask.sum()} targets with N≥2:")
        for method in ['EFX', 'SFX', 'NCOB', 'ECOB']:
            count = np.sum((n2_plus_mask) & (assigned_metrics == method)) # type: ignore
            print(f"  {method}: {count}") # type: ignore
    
    return (effective_efficiency, assigned_metrics, metric_counts, n_counts, method_efficiencies, total_fps, detected_fps, n_targets, count_ncob_single, count_ncob_multi,
            ncob_single_pct, ncob_multi_pct, ncob_percentage_of_n1, efx_sfx_used_pct, ncob_used_pct, ecob_used_pct, within_limits)

# Call the function
results = calculate_effective_efficiency_eb_rate(
    data, data_sec, data_ext,
    eta_nom_bt_24_cameras, eta_ext_bt_24_cameras,
    eta_c, eta_cob_nom_10first_24_cameras, eta_cob_ext_10first_24_cameras,
    fp_single_contaminant_24_cameras, secondary_mask_conditions_24_cameras,
    flux_thresh_nom_mask, flux_thresh_ext_mask, flux_thresh_sec_mask, cob_thresh,
    delta_obs, delta_obs_ext, delta_obs_t, sig_depth_24_cameras_10first, 
    depth_sig_scaling, 100, valid_eb_mask
)

print("\n" + "="*60)
print("Processing complete!")
print(f"Results saved to: {DIRout}")
print("="*60)