import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.ticker import PercentFormatter # type: ignore

from imagette import ran_unique_int

#dataDIR = '/home/fercho/double-aperture-photometry/test_results/'
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr_Leopold_PSF/'
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


data_catalogue = np.load(cataDIR + 'SFP_DR3_20230101.npy') # star catalogue from GAIA
magnitude_all_stars = data_catalogue[:, 2]
data = np.load(dataDIR + 'targets_P5.npy')
#0: ID_t
#1: P_t
#2: n_c
#3: magnitude of the most sinificant contaminant
#4: distance from the target to the most significant contaminant
#5: w_t_key
#6: w_t_size 
#7: NSR1h
#8: n_bad
#9: SPR_crit
#10: sprk[ind_sprk]
#11: SPR_tot
#12: eta_t
#13: delta_obs_t
#14: abs_cob
#15: eta_cob
#16: sigma_1_24
#17-26: SPRK_10first
#27-36: eta_10first
#37: eta_true_positive_corrected_24_cameras
#38: eta_true_positive_corrected_6_cameras
#39: eta_true_positive_corrected_24_cameras_super-earth
#40: eta_true_positive_corrected_6_cameras_super-earth
#41: eta_cob_6_cameras
#42: abs_cob_6_cameras
#43: sigma_1_6_cameras
#44: SPR_crit_6_cameras
#45: eta_t_6_cameras
#46-55: eta_cob_10first
#56-65: sigma_cob_10first
#66-75: abs_cob_shift_10first
#76-85: eta_cob_10first_6_cameras
#86-95: sigma_cob_10first_6_cameras
#96-105: abs_cob_shift_10first_6_cameras
#106-115: gamma_nom_10first
#116-125: gamma_nom_10first_6_cameras
#126-135: td_10first
#136-145: dback_10first
#146: gamma_nom
#147: gamma_nom_6_cameras
#148: nsr1h_6cameras
#149-158: x_coordinate_contaminant_stars_10first
#159-168: y_coordinate_contaminant_stars_10first
#169-178: magnitude_10first_contaminants
#179-188: distance from target to 10 first contaminants
#189-198: IDs_from_the_10first_contaminants
#199-208: delta_x_from_target_to_10first_contaminants
#209-218: delta_y_from_target_to_10first_contaminants
#219: x_tar
#220: y_tar
 
data_sec = np.load(dataDIR + 'targets_P5_secondary.npy')
#0: ID_t
#1: P_t
#2: secondary_mask_key
#3: secondary_mask_size
#4: nsr_1h_24_cameras_secondary_mask
#5: spr_tot_secondary_mask
#6: eta_c
#7: delta_obs_secondary_mask
#8: abs_cob_c
#9: eta_cob_c
#10: sigma_1_24_c
#11: eta_c_6_cameras
#12: eta_cob_c_6_cameras
#13: abs_cob_c_6_cameras
#14: sigma_1_6_cameras_c
#15: gamma_cob_c_6_cameras
#16: SPR_tot_sec_6_cameras
#17: gamma_cob_c
#18: nsr_1h_6_cameras_secondary_mask
#19: ID_contaminant_highest_sprk

data_ext = np.load(dataDIR + 'targets_P5_extended.npy')
#0: ID_t
#1: P_t
#2: extended mask_key
#3: extended_mask_size
#4: NSR_ext_1h_24_cameras
#5: sprk_ext[ind_sprk]
#6: SPR_crit_ext
#7: eta_ext
#8: delta_obs_ext
#9: abs_cob_ext
#10: eta_cob_ext
#11: sigma_1_24_ext
#12: n_bad_ext
#13: SPR_tot_ext
#14-23: SPRK_ext_10first
#24-33: eta_ext_10first
#34-43: eta_ext_10first_6_cameras
#44: NSR_ext_1h_6_cameras
#45-54: eta_cob_ext_10first
#55-64: sigma_cob_ext_10first
#65-74: abs_cob_shift_ext_10first
#75-84: eta_cob_ext_10first_6_cameras
#85-94: sigma_cob_ext_10first_6_cameras
#95-104: abs_cob_shift_ext_10first_6_cameras
#105-114: gamma_ext_10first
#115-124: gamma_ext_10first_6_cameras
#125: gamma_ext
#126-135: SPRK_ext_2_pix_10first
#136-145: eta_ext_2_pix_10first
#146-155: eta_ext_2_pix_10first_6_cameras
#156: NSR_ext_2_pix_1h_6_cameras
#157-166: eta_cob_ext_2_pix_10first
#167-176: sigma_cob_ext_2_pix_10first
#177-186: abs_cob_shift_ext_2_pix_10first
#187-196: eta_cob_ext_2_pix_10first_6_cameras
#197-206: sigma_cob_ext_2_pix_10first_6_cameras
#207-216: abs_cob_shift_ext_2_pix_10first_6_cameras
#217-226: gamma_ext_2_pix_10first
#227-236: gamma_ext_2_pix_10first_6_cameras
#237: gamma_ext_2_pix
#238: NSR_ext_2_pix_1h_24_cameras
#239: extended_mask_key_2_pix
#240: extended_mask_size_2_pix
#241: delta_obs_ext_2_pix
#242: eta_ext_2_pix
#243: n_bad_ext_2_pix
#244: eta_cob_ext_2_pix
#245: sigma_1_24_ext_2_pix
#246: abs_cob_ext_2_pix
#247: SPR_crit_ext_2_pix
#248: SPR_to_ext_2_pix

data_bray = np.load(dataDIR + 'targets_P5_bray.npy')
#0: ID_t
#1: P_t
#2: n_c
#3: NSR_bray_1h
#4: n_bad_bray
#5: SPR_crit_bray
#6: sprk_bray[ind_sprk]
#7: SPR_tot_bray

mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])
ntr = 3        # number of transits in one hour
n = data.shape[0]
td = data[:, 126:136]
dback = data[:, 136:146]
seed = 123434434

plt.figure(0)
plt.plot(mag_value, star_count, 'o')
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel('Numb. of stars', fontsize=fsize)

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
spr_crit = data[:, 8]
nsr1h_ext = data_ext[:, 4]
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]
#eta_cob_ext_2_pix = data_ext[:, 244]
#sigma_1_24_ext_2_pix = data_ext[:, 245]
#delta_cob_ext_2_pix = data_ext[:, 246]

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
SPRK10_first = data_ext[:, 14:24]

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
# We get now the 10 first values for each param. of the nominal cob shift

# We get now the 10 first values for each param. of the extended cob shift
SPRK10_first_ext = data_ext[:,14:24]
eta_cob_ext_10first_24_cameras = data_ext[:, 45:55]
sigma_cob_ext_10first_24_cameras = data_ext[:, 55:65]
delta_cob_ext_10first_24_cameras = data_ext[:, 65:75]
eta_cob_ext_10first_6_cameras = data_ext[:, 75:85]
sigma_cob_ext_10first_6_cameras = data_ext[:, 85:95]
delta_cob_ext_10first_6_cameras = data_ext[:, 95:105]
gamma_cob_ext_10first_24_cameras = data_ext[:, 105:115]
#SPRK10_first_ext_2_pix = data_ext[:,126:136]
#eta_cob_ext_2_pix_10first_24_cameras = data_ext[:, 157:167]
#sigma_cob_ext_2_pix_10first_24_cameras = data_ext[:, 167:177]
#delta_cob_ext_2_pix_10first_24_cameras = data_ext[:, 177:187]
#eta_cob_ext_2_pix_10first_6_cameras = data_ext[:, 187:197]
#sigma_cob_ext_2_pix_10first_6_cameras = data_ext[:, 197:207]
#delta_cob_ext_2_pix_10first_6_cameras = data_ext[:, 207:217]
#gamma_cob_ext_2_pix_10first_24_cameras = data_ext[:, 217:227]
#gamma_cob_ext_2_pix_10first_24_cameras = data_ext[:, 227:237]
#gamma_ext_2_pix = data_ext[:, 237]



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
mag_2d = np.repeat(mag[:,np.newaxis], 10, axis=1)
eta_ext_2_pix_bt_24_cameras = np.zeros((n, 10))
eta_ext_2_pix_bt_6_cameras = np.zeros((n, 10))

for i in range(n):
    dback = np.ones(10)*dback_ref
    td = np.ones(10)*td_ref
    eta_nom_bt_24_cameras[i, :] = gamma_factor_significance*dback*data[i, 17:27]*np.sqrt(td*ntr)/(data[i, 7]*(1 - data[i, 11])) # Eq.(12) from the paper
    eta_nom_bt_6_cameras[i, :] = gamma_factor_significance*dback*data[i, 17:27]*np.sqrt(td*ntr)/(data[i, 148]*(1 - data[i, 11]))
    eta_ext_bt_24_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 4] * (1 - data_ext[i, 13])) # Eq.(18) from the paper
    eta_ext_bt_6_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 44] * (1 - data_ext[i, 13]))
    #eta_ext_2_pix_bt_24_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 126:136]*np.sqrt(td*ntr)/(data_ext[i, 238] * (1 - data_ext[i, 248]))
    #eta_ext_2_pix_bt_6_cameras[i, :] = gamma_factor_significance*dback*data_ext[i, 126:136]*np.sqrt(td*ntr)/(data_ext[i, 156] * (1 - data_ext[i, 248]))
    delta_obs[i,:] = dback*SPRK10_first[i,:] # observed transit depth
    #delta_int = delta_obs[i,:]/(1. -data_nommask[i,9] ) # inferred intrinsic transit depth
    delta_obs_ext[i,:] = dback*data_ext[i,14:24] # observed transit depth
    delta_obs_ext_6_cameras[i, :] = dback*data_ext[i,14:24] # observed transit depth with 6 cameras
    #delta_obs_ext_2_pix_24_cameras[i, :] = dback*data_ext[i, 126:136]
    #delta_obs_ext_2_pix_6_cameras[i, :] = dback*data_ext[i, 126:136]
    delta_int = delta_obs_t[i]/ (1 - data[i, 11])
    nbad_sp[i] = np.sum( (eta_nom_bt_24_cameras[i,:]>7.1) & (delta_int<4*84. ))

# Save the eta_nom_bt_24_camerasarray
np.save(dataDIR+'eta_bt_24_cameras.npy', eta_nom_bt_24_cameras)
np.save(dataDIR+'eta_nom_bt_6_cameras', eta_nom_bt_6_cameras)
np.save(dataDIR+'eta_ext_bt_24_cameras.npy', eta_ext_bt_24_cameras)
np.save(dataDIR+'eta_ext_bt_6_cameras.npy', eta_ext_bt_6_cameras)

"""
Now we will obtain a plot showing the amount of false positives detected by Bray et al 2 x 2 mask in comparison with the
amount of false positives detected by our Nominal mask
"""
n_bad_bray_p5 = n_bad_bray[mask_p5]

plt.figure(0)
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Number of potential N_bad', fontsize=fsize)
plt.ylabel('Fraction of targets [%] ', fontsize=fsize)

plt.figure(1)
# Now we plot a percentage histogram like the one presented by Marchiori
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8,
         label='Nominal Mask', alpha=0.5)
plt.hist(n_bad_bray_p5, bins=bins, weights=[1 / len(n_bad_bray_p5)] * len(n_bad_bray_p5), edgecolor='black', rwidth=0.8,
         label='Bray 2 x 2 Mask', alpha=0.5)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Number of potential N_bad', fontsize=fsize)
plt.ylabel('Fraction of targets [%]', fontsize=fsize)
plt.legend()

"""
Now we obtain the NSR for both Bray et al 2 x 2 and our nominal masks as a function of the Target magnitude
"""
plt.figure(2)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    nsr_nominal = np.median(nsr1h[m])
    nsr_bray = np.median(nsr1h_bray[m])
    plt.scatter(Pi, nsr_nominal, color='black')
    plt.scatter(Pi, nsr_bray, color='orange')

plt.xlabel(" P Magnitude", fontsize=fsize)
plt.ylabel(r"$NSR_{1hr}[ppm \sqrt{hr}]$", fontsize=fsize)
plt.title("NSR for Bray and Marchiori nominal masks")

"""
Now we obtain several plots for showing the average size of the Nominal, Secondary and Extended Masks as a function of 
the target magnitude
"""
plt.figure(3)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    size_nominal_mask = np.mean(size_nom[m])
    size_extended_mask = np.mean(size_e[m])
    size_secondary_mask = np.mean(size_sec[m])
    plt.plot(Pi, size_nominal_mask, 'ko', markersize=8)
    plt.plot(Pi, size_secondary_mask, 'rP', markersize=8)
    plt.plot(Pi, size_extended_mask, 'b^', markersize=8)
    plt.legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'])

plt.xlabel(" P Magnitude", fontsize=fsize)
plt.ylabel(r"Average mask size [pixels]", fontsize=fsize)
#plt.title("Average mask size")

"""
Now we plot the size of every mask as a function of the target P magnitude
"""
plt.figure(4)
plt.plot(mag, size_nom, 'k+', label='Nominal Mask')
plt.plot(mag, size_sec, 'r+', label='Secondary Mask')
#plt.plot(mag, size_sec_2, 'g+', label='Secondary Mask (1 pixel ring)')
plt.plot(mag, size_e, 'b+', label='Extended mask')
#plot(mag, size_e_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, size_e_3, 'm+', label='extended mask (3)')
plt.legend()
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel(r'Mask size', fontsize=fsize)


"""
Now we plot the cumulative or total number of nominal mask shapes to address the total number of target stars 
"""
plt.figure(5)
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
    
    # Plot individual mask shapes (use their lengths to plot counts)
    plt.plot(Pi, len(unique_nominal), 'ko', label='Nominal Mask' if i == 0 else "")
    plt.plot(Pi, len(unique_sec), 'rP', label='Secondary Mask' if i == 0 else "")
    plt.plot(Pi, len(unique_ext), 'b^', label='Extended Mask' if i == 0 else "")

# Plot the cumulative total without connecting the dots
plt.plot(P_values, cumulative_total, marker='s', linestyle='', color='orange', label='Three masks combined')

# Add legend and labels
plt.legend(loc='best')
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Cum. count of mask shapes', fontsize=fsize)
print('Cumulative count of unique mask shapes:', cumulative_total)

"""
Now we plot the NSR over 1h  for every mask as a function of the target P magnitude
"""
plt.figure(6)
plt.plot(mag, nsr1h, 'k+', label='Nominal mask')
plt.plot(mag, nsr1h_sec, 'r+', label='Secondary Mask')
plt.plot(mag, nsr1h_ext, 'b+', label='Extended Mask')
plt.semilogy()
plt.legend()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$NSR_{1hr} [ppm \sqrt{hr}]$', fontsize=fsize)
plt.title(r'$NSR_{1hr}$ for every mask')

"""
Now we plot the COB shift as a function of the target P magnitude
"""
plt.figure(7)
plt.plot(mag_2d, sigma_cob_10first_24_cameras, 'ko', markersize=2, alpha=0.5)
plt.plot(mag_2d, sigma_cob_ext_10first_24_cameras, 'b^', markersize=2, alpha=0.5)
plt.plot(mag, sigma_cob_sec_24_cameras, 'rP', markersize = 2, alpha = 0.5)
plt.plot([], [], 'ko', label='nominal mask')
plt.plot([], [], 'b^', label='extended mask')
plt.plot([], [], 'rP', label='secondary mask')
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma^{1h, N_{T}}$', fontsize=fsize)
plt.semilogy()
plt.legend()


"""
Now we plot the comparison between the flux and COB shift methods as a function of the target P magnitude
"""
plt.figure(8)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (cd[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    ext = (eficiency_extended_mask_highest_spr_contaminant[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    plt.plot(Pi, ext, 'b+')
    plt.plot(Pi, eff_cob, 'k^')
    plt.plot(Pi, eff_cob_sec, 'r^')
    plt.plot(Pi, eff_cob_ext, 'b^')
    plt.legend(['Sec. Mask Flux', 'Ext. Mask Flux', 'Nom. Mask COB shift',
            'Sec. Mask COB shift', 'Ext. Mask COB shift'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison of Every Method for Every Mask', fontsize=fsize)


"""
Now we plot the two expressions for the error of the COB (the ones )
"""
plt.figure(9)
plt.plot(mag, sigma_cob, 'ko')
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma_{\Delta C} [pixel]$', fontsize=fsize)
plt.title(r'COB error expression that does not dependend on $\delta_{back}$')



""" 
Comparing double-aperture photometry with (nominal) COB shift
"""
plt.figure(10)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (cd[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    ext = (eficiency_extended_mask_highest_spr_contaminant[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    plt.plot(Pi, ext, 'b+')
    plt.legend(['Sec. Mask Flux', 'Ext. Mask Flux'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison for the two types of DAP', fontsize=fsize)

"""
Now we plot the comparison between extended mask and the correct version of it as a function of the target P magnitude
"""
plt.figure(11)
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
    plt.fill_between([9, 11.7], [60, 60], [100, 100], color='aqua', alpha=0.1)
    plt.fill_between([11, 13.4], [60, 60], [100, 100], color='plum', alpha=0.1)

    # Update label tracking
    labels_added['sec_24'] = True
    labels_added['sec_6'] = True
    labels_added['ext_24'] = True
    labels_added['ext_6'] = True
    labels_added['ext_2_pix_24'] = True
    labels_added['ext_2_pix_6'] = True

    # Connecting lines for the previous data points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_sec, eff_sec], color='purple', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_sec_6_cameras, eff_sec_6_cameras], color='green', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall, eff_ext_overall_24_cameras], color='blue', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall_6_cameras, eff_ext_overall_6_cameras], color='red', linestyle='-', markersize=0)
    
    # Update previous values
    prev_Pi, prev_eff_sec, prev_eff_sec_6_cameras, prev_eff_ext_overall, prev_eff_ext_overall_6_cameras = Pi, eff_sec, eff_sec_6_cameras, eff_ext_overall_24_cameras, eff_ext_overall_6_cameras

# Additional plot settings
plt.vlines(11.7, ymin=60, ymax=100, linestyles='dashed', colors='green')
plt.vlines(11, ymin=60, ymax=100, linestyles='dashdot', colors='red')
plt.ylim(60, 100)
plt.xlim(9.9, 13.1)
plt.text(10, 78.1, 'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold')
plt.text(11.2, 75, 'On-board light curve processing region', color='red', weight='bold')

# Display legend below the plot
plt.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', borderaxespad=0., ncol=2)
plt.ylim(60, 100)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)

# Adjust layout to accommodate the legend below
plt.tight_layout(rect=[0, 0, 1, 0.88])
#plt.title('Double-Aperture Photometry Comparison', fontsize=fsize)

"""
Now we plot the statistical significance for every mask as a function of the target P magnitude
"""
plt.figure(12)

plt.plot(mag, eta_t_6_cameras_earth_like, 'd', markersize=4, label='Jovian planets (6 cameras)')  # Diamond first
plt.plot(mag, eta_t_6_cameras_superearth, '*', markersize=4, label='Super-Earth (6 cameras)')      # Star second
plt.plot(mag, eta_t_24_cameras, 'o', markersize=4, label='Earth-like planets (24 cameras)')        # Circle third

# Horizontal line and other settings
plt.hlines(7.1, xmin=10, xmax=13, linestyles='dashed', colors='red')
plt.semilogy()
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel(r'$ \eta $', fontsize=fsize)
plt.legend()
plt.xlim(10, 13)

"""
Now we plot the efficiency of the COB shift (all contaminants)
    
"""
plt.figure(13)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    s_24_cameras = (eta_nom_bt_24_cameras>flux_thresh_nom_mask)[m,:].sum()
    s_6_cameras = (eta_nom_bt_6_cameras>flux_thresh_nom_mask)[m,:].sum()
    eff_ext_cob_overall = ((eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras>flux_thresh_nom_mask))[m,:].sum() / s_24_cameras * 100.
    eff_ext_cob_overall_6_cameras = ((eta_cob_ext_10first_6_cameras > cob_thresh) &  (eta_nom_bt_6_cameras>flux_thresh_nom_mask))[m,:].sum() / s_6_cameras * 100.
    eff_cob = ((eta_cob_nom_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras>flux_thresh_nom_mask))[m,:].sum() / s_24_cameras * 100.
    eff_cob_6_cameras = ((eta_cob_nom_10first_6_cameras > cob_thresh) & (eta_nom_bt_6_cameras>flux_thresh_nom_mask))[m,:].sum() / s_6_cameras * 100.
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_single_contaminant_24_cameras[m].sum()) * 100
    eff_cob_sec_6_cameras = (secondary_mask_conditions_cob_6_cameras[m].sum() / fp_single_contaminant_6_cameras[m].sum()) * 100
    
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
    plt.fill_between([9, 11.7], [60, 60], [100, 100], color='aqua', alpha=0.1)
    plt.fill_between([11, 13.4], [60,60], [100, 100], color='plum', alpha=0.1)
    
    # Plot lines connecting the points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall, eff_ext_cob_overall], color='blue', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob, eff_cob], color='orange', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec, eff_cob_sec], color='purple', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall_6_cameras, eff_ext_cob_overall_6_cameras], color='red', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_6_cameras, eff_cob_6_cameras], color='olive', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec_6_cameras, eff_cob_sec_6_cameras], color='green', linestyle='-')   
    # Update previous values
    prev_Pi, prev_eff_ext_cob_overall, prev_eff_cob, prev_eff_cob_sec, prev_eff_ext_cob_overall_6_cameras, prev_eff_cob_6_cameras, prev_eff_cob_sec_6_cameras = Pi, eff_ext_cob_overall, eff_cob, eff_cob_sec, eff_ext_cob_overall_6_cameras, eff_cob_6_cameras, eff_cob_sec_6_cameras


    #plt.legend(['Ext. Mask (10 contaminants)'], loc='best')
    plt.vlines(11.7, ymin=60, ymax = 100, linestyles='dashed', colors='green')
    plt.vlines(11, ymin=60, ymax=100, linestyles='dashdot', colors='red')
    plt.ylim(60, 100)
    plt.xlim(9.9, 13.1)
    plt.text(10, 76.1,'Earth-like planet detection \nregion (24 cameras)', color='green', weight='bold')
    plt.text(11.2, 72, 'On-board light curve processing region', color='red', weight='bold')
       
# Move the legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.21), borderaxespad=0., fancybox=True, ncol=2)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)
# Adjust layout to accommodate the legend below
plt.tight_layout(rect=[0, 0, 1, 0.88])


nfp = (eta_nom_bt_24_cameras > flux_thresh_nom_mask)
nfp_ext_mask = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) & (delta_obs_ext > delta_obs + depth_sig_scaling*sig_depth_24_cameras_10first)
nfp_ext_flux_without_significant_transit_depth_condition = (eta_ext_bt_24_cameras> flux_thresh_ext_mask) &  (delta_obs_ext > delta_obs)
nfp_ext_mask_single_contaminant = (eficiency_extended_mask_highest_spr_contaminant)
nfp_nom_cob = (eta_cob_nom_10first_24_cameras > cob_thresh)
nfp_ext_cob = (eta_cob_ext_10first_24_cameras > cob_thresh)
nfp_sec_mask = (secondary_mask_conditions_24_cameras)

nfp_highest_contaminant = (eta_nom_bt_24_cameras[:,0]>flux_thresh_nom_mask)
nfp_ext_mask_highest_contaminant = (eta_ext_bt_24_cameras[:,0]>flux_thresh_ext_mask) & (delta_obs_ext[:,0]>delta_obs[:,0]+depth_sig_scaling*sig_depth_24_cameras_10first[:,0])

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

# Initialize errors for fractions
weighted_variance_fraction_fp_ext_cob_no_ext_flux = 0
weighted_variance_fraction_fp_nom_cob_no_ext_flux = 0
weighted_variance_fraction_fp_ext_flux_no_nom_cob = 0
weighted_variance_fraction_fp_ext_cob_no_nom_cob = 0

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
    eff_ext_flux = (nfp[m] & nfp_ext_mask[m]).sum() / nfp[m].sum() * 100.
    eff_sec_flux = (nfp_sec_mask[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100
    eff_ext_flux_single_contaminant = (nfp_ext_mask_single_contaminant[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100.
    eff_ext_flux_without_significant_transit_depth_condition = (nfp[m] & nfp_ext_flux_without_significant_transit_depth_condition[m]).sum() / nfp[m].sum() * 100.
    eff_nom_cob = (nfp[m] & nfp_nom_cob[m]).sum() / nfp[m].sum() * 100.
    eff_ext_cob = (nfp[m] & nfp_ext_cob[m]).sum() / nfp[m].sum() * 100.
    eff_sec_cob = (secondary_mask_conditions_cob_24_cameras[m]).sum() / fp_single_contaminant_24_cameras[m].sum() * 100.
    fraction_fp_ext_cob_no_ext_flux = ((nfp_ext_mask[m]==False) & nfp_ext_cob[m] & nfp[m] ).sum()/nfp[m].sum()
    fraction_fp_nom_cob_no_ext_flux = ((nfp_ext_mask[m]==False) & nfp_nom_cob[m] & nfp[m] ).sum()/nfp[m].sum()
    fraction_fp_ext_flux_no_nom_cob = ((nfp_nom_cob[m]==False) & nfp_ext_mask[m] & nfp[m] ).sum()/nfp[m].sum()
    fraction_fp_ext_cob_no_nom_cob = ((nfp_nom_cob[m]==False) & nfp_ext_cob[m] & nfp[m] ).sum()/nfp[m].sum()

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


    # Accumulate weighted variance
    weighted_variance_fraction_fp_ext_cob_no_ext_flux += weights[i] * variance_fraction_fp_ext_cob_no_ext_flux
    weighted_variance_fraction_fp_nom_cob_no_ext_flux += weights[i] * variance_fraction_fp_nom_cob_no_ext_flux
    weighted_variance_fraction_fp_ext_flux_no_nom_cob += weights[i] * variance_fraction_fp_ext_flux_no_nom_cob
    weighted_variance_fraction_fp_ext_cob_no_nom_cob += weights[i] * variance_fraction_fp_ext_cob_no_nom_cob

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

# Print weighted efficiency results with errors
print(f'Weighted extended flux efficiency: {weighted_eff_ext_flux:.2f}% ± {weighted_error_eff_ext_flux:.2f}%')
print(f'Weighted secondary flux efficiency: {weighted_eff_sec_flux:.2f}% ± {weighted_error_eff_sec_flux:.2f}%')
print(f'Weighted extended flux for a single contaminant efficiency: {weighted_eff_ext_flux_single_contaminant:.2f}% ± {weighted_error_eff_ext_flux_single_contaminant:.2f}%')
print(f'Weighted extended flux withouth the significant transit condition: {weighted_eff_ext_flux_without_significant_transit_depth_condition:.2f}% ± {weighted_error_eff_ext_flux_without_significant_transit_depth_condition:.2f}%')
print(f'Weighted nominal COB efficiency: {weighted_eff_nom_cob:.2f}% ± {weighted_error_eff_nom_cob:.2f}%')
print(f'Weighted extended COB efficiency: {weighted_eff_ext_cob:.2f}% ± {weighted_error_eff_ext_cob:.2f}%')
print(f'Weighted secondary mask COB efficiency: {weighted_eff_sec_cob:.2f}% ± {weighted_error_eff_sec_cob:.2f}%')

# Print weighted fractions and their errors
print(f'Weighted fraction of FPs detected by ECOB but not by EFX: {weighted_fraction_fp_ext_cob_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by NCOB but not by EFX: {weighted_fraction_fp_nom_cob_no_ext_flux * 100:.2f}% ± {weighted_error_fraction_fp_nom_cob_no_ext_flux * 100:.2f}%')
print(f'Weighted fraction of FPs detected by EFX but not by NCOB: {weighted_fraction_fp_ext_flux_no_nom_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_flux_no_nom_cob * 100:.2f}%')
print(f'Weighted fraction of FPs detected by ECOB but not by NCOB: {weighted_fraction_fp_ext_cob_no_nom_cob * 100:.2f}% ± {weighted_error_fraction_fp_ext_cob_no_nom_cob * 100:.2f}%')
"""
The condition of eta_ext
"""
plt.figure(14)
mag_value = 12.5

# Lists to store data for plotting histograms
count_eta_ext = []
count_eta_ext_6_cameras = []
magnitude_ranges = []
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #if Pi == mag_value:
    condition_eta_ext = (delta_obs_ext>delta_obs) & (eta_ext_bt_24_cameras>flux_thresh_nom_mask)
    filtered_eta_ext = condition_eta_ext[m, :]
    count_ext = filtered_eta_ext.sum()
    count_eta_ext.append(count_ext) 
    condition_eta_ext_6_cameras = (delta_obs_ext>delta_obs) & (eta_ext_bt_6_cameras>flux_thresh_nom_mask)
    filtered_eta_ext_6_cameras = condition_eta_ext_6_cameras[m, :]
    count_ext_6_cameras = filtered_eta_ext_6_cameras.sum()
    count_eta_ext_6_cameras.append(count_ext_6_cameras)
    magnitude_ranges.append(str(Pi))

plt.bar(magnitude_ranges, count_eta_ext, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{ext} \geq \eta_{min})$ (24 cameras)')
plt.bar(magnitude_ranges, count_eta_ext_6_cameras, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{ext} \geq \eta_{min})$ (6 cameras)')
#plt.hist(condition_eta_ext_list, label='eta_ext', alpha=0.5)
#plt.hist(condition_eta_ext_6_cameras_list, label='eta_ext_6_cameras', alpha=0.5)
plt.legend()  
plt.xlabel('Magnitude Bin')
plt.ylabel('Counts')
plt.legend()


plt.figure(15)
mag_value = 12.5

# Lists to store data for plotting histograms
count_eta_nom = []
count_eta_nom_6_cameras = []
magnitude_ranges = []
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #if Pi == mag_value:
    condition_eta_nom = (delta_obs_ext>delta_obs) & (eta_nom_bt_24_cameras>flux_thresh_nom_mask)
    filtered_eta_nom = condition_eta_nom[m, :]
    count_nom = filtered_eta_nom.sum()
    count_eta_nom.append(count_nom) 
    condition_eta_nom_6_cameras = (delta_obs_ext>delta_obs) & (eta_nom_bt_6_cameras>flux_thresh_nom_mask)
    filtered_eta_nom_6_cameras = condition_eta_nom_6_cameras[m, :]
    count_nom_6_cameras = filtered_eta_nom_6_cameras.sum()
    count_eta_nom_6_cameras.append(count_nom_6_cameras)
    magnitude_ranges.append(str(Pi))
        
  

plt.bar(magnitude_ranges, count_eta_nom, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$ (24 cameras)')
plt.bar(magnitude_ranges, count_eta_nom_6_cameras, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$ (6 cameras)')
#plt.hist(condition_eta_ext_list, label='eta_ext', alpha=0.5)
#plt.hist(condition_eta_ext_6_cameras_list, label='eta_ext_6_cameras', alpha=0.5)
plt.legend()  
plt.xlabel('Magnitude Bin')
plt.ylabel('Counts')
plt.legend()

"""
The histogram of eta_ext
"""
plt.figure(16)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt_24_cameras[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

eta_and_delta_mask = (eta_nom_bt_24_cameras[mag_mask] >= flux_thresh_nom_mask) & (delta_obs_ext[mag_mask] >= delta_obs[mag_mask])
eta_and_delta_mask_6_cameras = (eta_nom_bt_6_cameras[mag_mask] >= flux_thresh_nom_mask) & (delta_obs_ext[mag_mask] >= delta_obs[mag_mask])

final_eta_ext = eta_ext[eta_and_delta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[eta_and_delta_mask_6_cameras]

number_below_threshold = len(np.where(final_eta_ext < flux_thresh_nom_mask)[0])
number_above_threshold = len(np.where(final_eta_ext > flux_thresh_nom_mask)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=150, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 100, number_below_threshold, color='blue', weight='bold')
plt.text(50, 100, number_above_threshold, color='blue', weight='bold')
plt.text(0, 70, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 70, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given $(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$')
plt.ylabel('Counts')
plt.ylim(0, 150)
plt.legend()


plt.figure(17)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_nom_24_cameras = eta_nom_bt_24_cameras[mag_mask]
eta_nom_6_cameras = eta_nom_bt_6_cameras[mag_mask]

eta_nom_24_cameras_flat = eta_nom_24_cameras.flatten()

eta_nom_6_cameras_flat = eta_nom_6_cameras.flatten()


number_below_threshold = len(np.where(eta_nom_24_cameras_flat < 7.1)[0])
number_above_threshold = len(np.where(eta_nom_24_cameras_flat > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(eta_nom_6_cameras_flat < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(eta_nom_6_cameras_flat > 7.1)[0])

plt.hist(eta_nom_24_cameras_flat, bins=15, range=(0, 50), alpha=0.5, label='24 cameras')
plt.hist(eta_nom_6_cameras_flat, bins=15, range=(0,50), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=19800, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
plt.text(0, 17500, number_below_threshold, color='blue', weight='bold')
plt.text(50, 17500, number_above_threshold, color='blue', weight='bold')
plt.text(0, 12500, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 12500, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{k}^{nom}$')
plt.legend()

"""
The histogram of eta_ext
"""
plt.figure(18)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt_24_cameras[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

eta_mask = (eta_nom_bt_24_cameras[mag_mask] >= flux_thresh_nom_mask)
eta_mask_6_cameras = (eta_nom_bt_6_cameras[mag_mask] >= flux_thresh_nom_mask)

final_eta_ext = eta_ext[eta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[eta_mask_6_cameras]

number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=150, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 100, number_below_threshold, color='blue', weight='bold')
plt.text(50, 100, number_above_threshold, color='blue', weight='bold')
plt.text(0, 70, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 70, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given $(\eta_{nom} \geq \eta_{min})$')
plt.ylabel('Counts')
plt.ylim(0, 150)
plt.legend()


"""
The histogram of eta_ext
"""
plt.figure(19)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt_24_cameras[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

delta_mask = (delta_obs_ext[mag_mask] >= delta_obs[mag_mask])

final_eta_ext = eta_ext[delta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[delta_mask]

number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=22000, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 1000, number_below_threshold, color='blue', weight='bold')
plt.text(50, 1000, number_above_threshold, color='blue', weight='bold')
plt.text(0, 7000, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 7000, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given  $(\delta_{back}^{ext} \geq \delta_{back}^{nom})$')
plt.ylabel('Counts')
#plt.ylim(0, 200)
plt.legend()

plt.figure(20)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt_24_cameras[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

final_eta_ext = eta_ext.flatten()
final_eta_ext_6_cameras = eta_ext_6_cameras.flatten()


number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=20210, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
plt.text(0, 17000, number_below_threshold, color='blue', weight='bold')
plt.text(50, 17000, number_above_threshold, color='blue', weight='bold')
plt.text(0, 16000, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 16000, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$')
plt.ylabel('Counts')
plt.yscale('log')
#plt.ylim(0, 200)
plt.legend()


plt.figure(21)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt_24_cameras[mag_mask]
eta_ext_cob = eta_cob_ext_10first_24_cameras[mag_mask]
eta_ext_flat = eta_ext.flatten()
eta_ext_cob_flat = eta_ext_cob.flatten()
eta_ext_ratio = eta_ext_flat / eta_ext_cob_flat

eta_ratio_above = len(np.where(eta_ext_ratio > 1)[0])
eta_ratio_below = len(np.where(eta_ext_ratio < 1)[0])

plt.hist(eta_ext_flat/eta_ext_cob_flat, bins=51, range=(0, 10), alpha=0.5)
plt.vlines(1, ymin=0, ymax=2500, linestyles='dashdot', colors='orange')
plt.text(0.1, 2000, eta_ratio_below, color='blue', weight='bold')
plt.text(6, 2000, eta_ratio_above, color='blue', weight='bold')
plt.xlabel(r'$ \eta_{k}^{ext} / \eta^{COB, ext}_{k}$')
plt.ylabel('Counts')

plt.figure(22)
mag_mask_10 = (mag >= 9.75) & (mag <= 10.25)
mag_mask_13 = (mag >= 12.75) & (mag <= 13.25)
eta_sec_mag_10 = eta_c[mag_mask_10]
eta_sec_mag_13 = eta_c[mag_mask_13]
delta_sec_mag_10 = delta_obs_c[mag_mask_10]

#eta_ext_mag_10_flat = eta_ext_mag_10.flatten()
#delta_ext_mag_10_flat = delta_ext_mag_10.flatten()
#eta_ext_mag_13_flat = eta_ext_mag_13.flatten()

number_below_mag_10 = len(np.where(eta_sec_mag_10 < 7.1)[0])
number_above_mag_10 = len(np.where(eta_sec_mag_10 > 7.1)[0])

number_below_mag_13 = len(np.where(eta_sec_mag_13 < 7.1)[0])
number_above_mag_13 = len(np.where(eta_sec_mag_13 > 7.1)[0])

plt.hist(eta_sec_mag_10, bins=51, range=(0, 10), alpha=0.5, label=r'$\eta_{sec}$ (P = 10)')
plt.hist(eta_sec_mag_13, bins=51, range=(0, 10), alpha=0.5, label=r'$\eta_{sec}$ (P = 13)')
plt.vlines(7.1, ymin=0, ymax=12000, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
plt.text(4, 6700, number_below_mag_10, color='blue', weight='bold')
plt.text(9, 6700, number_above_mag_10, color='blue', weight='bold')
plt.text(4, 5600, number_below_mag_13, color='orange', weight='bold')
plt.text(9, 5600, number_above_mag_13, color='orange', weight='bold')
plt.yscale('log')
plt.xlabel(r'$\eta_{ext}$')
plt.ylabel('Counts')
plt.legend()
#plt.hist(delta_ext_mag_10_flat, bins=51, range=(0, 10), alpha=0.5, label=r'$\delta_{ext}$')

delta_obs_sec_mag_10 = delta_obs_c[mag_mask_10]
delta_obs_sec_mag_13 = delta_obs_c[mag_mask_13]


#delta_ext_mag_10_flat = delta_ext_mag_10.flatten()
#delta_ext_mag_13_flat = delta_ext_mag_13.flatten()

l1 = (delta_obs_ext[mag_mask] >= delta_obs[mag_mask]) & np.logical_and(flux_thresh_nom_mask <= eta_ext, eta_ext <= 2*flux_thresh_nom_mask) & (eta_nom_bt_24_cameras[mag_mask] >= 2 *flux_thresh_nom_mask)
l2 = (delta_obs_ext[mag_mask] >= delta_obs[mag_mask]) & np.logical_and(flux_thresh_nom_mask <= eta_nom_bt_24_cameras[mag_mask], eta_nom_bt_24_cameras[mag_mask] <= 2*flux_thresh_nom_mask) & (eta_ext >= 2*flux_thresh_nom_mask)
l3 = (delta_obs_ext[mag_mask] >= delta_obs[mag_mask]) & np.logical_and(flux_thresh_nom_mask <= eta_nom_bt_24_cameras[mag_mask], eta_nom_bt_24_cameras[mag_mask] <= 2*flux_thresh_nom_mask) & np.logical_and(flux_thresh_nom_mask <= eta_ext, eta_ext <= (2*flux_thresh_nom_mask))
l4 = np.logical_and(flux_thresh_nom_mask < eta_nom_bt_24_cameras[mag_mask], eta_nom_bt_24_cameras[mag_mask] < 2*flux_thresh_nom_mask)


print('L1 is:', np.sum(l1))
print('L2 is:', np.sum(l2))
print('L3 is:',np.sum(l3))
print('L4 is:',np.sum(l4))
print('L1 + L2 + L3 is:', np.sum(l1 + l2 + l3))


mask_ext=  (eta_ext_bt_24_cameras > flux_thresh_ext_mask) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (eta_cob_ext_10first_24_cameras > cob_thresh) 
#& (delta_obs_ext > delta_obs)  

mask_nom = (eta_ext_bt_24_cameras > flux_thresh_ext_mask) & (eta_nom_bt_24_cameras > flux_thresh_nom_mask) & (eta_cob_nom_10first_24_cameras > cob_thresh) 
#& (delta_obs_ext > delta_obs)

cob_mask_cases =  (eta_cob_ext_10first_24_cameras > eta_ext_bt_24_cameras) & (eta_ext_bt_24_cameras > flux_thresh_nom_mask) & (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras> flux_thresh_nom_mask)
cob_ext = (eta_ext_bt_24_cameras > eta_cob_ext_10first_24_cameras ) & (eta_ext_bt_24_cameras > flux_thresh_ext_mask) & (eta_cob_ext_10first_24_cameras > cob_thresh) & (eta_nom_bt_24_cameras> flux_thresh_nom_mask)
cob_ext_mask = (eta_cob_ext_10first_24_cameras > cob_thresh)
cob_nom_mask = (eta_cob_nom_10first_24_cameras > cob_thresh)

print('Number of times eta_ext > eta_ext^COB given (eta_ext > eta_min) and (eta_ext^COB > eta^COB_min) and (eta_nom > eta_min):', np.sum(cob_ext))
print('Number of times eta_ext^COB > eta_ext given (eta_ext > eta_min) and (eta_ext^COB > eta^COB_min) and (eta_nom > eta_min):', np.sum(cob_mask_cases))
print('Number of times eta_ext^COB > eta^COB_min:', np.sum(cob_ext_mask))
print('Number of times eta_nom^COB > eta^COB_min:', np.sum(cob_nom_mask))

etas_ext_nom = eta_ext_bt_24_cameras[mask_ext]/eta_nom_bt_24_cameras[mask_ext]
etas =  eta_ext_bt_24_cameras[mask_ext]/eta_cob_ext_10first_24_cameras[mask_ext]
etas_nom = eta_ext_bt_24_cameras[mask_nom] / eta_cob_nom_10first_24_cameras[mask_nom]


# Find the indices where eta_cob_ext is higher than eta_ext
indices_for_when_eta_cob_ext_is_higher_than_eta_ext = np.where(etas < 1)[0]
indices_for_when_eta_cob_ext_is_smaller_than_eta_ext = np.where(etas > 1)[0]

# Find the indices where eta_cob_nom is higher than eta_ext
indices_for_when_eta_cob_nom_is_higher_than_eta_ext = np.where(etas_nom < 1)[0]
indices_for_when_eta_cob_nom_is_smaller_than_eta_ext = np.where(etas_nom > 1)[0]

print(len(indices_for_when_eta_cob_ext_is_higher_than_eta_ext))
print(len(indices_for_when_eta_cob_ext_is_smaller_than_eta_ext))
print(len(indices_for_when_eta_cob_nom_is_higher_than_eta_ext))
print(len(indices_for_when_eta_cob_nom_is_smaller_than_eta_ext))

eta_nom_mag_10 = eta_nom_bt_24_cameras[mag_mask_10]
eta_nom_mag_13 = eta_nom_bt_24_cameras[mag_mask_13]

eta_nom_mag_10_flat = eta_nom_mag_10.flatten()
eta_nom_mag_13_flat = eta_nom_mag_13.flatten()

number_of_stars_below_eta_min_for_eta_nom_mag_10 = len(np.where(eta_nom_mag_10_flat < 7.1)[0])
number_of_stars_below_eta_min_for_eta_nom_mag_13 = len(np.where(eta_nom_mag_13_flat < 7.1)[0])

number_of_stars_above_eta_min_for_eta_nom_mag_10 = len(np.where(eta_nom_mag_10_flat > 7.1)[0])
number_of_stars_above_eta_min_for_eta_nom_mag_13 = len(np.where(eta_nom_mag_13_flat > 7.1)[0])


plt.figure(23)
n_sec_fp_mag_10 = np.sum((eta_c[mag_mask_10] > flux_thresh_sec_mask) & (delta_obs_c[mag_mask_10]>delta_obs_t[mag_mask_10] + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_24_cameras[mag_mask_10]**2 + sig_depth_nominal_mask_24_cameras[mag_mask_10]**2)) & (eta_t[mag_mask_10]>flux_thresh_nom_mask))
n_fp_mag_10 = np.sum((eta_t[mag_mask_10] > flux_thresh_nom_mask))

n_sec_fp_mag_13 = np.sum((eta_c[mag_mask_13] > flux_thresh_sec_mask) & (delta_obs_c[mag_mask_13]>delta_obs_t[mag_mask_13] + depth_sig_scaling*np.sqrt(sig_depth_secondary_mask_6_cameras[mag_mask_10]**2 + sig_depth_nominal_mask_6_cameras[mag_mask_10]**2)) & (eta_t[mag_mask_13]>flux_thresh_nom_mask))
n_fp_mag_13 = np.sum((eta_t[mag_mask_13] > flux_thresh_nom_mask))

# Values to plot
colors = ['blue', 'green', 'red', 'purple']
values = [n_sec_fp_mag_10, n_sec_fp_mag_13, n_fp_mag_10, n_fp_mag_13]
labels = [r'$N_{FP}^{sec}(mag = 10)$', r'$N_{FP}^{sec} (mag = 13)$', r'$N_{FP}(mag=10)$', r'$N_{FP}(mag=13)$']

# Plotting the histogram
bars = plt.bar(labels, values, color=['blue', 'green', 'red', 'purple'])
# Adding the value of each bar on top
for bar, value, color in zip(bars, values, colors):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.05, value, ha='center', va='bottom', color=color)
plt.xlabel('Conditions')
plt.ylabel('Count')
plt.ylim(0, 3700)
plt.title(r'Histogram of $N_{sec}^{FP}$ and $N_{FP}$')


plt.figure(24)

mask_for_eta_ext_given_delta_and_eta_nom_conditions = (delta_obs_ext > delta_obs) & (eta_nom_bt_24_cameras> flux_thresh_nom_mask)
eta_ext_given_delta_and_eta_nom_conditions = eta_ext_bt_24_cameras[mask_for_eta_ext_given_delta_and_eta_nom_conditions]
eta_bt_given_delta_and_eta_nom_conditions = eta_nom_bt_24_cameras[mask_for_eta_ext_given_delta_and_eta_nom_conditions]

eta_ratio = eta_ext_given_delta_and_eta_nom_conditions/eta_bt_given_delta_and_eta_nom_conditions

above_one = np.sum(eta_ratio > 1)
below_one = np.sum(eta_ratio < 1)
print(above_one)
print(below_one)

plt.plot(mag_2d[mask_for_eta_ext_given_delta_and_eta_nom_conditions], eta_ratio, 'bo', markersize=2, alpha=0.65)
plt.hlines(1, xmin=10,xmax=13, linestyles='dashdot', colors='red')
plt.text(11, 950, f'Above one: {above_one}', color='red', fontsize=12, weight='bold')
plt.text(11, 0.02, f'Below one: {below_one}', color='red', fontsize=12, weight='bold')
plt.xlabel('P magnitude', fontsize=12)
plt.ylabel(r'$ \eta_{ext} / \eta_{nom} $ given $(\eta_{nom} > \eta_{min})$ and $(\delta_{ext} > \delta_{nom})$')
#plt.plot(mag_2d, eta_bt)
plt.xlim(10,13)
plt.yscale('log')
plt.legend()

mask_for_delta_ext_given_eta_ext_and_eta_nom_conditions = (eta_ext_bt_24_cameras> flux_thresh_nom_mask) & (eta_nom_bt_24_cameras> flux_thresh_nom_mask)
delta_ext_given_eta_ext_and_eta_nom_conditions = delta_obs_ext[mask_for_delta_ext_given_eta_ext_and_eta_nom_conditions]

mask_for_eta_nom_given_eta_ext_and_delta_ext_conditions = (eta_ext_bt_24_cameras> flux_thresh_nom_mask) & (delta_obs_ext > delta_obs)
eta_nom_given_eta_ext_and_delta_ext_conditions = eta_nom_bt_24_cameras[mask_for_eta_nom_given_eta_ext_and_delta_ext_conditions]

mask_for_eta_nom_higher_than_eta_min = (eta_nom_bt_24_cameras> flux_thresh_nom_mask)
eta_nom_higher_than_eta_min = eta_nom_bt_24_cameras[mask_for_eta_nom_higher_than_eta_min]


eta_nom_flat = eta_nom_bt_24_cameras.flatten()
eta_ext_flat = eta_ext_bt_24_cameras.flatten()
eta_ext_cob_flat = eta_cob_ext_10first_24_cameras.flatten()
eta_nom_cob_flat = eta_cob_nom_10first_24_cameras.flatten()
mag_2d_flat = mag_2d.flatten()



indices_fp_eta_cob_nom = np.where(eta_nom_bt_24_cameras > 7.1)
eta_cob_nom_fp = eta_cob_nom_10first_24_cameras[indices_fp_eta_cob_nom]
#eta_cob_nom_fp = eta_cob_nom_10first_24_cameras[indice_eta_cob_nom_below_3]
#eta_cob_nom_fp_cleaned = eta_cob_nom_fp[~np.isnan(eta_cob_nom_fp)]
index_below_3_nom  = (eta_cob_nom_fp < 3)
eta_cob_nom_fp_cleaned = eta_cob_nom_fp[index_below_3_nom]

indices_fp_eta_cob_ext = np.where(eta_nom_bt_24_cameras > 7.1)
#eta_cob_ext_fp = eta_cob_ext_10first_24_cameras[indices_eta_cob_ext_below_3]
eta_cob_ext_fp = eta_cob_ext_10first_24_cameras[indices_fp_eta_cob_ext]
index_below_3_ext = (eta_cob_ext_fp < 3)
eta_cob_ext_fp = eta_cob_ext_fp[index_below_3_ext]
mag_ext = mag_2d[indices_fp_eta_cob_nom]
mag_etas_ratio = mag_ext[index_below_3_nom]

trim_length_for_ratio_etas = min(len(mag_etas_ratio), len(eta_cob_nom_fp_cleaned), len(eta_cob_ext_fp))

eta_cob_nom_fp_trimmed = eta_cob_nom_fp_cleaned[:trim_length_for_ratio_etas]
eta_cob_ext_fp_trimmed = eta_cob_ext_fp[:trim_length_for_ratio_etas]
mag_trimmed = mag_etas_ratio[:trim_length_for_ratio_etas]

sigma_ext = sigma_cob_ext_10first_24_cameras[indices_fp_eta_cob_ext]
sigma_nom = sigma_cob_10first_24_cameras[indices_fp_eta_cob_nom]

gamma_ext = gamma_cob_ext_10first_24_cameras[indices_fp_eta_cob_ext]
gamma_nom = gamma_cob_nom_10first_24_cameras[indices_fp_eta_cob_nom]
mag_fp = mag_2d[indices_fp_eta_cob_nom]

# NowTrim the arrays to match the length of the shortest array (mag in this case)
trim_length = min(len(mag), len(eta_cob_nom_fp_cleaned), len(eta_cob_ext_fp))

count_below_3_nom = np.sum(eta_cob_nom_fp < 3)
count_below_3_ext = np.sum(eta_cob_ext_fp < 3)


indices_eta_cob_ext_below_3 = np.where((eta_cob_ext_10first_24_cameras < 3) & (eta_nom_bt_24_cameras > 7.1))
indice_eta_cob_nom_below_3 = np.where((eta_cob_nom_10first_24_cameras < 3) & (eta_nom_bt_24_cameras > 7.1))

sigma_cob_ext = sigma_cob_ext_10first_24_cameras[indices_eta_cob_ext_below_3]
sigma_cob_nom = sigma_cob_10first_24_cameras[indice_eta_cob_nom_below_3]

gamma_cob_ext = gamma_cob_ext_10first_24_cameras[indices_eta_cob_ext_below_3]
gamma_cob_nom = gamma_cob_nom_10first_24_cameras[indice_eta_cob_nom_below_3]

trim_length_sigma = min(len(mag), len(sigma_cob_nom), len(sigma_cob_ext))

sigma_cob_nom_trimmed = sigma_cob_nom[:trim_length_sigma]
sigma_cob_ext_trimmed = sigma_cob_ext[:trim_length_sigma]
mag_trimmed_sigma = mag[:trim_length_sigma]

trim_length_gamma = min(len(mag), len(gamma_cob_nom), len(gamma_cob_ext))
gamma_cob_ext_trimmed = gamma_cob_ext[:trim_length_gamma]
gamma_cob_nom_trimmed = gamma_cob_nom[:trim_length_gamma]

# Printing the results
print(f'Number of elements in eta_cob_nom_fp_trimmed below 3: {count_below_3_nom}')
print(f'Number of elements in eta_cob_ext_fp_trimmed below 3: {count_below_3_ext}')

plt.figure(25)
#plt.plot(mag_etas_ratio, eta_cob_ext_fp_trimmed/eta_cob_nom_fp_trimmed, 'bo', markersize = 5)
plt.plot(mag_2d[eta_nom_bt_24_cameras > 7.1], eta_cob_ext_10first_24_cameras[eta_nom_bt_24_cameras > 7.1]/ eta_cob_nom_10first_24_cameras[eta_nom_bt_24_cameras> 7.1], 'bo', markersize=5)
#plt.plot(mag_trimmed, eta_cob_nom_fp_trimmed, 'k+')
plt.hlines(1, xmin=9, xmax=14, linestyles='dashdot', colors='red')
#plt.ylim(0.1, 10)
plt.xlim(9.5, 13.5)
plt.semilogy()
plt.ylabel(r'$ \eta_{k}^{ext, \Delta C} / \eta_{k}^{nom, \Delta C}$', fontsize=fsize)
plt.xlabel('P mag', fontsize=fsize)
# Plot a single point for the legend with labels
#plt.plot([], [], 'b^', label='Ext. mask')
#plt.plot([], [], 'k+', label='Nom. mask')
print('the median of eta_cob_ext is:', np.median(eta_cob_ext_fp_trimmed))
print('the median of eta_cob_nom is:', np.median(eta_cob_nom_fp_trimmed))
plt.legend()

plt.figure(26)
plt.plot(mag_2d, nsr1h_ext, 'ko')
plt.plot(mag_2d, nsr1h, 'gP')

plt.plot([], [], 'ko', label='nsr_ext')
plt.plot([], [], 'gP', label='nsr_nom')

plt.legend()

plt.figure(27)
plt.plot(mag_2d, sigma_cob_ext_10first_24_cameras, 'b^', markersize = 3)
plt.plot(mag_2d, sigma_cob_10first_24_cameras, 'ko', markersize=3)
#plt.ylim(0.1, 10)
#plt.xlim(9.5, 13.5)
plt.semilogy()
#plt.hlines(1, xmin=9, xmax=14, linestyles='dashdot', colors='red')
plt.ylabel(r'$ \sigma_{k}^{1h, N_{T}} [pix]$', fontsize=fsize)
plt.xlabel('P Magnitude', fontsize=fsize)
plt.plot([], [], 'b^', label='Ext. mask')
plt.plot([], [], 'ko', label='Nom. mask')
plt.legend()


plt.figure(29)
plt.plot(mag_fp, gamma_ext/sigma_ext, 'b^')
plt.plot(mag_fp, gamma_nom/sigma_nom, 'k+')
plt.hlines(3, xmin=9, xmax=14, linestyles='dashdot', colors='red', label=r'$\eta_{min}^{\Delta C} = 3$')
# Plot a single point for the legend with labels
plt.plot([], [], 'b^', label='Ext. mask')
plt.plot([], [], 'k+', label='Nom. mask')
plt.xlim(9.5, 13.5)
plt.semilogy()
plt.ylabel(r'$ \Gamma_{k} / \sigma^{1h, NT}$', fontsize=fsize)
plt.legend()


plt.figure(30)
plt.plot(mag[eta_t > flux_thresh_nom_mask], eta_c[eta_t > flux_thresh_nom_mask], 'rP', markersize=5)
#plt.plot(mag[eta_t > flux_thresh_nom_mask], nsr1h_sec[eta_t > flux_thresh_nom_mask], 'ko', markersize=5)
plt.semilogy()
plt.ylabel(r'$\eta_{kmax}^{sec}$ for cases where $(\eta_{nom} > \eta_{min})$', fontsize=fsize)
plt.xlabel('P mag', fontsize=fsize)

plt.show()