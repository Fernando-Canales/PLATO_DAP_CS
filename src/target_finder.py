"""
Script that finds a target given specific conditions for the
metrics values. The targets are in the .npy files saved from
dap_metrics.py and plotting_efficiency.py

Fernando Oct. 23th 2024
"""
import numpy as np # type: ignore

# Load the targets.npy file
file_path = '/home/fercho/double-aperture-photometry/test_results/'
data_nominal_mask = np.load(file_path + 'targets_P5.npy')
data_extended_mask = np.load(file_path + 'targets_P5_extended.npy')
data_eta_nominal_mask = np.load(file_path + 'eta_bt_24_cameras.npy')
data_eta_extended_mask = np.load(file_path + 'eta_ext_bt_24_cameras.npy')

# Some statistical thresholds
flux_thrsh = 7.1
cob_thresh = 3
etx_flux_thrsh = 3

# Let's get all the target IDs and magnitudes
target_IDs = data_nominal_mask[:, 0]
target_magnitudes = data_nominal_mask[:, 1]
target_magnitudes_ext = data_extended_mask[:, 1]

# Let's get all the 10-element eta_k_ext values coming from plottiing_efficiency.py
eta_ext_flux = data_eta_extended_mask[:, :]

# Let's get all the 10-element eta_k_ext values coming from plottiing_efficiency.py
eta_nom_flux = data_eta_nominal_mask[:, :]

# Let's get all the 10-element eta_k_nom_cob values coming from dap_metrics.py
eta_nom_cob = data_nominal_mask[:, 46:56]

# Let's get the centroid shift noise over 1h for the nominal mask
centroid_shift_error_nom = data_nominal_mask[:, 56:66]

# Option to choose condition set: Set to 1 or 2 depending on what you want
# Set 1 of conditions (ETFX but no NCOB) 

# Now we have to find the index of targets for which:
# 1. eta_nom_flux > flux_thrsh (7.1)
# 2. eta_ext_flux > eta_cob_nom
# 3. eta_ext_flux > ext_flux_thrsh (3)
# 4. eta_nom_cob  < cob_thrsh (3)

# Set 2 for the new ones (NCOB but no ETFX)

# Now we have to find the index of targets for which:
# 1. eta_nom_flux > flux_thrsh (7.1)
# 2. eta_nom_cob > eta_ext_flux
# 3. eta_ext_flux < ext_flux_thresh (3)
# 4. eta_nom_cob > cob_thrsh (3)

# Set a condition
condition_set = 1  

if condition_set == 1:

    # Applying the conditions element-wise and finding the indices of the matching targets
    condition = (eta_nom_flux > flux_thrsh) & (eta_ext_flux > eta_nom_cob) & (eta_ext_flux > etx_flux_thrsh) & (eta_nom_cob < cob_thresh)

elif condition_set == 2:

    # Applying the conditions element-wise and finding the indices of the matching targets
    condition = (eta_nom_flux > flux_thrsh) & (eta_nom_cob > eta_ext_flux)  & (eta_ext_flux < etx_flux_thrsh) & (eta_nom_cob > cob_thresh)

# Now, find the target index and contaminant index where the condition holds true
index_of_matching_targets = []
index_of_matching_contaminants = []

# Loop through each target and contaminant
for target_index in range(condition.shape[0]):  # Loop over 7000 targets
    for contaminant_index in range(condition.shape[1]):  # Loop over 10 contaminant stars
        if condition[target_index, contaminant_index]:
            index_of_matching_targets.append(target_index)
            index_of_matching_contaminants.append(contaminant_index)

# Convert lists to arrays
index_of_matching_targets = np.array(index_of_matching_targets)
index_of_matching_contaminants = np.array(index_of_matching_contaminants)

# Print the results
print(f"Total matches: {len(index_of_matching_targets)}")
for i in range(len(index_of_matching_targets)):
    print(f"Target Index: {index_of_matching_targets[i]}, Contaminant Index: {index_of_matching_contaminants[i]}")

#Let's get a target that fulfills the conditons
good_target_ID = target_IDs[index_of_matching_targets[0]]
good_target_magnitude = target_magnitudes[index_of_matching_targets[0]]
good_target_centroid_shift_error_nom = centroid_shift_error_nom[index_of_matching_targets[0]]

# Let's get the NSR_1h for the extended mask
good_target_nsr_1h_extended_mask = data_extended_mask[index_of_matching_targets[0], 4]
good_target_spr_tot_ext = data_extended_mask[index_of_matching_targets[0], 13]
good_target_sprk_ext = data_extended_mask[index_of_matching_targets[0], 14:24]
good_target_eta_nom_cob_error = data_nominal_mask[index_of_matching_targets[0], 56:66]
good_target_eta_nom_cob_shifts = data_nominal_mask[index_of_matching_targets[0], 66:76]


print('Target ID:', good_target_ID)
print('Target NSR_ext_1h:', good_target_nsr_1h_extended_mask)
print('Target (1 - SPRtot_ext):', 1 - good_target_spr_tot_ext)
print('Target SPRk_ext', good_target_sprk_ext)
print('Target magnitude:', good_target_magnitude)
print('Target cob errors:', good_target_eta_nom_cob_error)
print('Target cob shifts:', good_target_eta_nom_cob_shifts)
