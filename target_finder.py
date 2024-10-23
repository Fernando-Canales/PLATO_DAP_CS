"""
Script that finds a target given specific conditions for the
metrics values. The targets are in the .npy files saved from
dap_metrics.py 

Fernando Oct. 23th 2024
"""
import numpy as np # type: ignore

# Load the targets.npy file
file_path = '/home/fercho/double-aperture-photometry/test_results/'
data_nominal_mask = np.load(file_path + 'targets_P5.npy')
data_eta_nominal_mask = np.load(file_path + 'eta_bt_24_cameras.npy')
data_eta_extended_mask = np.load(file_path + 'eta_ext_bt_24_cameras.npy')

# Some statistical thresholds
flux_thrsh = 7.1
cob_thresh = 3
etx_flux_thrsh = 3

# Let's get all the target IDs
target_IDs = data_nominal_mask[:, 0]

# Let's get all the 10-element eta_k_ext values coming from plottiing_efficiency.py
eta_ext_flux = data_eta_extended_mask[:, :]

# Let's get all the 10-element eta_k_ext values coming from plottiing_efficiency.py
eta_nom_flux = data_eta_nominal_mask[:, :]

# Let's get all the 10-element eta_k_nom_cob values coming from dap_metrics.py
eta_nom_cob = data_nominal_mask[:, 46:56]

print(target_IDs.shape)
print(eta_ext_flux.shape)
print(eta_nom_flux.shape)
print(eta_nom_cob.shape)

# Now we have to find the index of targets for which:
# 1. eta_nom_flux > flux_thrsh (7.1)
# 2. eta_ext_flux > eta_nom_cob
# 3. eta_ext_flux > ext_flux_thrsh
# 4. eta_nom_cob > cob_thrsh

# Applying the conditions element-wise and finding the indices of the matching targets
condition = (eta_nom_flux > flux_thrsh) & (eta_ext_flux > eta_nom_cob) & (eta_nom_cob > cob_thresh) & (eta_ext_flux > etx_flux_thrsh) 

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
print(index_of_matching_targets[0])
good_target = target_IDs[index_of_matching_targets[0]]
print(good_target)
