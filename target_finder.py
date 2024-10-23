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

# Now we have to find the index of a target for which eta_nom > 7.1 and eta_ext_flux > eta_nom_cob