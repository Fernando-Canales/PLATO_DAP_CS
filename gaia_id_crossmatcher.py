import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

# I want to get the 'local' ID of a target for which eta_nom is > 7.1
# Target ID to find (assuming it corresponds to the row index)
#target_id = 627990
#target_id = 2805234
# Load the numpy array from the .npy file
gaia_catalogue_file_path = '/home/fercho/double-aperture-photometry/catalogues_stars/SFP_DR3_20230101.npy'
simulation_results_path = '/home/fercho/double-aperture-photometry/test_results/'

data_gaia_catalogue = np.load(gaia_catalogue_file_path, allow_pickle=True)
data_nominal_mask = np.load(simulation_results_path+'targets_P5.npy')

#data_gaia_catalogue = np.load(gaia_catalogue_file_path, allow_pickle=True)
#data_nominal_mask = np.load(simulation_results_path+'targets_P5.npy')

target_ids = data_nominal_mask[:, 0]

# Target ID to find
#target_id = 1911666
#target_id = 1135002
target_id = 2427
index_target_id = np.where(target_ids == target_id)[0]

contaminant_ids = data_nominal_mask[index_target_id, 189:199].flatten().astype(int)  # Ensure 1-D and convert to int


# Function to retrieve and print information based on a given ID
def print_gaia_info_by_id(given_id, label="Target"):
    if 0 <= given_id < data_gaia_catalogue.shape[0]:
        entry = data_gaia_catalogue[given_id]
        
        # Extract GAIAID1 and GAIAID2 from the corresponding columns
        gaiaid1 = entry[8]
        gaiaid2 = entry[9]
        
        # Compute the full GAIA ID
        full_gaia_id = int(gaiaid1 * 2**32) + int(gaiaid2)
        # Extract RA and Dec from the corresponding columns
        ra = entry[0]
        dec = entry[1]
        
        print(f"Full GAIA ID for {label} {given_id}: {full_gaia_id}")
        print(f"RA for {label} {given_id}: {ra} deg")
        print(f"Dec for {label} {given_id}: {dec} deg")
        print("-" * 50)
    else:
        print(f"{label} ID {given_id} is out of range.")
        print("-" * 50)

# Print information for the main target
print_gaia_info_by_id(target_id, label="Target")

# Loop through all contaminant IDs and print their information
for i in contaminant_ids:
    print_gaia_info_by_id(i, label="Contaminant")

# Extract RA and Dec for the target
target_entry = data_gaia_catalogue[target_id]
target_ra = target_entry[0]
target_dec = target_entry[1]

# Extract RA and Dec for contaminants
contaminant_ras = []
contaminant_decs = []
for contaminant_id in contaminant_ids:
    entry = data_gaia_catalogue[contaminant_id]
    contaminant_ras.append(entry[0])
    contaminant_decs.append(entry[1])

# Plotting
plt.figure(figsize=(8, 6))
plt.scatter(contaminant_ras, contaminant_decs, label="Contaminants", color='orange')
plt.scatter(target_ra, target_dec, label="Target", color='blue', marker='x', s=100)
plt.xlabel("RA (degrees)")
plt.ylabel("Dec (degrees)")
plt.title("Target and Contaminants in Focal Plane")
plt.legend()
plt.grid()
plt.show()