import numpy as np # type: ignore

# Load the targets.npy file
file_path = '/home/fercho/double-aperture-photometry/test_results/'  # Replace with your actual path
data = np.load(file_path + 'eta_bt_24_cameras.npy')

# Find the indices of rows where the first value is greater than 7.1
row_indices = np.where(data[:, 0] > 7.1)[0]

# Specify the index you want to check
specific_index = 117

# Check if the specific index is in row_indices
if specific_index in row_indices:
    print(f"Index {specific_index} is in the row_indices where the first value is > 7.1.")
else:
    print(f"Index {specific_index} is NOT in the row_indices where the first value is > 7.1.")