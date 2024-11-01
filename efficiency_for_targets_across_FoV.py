"""
Script for computing the efficiency for detecting FPs
with each method
Fernando 28th October 2024
"""
import numpy as np              # type:ignore
import matplotlib.pyplot as plt # type:ignore

# We load the data file
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr_Noblesse_PSF/'
data_nominal_mask = np.load(dataDIR + 'targets_P5.npy')


# We load the x and y positions of the targets across the FoV
x_fp = data_nominal_mask[:, 219]
y_fp = data_nominal_mask[:, 220]

#  Plot the target stars
plt.plot(x_fp, y_fp, 'o')

# Set the x and y limits
plt.xlim([-4810, 4810])
plt.ylim([-4810, 4810])

# Add horizontal and vertical lines at x=0 and y=0 (central axes)
plt.axhline(0, color='black',linewidth=1)  # horizontal line
plt.axvline(0, color='black',linewidth=1)  # vertical line

# Show the plot
plt.show()

# Now we define the rings
star_distances = np.sqrt(x_fp**2 + y_fp**2)  # Radial distances from the center

# We now define a radius
R = 1000
N = 7

# Filter stars to only include those within the disk radius
stars_in_disk = star_distances <= R
star_distances = star_distances[stars_in_disk]

# Calculate the boundaries for each ring and select stars within them
for i in range(N):
    r_i = np.sqrt(i / N) * R
    r_next = np.sqrt((i + 1) / N) * R
