"""
Script for computing the efficiency for detecting FPs
with each method
Fernando 28th October 2024
"""
import numpy as np              # type:ignore
import matplotlib.pyplot as plt # type:ignore
from fitting_psf import from_pix_2_mm 

# Load the data file
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr_Noblesse_PSF/'
data_nominal_mask = np.load(dataDIR + 'targets_P5.npy')
fsize=10

# Load the x and y positions of the targets across the FoV
x_fp = data_nominal_mask[:, 219]
y_fp = data_nominal_mask[:, 220]

# Convert the values from pixels to mm
x_fp_mm, y_fp_mm = from_pix_2_mm(x_star=x_fp, y_star=y_fp)

# First Plot: Only the target stars in the focal plane
plt.figure(figsize=(6, 6), dpi=120)  # Square plot for focal plane targets
plt.plot(x_fp_mm, y_fp_mm, 'o', markersize=2, label='Target stars')
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("mm", fontsize=fsize)
plt.ylabel("mm", fontsize=fsize)
#plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# Second Plot: Current plot with concentric circles and labels
plt.figure(figsize=(6, 12), dpi=120)  # Taller figure for the circles plot
plt.plot(x_fp_mm, y_fp_mm, 'o', markersize=2, label='Target stars')
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)

# Define parameters for the rings
R = 86  # Total radius of the disk in mm
N = 7   # Number of circles

# Create six circles labeled from Circle=1 to Circle=6
for i in range(N):
    r_i = (i + 1) / N * R  # Linear spacing
    circle = plt.Circle((0, 0), r_i, color='orange', linestyle='-', fill=True, linewidth=4, alpha=0.1)
    plt.gca().add_artist(circle)

    # Highlighted circumference with a thicker line
    edge_circle = plt.Circle((0, 0), r_i, color='red', linestyle='--', fill=False, linewidth=1.5)
    plt.gca().add_artist(edge_circle)

    # Add bold text label for each circle
    plt.text(0, r_i, f"Circle = {i+1}", 
             color='red', fontsize=10, ha='center', va='center', fontweight='bold',
             bbox=dict(facecolor='white', edgecolor='red', boxstyle='round,pad=0.3'))

# Set plot limits to fit the circles and targets comfortably
plt.xlim(-R, R)
plt.ylim(-R, R)

# Add labels, grid, and title
plt.xlabel("mm", fontsize=fsize)
plt.ylabel("mm", fontsize=fsize)
#plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.show()