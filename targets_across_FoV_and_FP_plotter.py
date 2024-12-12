"""
Script for computing the efficiency for detecting FPs
with each method
Fernando 28th October 2024
"""
import numpy as np              # type:ignore
import matplotlib.pyplot as plt # type:ignore
from matplotlib.colors import LogNorm #type: ignore
from fitting_psf import from_pix_2_mm 
import matplotlib.patheffects as PathEffects # type: ignore

# Load the data file
dataDIR = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_132000ppm_and_td_1_422_hr_Noblesse_PSF/'
data_nominal_mask = np.load(dataDIR + 'targets_P5.npy')
fsize=10
# Define parameters for the rings
R = 86  # Total radius of the disk in mm
N = 7   # Number of circles
# Parameters
n = data_nominal_mask.shape[0]  # Number of targets
ntr = 3  # Number of transits
gamma_factor_significance = 1
td_ref = 6.72 * 0.46**2
dback_ref = 132000

# Load the x and y positions of the targets across the FoV
x_fp = data_nominal_mask[:, 219]  # Adjust index if necessary
y_fp = data_nominal_mask[:, 220]  # Adjust index if necessary

# Convert the values from pixels to mm
x_fp_mm, y_fp_mm = from_pix_2_mm(x_star=x_fp, y_star=y_fp)

# Compute eta_nom_bt_24_cameras
eta_nom_bt_24_cameras = np.zeros((n, 10))

for i in range(n):
    dback = np.ones(10) * dback_ref
    td = np.ones(10) * td_ref
    # Ensure indices correspond to your data structure
    sprk_10first = data_nominal_mask[i, 17:27]   # Adjust indices if necessary
    nsr_1h_24_cameras_nominal_mask = data_nominal_mask[i, 7]  # Adjust index if necessary
    SPR_tot = data_nominal_mask[i, 11]  # Adjust index if necessary

    eta_nom_bt_24_cameras[i, :] = (
        gamma_factor_significance * dback * sprk_10first * np.sqrt(td * ntr)
        / (nsr_1h_24_cameras_nominal_mask * (1 - SPR_tot))
    )

# Identify targets where any eta_nom_bt_24_cameras > 7.1
eta_nom_above_threshold = eta_nom_bt_24_cameras > 7.1
targets_above_threshold = np.any(eta_nom_above_threshold, axis=1)
num_targets_above_threshold = np.sum(targets_above_threshold)
print(f"Number of targets where any eta_nom_bt_24_cameras > 7.1: {num_targets_above_threshold}")

# Get the x and y coordinates of all targets
x_all = x_fp_mm
y_all = y_fp_mm

# Get the x and y coordinates of the selected targets
x_selected = x_fp_mm[targets_above_threshold]
y_selected = y_fp_mm[targets_above_threshold]

# Get the SPR_tot values for all targets
spr_tot_all = data_nominal_mask[:, 11]  # Adjust index if necessary

# For selected targets, get SPR_tot values
spr_tot_selected = spr_tot_all[targets_above_threshold]

# Define the concentric circles parameters
R = 86  # Total radius of the disk in mm
N = 7   # Number of circles

# Function to add concentric circles to the plot
def add_concentric_circles(ax, R, N):
    for i in range(N):
        r_i = (i + 1) / N * R  # Linear spacing for circle radii
        circle = plt.Circle((0, 0), r_i, color='darkorange', linestyle='--', fill=False, linewidth=2, zorder=3)
        ax.add_artist(circle)

# Plotting all targets, colored by SPR_tot, with concentric circles
fig, ax = plt.subplots(figsize=(6, 6), dpi=120)

# Add concentric circles
add_concentric_circles(ax, R, N)

# Plot the targets
sc_all = ax.scatter(x_all, y_all, c=spr_tot_all, cmap='viridis', s=20, zorder=2, norm=LogNorm())

# Add horizontal and vertical axes
ax.axhline(0, color='black', linewidth=2)
ax.axvline(0, color='black', linewidth=2)

# Add Roman numerals for the quadrants with a white font and black outline
txt_I = plt.text(R/2, R/2, 'I', color='white', fontsize=24, fontweight='bold', 
                 ha='center', va='center')
txt_II = plt.text(-R/2, R/2, 'II', color='white', fontsize=24, fontweight='bold', 
                  ha='center', va='center')
txt_III = plt.text(-R/2, -R/2, 'III', color='white', fontsize=24, fontweight='bold', 
                   ha='center', va='center')
txt_IV = plt.text(R/2, -R/2, 'IV', color='white', fontsize=24, fontweight='bold', 
                  ha='center', va='center')

# Apply path effects to outline the text in black
for txt in [txt_I, txt_II, txt_III, txt_IV]:
    txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='black')])

ax.set_xlim(-(R + 5), R + 5)
ax.set_ylim(-(R + 5), R + 5)

ax.set_xlabel("X Position (mm)", fontsize=fsize)
ax.set_ylabel("Y Position (mm)", fontsize=fsize)
ax.set_aspect('equal', adjustable='box')

cbar = plt.colorbar(sc_all, ax=ax)
cbar.set_label(r'$\rm SPR_{tot}$')

plt.savefig("targets_SPR_tot.pdf", format='pdf', bbox_inches='tight')
plt.show()

# Plotting only targets where any eta_nom_bt_24_cameras > 7.1, colored by SPR_tot, with concentric circles
fig, ax = plt.subplots(figsize=(6, 6), dpi=120)

# Add concentric circles
add_concentric_circles(ax, R, N)

# Plot the selected targets
sc_selected = ax.scatter(x_selected, y_selected, c=spr_tot_selected, cmap='viridis', s=20, zorder=2, norm=LogNorm())

# Add horizontal and vertical axes
ax.axhline(0, color='black', linewidth=2)
ax.axvline(0, color='black', linewidth=2)

ax.set_xlim(-(R + 5), R + 5)
ax.set_ylim(-(R + 5), R + 5)

ax.set_xlabel("X Position (mm)", fontsize=fsize)
ax.set_ylabel("Y Position (mm)", fontsize=fsize)
ax.set_aspect('equal', adjustable='box')

cbar = plt.colorbar(sc_selected, ax=ax)
cbar.set_label(r'$\rm SPR_{tot}$')
plt.savefig("targets_mask_SPR_tot.pdf", format='pdf', bbox_inches='tight')

plt.show()