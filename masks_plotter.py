"""
Scritp for making the schematics for the nominal, secondary and extended masks

28th of August 2024
Fernando Canales
"""
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
# Define center and contaminants
center_x, center_y = 2.5, 2.5  # Center of the 6x6 grid
contaminants = [
    (4.5, 4.5),
    (1.5, 4.5),
    (0.25, 0.75),
    (1.0, 1.5),
    (2.5, 1.5)
]

# Define the nominal mask
nominal_mask = np.zeros(36)
indices_for_the_nominal_mask = [8, 13, 14, 15, 20]
nominal_mask[indices_for_the_nominal_mask] = 1
nominal_mask = np.reshape(nominal_mask, (6, 6))

# Define the secondary mask
secondary_mask = np.zeros(36)
indices_for_the_secondary_mask = [28]
secondary_mask[indices_for_the_secondary_mask] = 1
secondary_mask = np.reshape(secondary_mask, (6, 6))

# Define the extended mask function
def extended_binary_mask(mask, W):
    ny, nx = mask.shape
    maske = np.zeros((ny, nx))
    maske[:, :] = mask[:, :]
    for j in range(ny):
        for i in range(nx):
            if mask[j, i] > 1.0 - 1e-5:
                for k in range(-W, W + 1):
                    if (j + k >= 0) & (j + k < ny):
                        for m in range(-W, W + 1):
                            if (i + m >= 0) & (i + m < nx):
                                maske[j + k, i + m] = 1
    return maske

extended_mask = extended_binary_mask(mask=nominal_mask, W=1)


# Define the extended mask
def extended_binary_mask(mask, W):
    ny, nx = mask.shape
    maske = np.zeros((ny, nx))
    maske[:, :] = mask[:, :]
    for j in range(ny):
        for i in range(nx):
            if mask[j, i] > 1.0 - 1e-5:
                for k in range(-W, W + 1):
                    if (j + k >= 0) & (j + k < ny):
                        for m in range(-W, W + 1):
                            if (i + m >= 0) & (i + m < nx):
                                maske[j + k, i + m] = 1
    return maske

extended_mask = extended_binary_mask(mask=nominal_mask, W=1)

# Create a single figure with three subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot the nominal mask
axes[0].imshow(nominal_mask, origin='lower', extent=(0, 6, 0, 6), cmap='plasma')
axes[0].plot(center_x, center_y, 'o', color='green', markersize=10, markeredgecolor='black')  # Center point
for (x, y) in contaminants:
    axes[0].plot(x, y, '^', color='cyan', markersize=10, markeredgecolor='black')  # Contaminants
axes[0].grid(True)  # Add a simple grid
axes[0].set_title('Nominal Mask')

# Plot the secondary mask with nominal mask overlaid in a different color
axes[1].imshow(secondary_mask + nominal_mask, origin='lower', extent=(0, 6, 0, 6), cmap='plasma')
axes[1].plot(center_x, center_y, 'o', color='green', markersize=10, markeredgecolor='black')  # Center point
for (x, y) in contaminants:
    axes[1].plot(x, y, '^', color='cyan', markersize=10, markeredgecolor='black')  # Contaminants
axes[1].grid(True)  # Add a simple grid
axes[1].set_title('Secondary Mask  (and Nominal Mask)')

# Plot the extended mask with nominal mask in dashed lines
axes[2].imshow(extended_mask, origin='lower', extent=(0, 6, 0, 6), cmap='plasma')
axes[2].plot(center_x, center_y, 'o', color='green', markersize=10, markeredgecolor='black')  # Center point
for (x, y) in contaminants:
    axes[2].plot(x, y, '^', color='cyan', markersize=10, markeredgecolor='black')  # Contaminants


axes[2].grid(True)  # Add a simple grid
axes[2].set_title('Extended Mask with Dashed Nominal Mask')

# Adjust layout and save the figure as a PDF
plt.tight_layout()
plt.savefig('corrected_masks_schematic.pdf')
plt.show()