"""
Script for making the schematics for the nominal, secondary and extended masks

28th of August 2024
Fernando Canales
"""
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
fsize = 14
#DIRout = '/home/fercho/Documents/PhD/Research_Papers/double-aperture-photometry/plots/'
DIRout = '/home/fercho/double-aperture-photometry/plots_pdfs/'
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

# Function to create individual plots or a combined plot
def plot_masks(save_individual=False):
    # Create a single figure with three subplots if save_individual is False
    if not save_individual:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Plot the nominal mask
        plot_single_mask(axes[0], nominal_mask, center_x, center_y, contaminants)
        axes[0].set_title('Nominal Mask')

        # Plot the extended mask with nominal mask in dashed lines
        plot_single_mask(axes[1], extended_mask, center_x, center_y, contaminants)
        axes[1].set_title('Extended Mask')
        # Manually add dashed lines for the nominal mask
        plot_nominal_mask_contour(axes[1], nominal_mask) # type: ignore

        # Plot the secondary mask with nominal mask in ddshed lines
        plot_single_mask(axes[2], secondary_mask + nominal_mask, center_x, center_y, contaminants)
        axes[2].set_title('Secondary Mask')
        plot_nominal_mask_contour(axes[2], nominal_mask) # type: ignore

        plt.tight_layout()
        plt.savefig(DIRout +'combined_masks_schematic.png', dpi=300)
        plt.show()
    else:
        # Plot and save each mask individually
        fig, ax = plt.subplots(figsize=(5, 5))
        plot_single_mask(ax, nominal_mask, center_x, center_y, contaminants)
        plt.tight_layout()
        plt.savefig(DIRout +'nominal_mask_schematic.png', dpi=300)
        plt.show()

        fig, ax = plt.subplots(figsize=(5, 5))
        plot_single_mask(ax, extended_mask, center_x, center_y, contaminants)
        plot_nominal_mask_contour(ax=ax, nominal_mask=nominal_mask, plot_color='black')  # Add dashed lines
        plt.tight_layout()
        plt.savefig(DIRout+'extended_mask_schematic.png', dpi=300)
        plt.show()

        fig, ax = plt.subplots(figsize=(5, 5))
        plot_single_mask(ax, secondary_mask, center_x, center_y, contaminants)
        plot_nominal_mask_contour(ax=ax, nominal_mask=nominal_mask, plot_color='white')
        plt.tight_layout()
        plt.savefig(DIRout+'secondary_mask_schematic.png', dpi=300)
        plt.show()

# Function to plot individual masks
def plot_single_mask(ax, mask, center_x, center_y, contaminants):
    ax.imshow(mask, origin='lower', extent=(0, 6, 0, 6), cmap='plasma')
    ax.plot(center_x, center_y, 'o', color='green', markersize=10, markeredgecolor='black')  # Center point
    for (x, y) in contaminants:
        ax.plot(x, y, '^', color='cyan', markersize=10, markeredgecolor='black')  # Contaminants
    ax.grid(True)
    ax.set_xlabel('x [pixels]', fontsize=12)
    ax.set_ylabel('y [pixels]', fontsize=12)

# Function to add dashed lines for the nominal mask
def plot_nominal_mask_contour(ax, nominal_mask, plot_color):
    for j in range(6):
        for i in range(6):
            if nominal_mask[j, i] == 1:
                ax.plot([i, i+1], [j, j], color=plot_color, linestyle='--', lw=2)  # Bottom edge
                ax.plot([i, i+1], [j+1, j+1], color=plot_color, linestyle='--', lw=2)  # Top edge
                ax.plot([i, i], [j, j+1], color=plot_color, linestyle='--', lw=2)  # Left edge
                ax.plot([i+1, i+1], [j, j+1], color=plot_color, linestyle='--', lw=2)  # Right edge


# Call the function with desired option
plot_masks(save_individual=True)  # Set to False to save combined plot