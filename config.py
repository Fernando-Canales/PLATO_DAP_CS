"""
Configuration file for DAP analysis.
Contains all paths, constants, and parameters organized in dataclasses.

Created by Fernando
21.08.2025
"""
from dataclasses import dataclass

@dataclass
class PathsConfig:
    """Files paths and directories"""
    catalogue_dir: str = '/home/fercho/double-aperture-photometry/catalogues_stars/'
    psf_file: str = '/home/fercho/double-aperture-photometry/plato_psfs/PSF_Focus_0mu_0.2pxdif.npz'
    output_dir: str = '/home/fercho/double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/1000_targets_per_magnitude_bin'
    
