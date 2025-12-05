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
    output_dir: str = '/home/fercho/double-aperture-photometry/simulation_results/Distribution_transit_depths_and_durations/1000_targets_per_magnitude_bin'
    plots_output_dir: str = '/home/fercho/double-aperture-photometry/plots_pdfs/all_contaminants_are_EBs/Distribution_transit_depth_and_durations/'
    data_dir_plots: str = output_dir
    output_dir_plots: str = '/home/fercho/double-aperture-photometry/plots_pdfs/all_contaminants_are_EBs/Distribution_transit_depth_and_durations/'

    # Catalogue filenames
    star_catalogue_file: str = 'SFP_DR3_20230101.npy'
    eb_catalogue_file: str = 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt'
    psf_file: str = '/home/fercho/double-aperture-photometry/plato_psfs/PSF_Focus_0mu_0.2pxdif.npz' 

@dataclass
class InstrumentConfig:
    """Instrument and imagette parameters."""
    # Imagette size
    size_im_x: int = 6
    size_im_y: int = 6
    
    # PSF parameters
    subres: int = 128
    bsres: int = 20
    
    # Noise parameters (in e-/px)
    sb: float = 45.0 * 21  # Background noise * integration time (21 sec)
    sd: float = 50.2       # Detector noise (readout, smearing, dark current)
    sq: float = 7.2        # Quantization noise


@dataclass
class EclipsingBinaryConfig:
    """Eclipsing binary simulation parameters."""
    # Transit parameters
    transit_depth: float = 132000              # ppm
    transit_duration: float = 6.72 * 0.46**2   # hours
    ntr: int = 3                               # number of transits per hour
    
    # Contaminant selection
    distance_max: float = 7.0                  # max distance in pixels
    delta_p_max: float = 15.0                  # max magnitude difference
    max_number_of_contaminants: int = 500
    
    # EB occurrence rate settings
    eb_occurrence_rate: float = 0.01           # 1% based on Prša et al.
    use_realistic_eb_rate: bool = False        # Set to True to use realistic EB rate
    use_fixed_values: bool = False             # Set to True for fixed transit params


@dataclass
class MagnitudeConfig:
    """Magnitude binning configuration."""
    n_targets_per_bin: int = 1000
    pmin: float = 10.0
    pmax: float = 13.0
    binsize: float = 0.5
    
    @property
    def n_bins(self) -> int:
        """Number of magnitude bins."""
        return int((self.pmax - self.pmin) / self.binsize + 1)


@dataclass
class ThresholdConfig:
    """Detection thresholds and significance levels."""
    # Flux thresholds
    flux_thresh_nom_mask: float = 7.1
    flux_thresh_ext_mask: float = 3.0
    flux_thresh_sec_mask: float = 3.0
    
    # Centroid shift threshold
    cob_thresh: float = 3.0
    
    # Significance scaling
    depth_sig_scaling: float = 3.0
    gamma_factor_significance: float = 1.0    ## other 0.46


@dataclass
class PlottingConfig:
    """Plotting parameters and styling."""
    # Font and styling
    font_size: int = 14
    
    # Reference values
    td_ref: float = 6.72 * 0.46**2
    dback_ref: float = 132000
    
    # Plot limits
    ylim_min: float = 40.0
    ylim_max: float = 100.0
    xlim_min: float = 9.9
    xlim_max: float = 13.1
    
    # Region boundaries
    earth_like_limit: float = 11.7
    onboard_processing_limit: float = 10.7

@dataclass 
class RandomSeedConfig:
    """Random seed configuration for reproducibility."""
    main_seed: int = 300
    target_processing_seed: int = 123434434