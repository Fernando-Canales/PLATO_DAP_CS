from dataclasses import dataclass

class ImagetteAndPsfDecompositionParameters:
    # Parameters for the imagette and PSF decomposition
    size_imagette_x = 6
    size_imagette_y = 6
    psf_sub_resolution = 128
    psf_bspline_resolution = 20
    
@dataclass
class TargetMagnitudeBinningParams:
    # Parameters for the magnitude intervals
    n_targets_per_magnitude_bin: int
    Pmin: float | int = 10                               # minimum magnitude
    Pmax: float | int = 13                               # maximum magnitude
    binsize: float = 0.5                           # binsize around every magnitude value
    
    def nP(self) -> int:
        return int((self.Pmax - self.Pmin) / self.binsize + 1)   # number of bins
    
magnitude_binning_params = TargetMagnitudeBinningParams(n_targets_per_magnitude_bin=10000)
print(magnitude_binning_params.nP())

magnitude_binning_params.n_targets_per_magnitude_bin # n_tar