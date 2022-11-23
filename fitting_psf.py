# This module will find the closest PSF to every random selected target from the GAIA Catalogue.

# Let's import the main libraries
import numpy as np

# I know that 1 pixel = 18 micro meters.

# From mm to pixels
def from_mm_2_pix(xpsf, ypsf):
    # First we transform from mm to micro meters
    xpsf_micro = xpsf * 1000
    ypsf_micro = ypsf * 1000
    
    # Now let's convert the micro meters to mm
    xpsf_pix = xpsf_micro / 18
    ypsf_pix = ypsf_micro / 18
    return xpsf_pix, ypsf_pix

# From pixels to mm
def from_pix_2_mm(x_star, y_star):
    # First we transform from pix to micro meters
    x_micro = x_star * 18
    y_micro = y_star * 18

    # Then we transform from micro meters to mm
    x_mm = x_micro / 1000
    y_mm = y_micro / 1000
    return x_mm, y_mm

# Now let's find the closest PSF to every target. This means that the chosen PSF is such that the distance between the
# PSF and a given target are minimal.

def closest_psf(x_star, y_star, xpsf_pix, ypsf_pix):
    # We start a list where the index of the closest PSF to each target are gonna be saved
    j = []
    # We start the loop over the randomly chosen targets. 
    for i in range(0, len(x_star)):
        # We calculate the square distance of every target with every PSF of the catalogue
        sd = (xpsf_pix - x_star[i]) ** 2 + (ypsf_pix - y_star[i]) ** 2
        # We save in the list the index of the PSF that has the smallest distance with respect with every target
        j.append(np.argmin(sd) + 1)
    return j

# Now let's define the function for obtaining the contaminants around every target
def contaminants(x_star, y_star, x_tar, y_tar):
    l = []
    for i in range(0, len(x_tar)):
        dist = np.sqrt((x_star - x_tar[i]) ** 2 + (y_star - y_tar[i]) ** 2)
        m = (dist > 0) & (dist <= 10)
        j = np.where(m)[0]
        l.append(j)
        return l

# Let's create a function for obtaining the reference flux for a target star in the PLATO band for a normal camera
def reference_flux_target(p):
    """
    :param p: Magnitude of the star in the PLATO band
    :param zp: Zero point of the PLATO band
    :param 21: Integration time (21 seconds)
    """
    fp = (10 ** (-0.4 * (p - 20.62))) * 21.
    # print('The reference flux for a 6000 K G0V star of magnitude P =', p ,'after the integration time, but without '
    # 'brightness attenuation is:', "%.2f" % fp,
    # 'e-')
    return fp


# Let's create a function for obtaining the reference flux for a contaminant star in the PLATO band for a normal camera
def reference_flux_contaminant(ft, mc, mt):
    """
    :param mt: Magnitude of the target star in the PLATO band
    :param mc: Magnitude of the contaminant star in the PLATO band
    :param 21: Integration time (21 seconds), but we don't need it, since we already took it into account in ft(fp)
    """
    fc = (ft * 10 ** (-0.4 * (mc - mt)))
    return fc

        

