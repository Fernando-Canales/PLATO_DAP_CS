import numpy as np
import sys


# Let's calculate  the critical SPR
def spr_crit(dback, SPR_tot, nsr, td, ntr):
    flux_trsh = 7.1
    sprcrit = flux_trsh * nsr * (1 - SPR_tot) / (np.sqrt(td * ntr) * dback)
    return sprcrit


# Let's calculate the NSR of the target
def NSRn(sb, sd, sq, ft, fc):
    n = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft
    return n


# Let's go Aussie mode, keegan, yo!
def nsr_AGG(x, y, sb, sd, sq):
    n = []
    for i in range(1, len(x) + 1):
        n.append(np.sqrt(np.sum(x[:i] + y[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(x[:i]))
    return n



def aperture_computation(ft, fc, sb, sd, sq):
    """ We are following the procedure described
    in subsection 4.6.3. of Marchiori paper
    (Binary mask)

    Args:
        ft (_float_): target flux
        fc (_float_): contaminant flux
        sb (_float_): background noise
        sd (_float_): detector noise
        sq (_float_): quantization noise
    """
    # First we compute the NSR of the system. Eq. (36) of Marchiori's paper
    nsr = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft
    
    # Then we flatten the nsr and also the fluxes
    nsr_1d = nsr.flatten()
    ft_1d = ft.flatten()
    fc_1d = fc.flatten()
    
    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting them
    nsr_1d_index = np.argsort(nsr_1d)
    
    # Then we obtain the target and contaminant flux for such indexes
    ft_1d = ft_1d[nsr_1d_index]
    fc_1d = fc_1d[nsr_1d_index]
    
    # Then we compute the aggregat noise-to-signal ratio. Eq. (37) in Marchiori's paper
    nsr_agg = np.zeros(len(ft_1d))
    for i in range(1, len(ft_1d) + 1):
        nsr_agg[i - 1] = np.sqrt(np.sum(ft_1d[:i] + fc_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i])
    
    # We create now a vector full with zeros
    aperture = np.zeros(36)
    
    # Then we create a boolean mask
    boolean_mask_for_the_aperture = nsr_1d_index[:np.argmin(nsr_agg) + 1]
    
    # Then we obtain the nominal mask
    aperture[boolean_mask_for_the_aperture] = 1
    
    # Then we reshape the aperture
    aperture = aperture.reshape((6,6))
    
    return aperture

def aperture(ft, fc, sb, sd, sq):
    # First we compute the NSR of the system
    nsr = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft

    # Then we flatten that nsr and the fluxes since it is easier this way
    nsr_1d = nsr.flatten()
    ft_1d = ft.flatten()
    fc_1d = fc.flatten()

    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting
    nsr_1d_index = np.argsort(nsr_1d)

    # Then we compute the intensity over those index for the target and the contaminant
    ft_1d = ft_1d[nsr_1d_index]
    fc_1d = fc_1d[nsr_1d_index]

    # Then we compute the aggregate noise-to-signal ratio
    n = np.zeros(len(ft_1d))
    for i in range(1, len(ft_1d) + 1):
        n[i - 1] = np.sqrt(np.sum(ft_1d[:i] + fc_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i])

    # We compute the aggregate noise-to-signal ratio over 1h and 24 cameras
    nsr1h_24 = ((10 ** 6) / (12 * np.sqrt(24))) * n
    # We compute the aggregate noise-to-signal ratio over 1h and 6 cameras
    nsr1h_6 = ((10 ** 6) / (12 * np.sqrt(6))) * n
    

    # First we create a vector with only zeroes
    w_24_cameras = np.zeros(36)
    w_6_cameras = np.zeros(36)

    # Then we create a vector with just the amount of indexes of the mask
    mask_24_cameras = nsr_1d_index[:np.argmin(nsr1h_24) + 1]
    mask_6_cameras = nsr_1d_index[: np.argmin(nsr1h_6) + 1]

    # Then we create our mask, we show the index where the mask vector has a value of 1
    w_24_cameras[mask_24_cameras] = 1
    w_6_cameras[mask_6_cameras] = 1

    # Then we reshape the mask
    w_24_cameras = w_24_cameras.reshape((6, 6))
    w_6_cameras = w_6_cameras.reshape((6,6))

    return min(nsr1h_24), min(nsr1h_6), w_24_cameras, w_6_cameras


# We define now a function that computes the value of the spr_k for every contaminant as well as the maximum value of
# sprk, SPR_tot and the total number of stars for which spr_k is above SPR_crit

# ---------------------- We introduce now Réza's correction to Marchiori's formulas --------------------
# def SPR(SPR_crit, n_c, f_contaminant, f_tot, w):
# Then we compute the sprk over the extended mask for all the contaminants for a given target
# sprk = np.zeros(n_c)
# for i in range(0, n_c):
# We compute the sprk of every contaminant
# sprk[i] = np.sum(f_contaminant[i] * w) / np.sum(f_tot * w)

# We compute SPR_tot now
# SPR_tot = np.sum(sprk)

# And now we store here the highest sprk value
# sprk_max = max(sprk)

# Now we obtain the index of all contaminants above SPR_crit for a given target
# j = np.where(sprk > SPR_crit)[0]
# Now we obtain the number of contaminants above SPR_crit for a given target
# n_bad = len(j)
# return sprk, sprk_max, SPR_tot, n_bad

def SPR(n_c, f_contaminant, f_tot, w):
    # First we create a numpy array to store the sprk of all contaminants for a given target
    sprk = np.zeros(n_c)
    # Then we start a for loop over all the contaminants for a given target
    for i in range(1, n_c + 1):
        # Then we compute the sprk of every contaminant for a given target
        sprk[i - 1] = np.sum(f_contaminant[i - 1] * w) / np.sum(f_tot * w)

    # Then we compute the total contribution of all the contaminants for a given target (SPR_tot)
    SPR_tot = np.sum(sprk)

    return sprk, SPR_tot


# We define now a function for creating a mask_key as performed by Emmanuel
def mask_to_bitmask(mask):
    # As we save data within a bit array of 64 cells, we are able to save a mask of maximum size 8*8
    if mask.shape[0] * mask.shape[1] > 64:
        print('ERROR: Mask size too big to be converted into a 64 bits unsigned integer')
        sys.exit(0)

    bitmask = np.uint64(0)
    flat_mask = mask.flatten()  # This line flattens the mask array
    for i in range(flat_mask.shape[0]):
        if (flat_mask[i]) == 1:
            bitmask += np.uint64(2 ** i)
    return bitmask


# Now we define a function for obtaining a mask from a mask key
def bitmask_to_mask(bitmask, mask_row_nb, mask_col_nb):
    # As we save data within a bit array of 64 celss, we are able to build a mask of maximum 8  * 8
    if mask_row_nb * mask_col_nb > 64:
        print("ERROR: Mask size too big to be converted from a 64 bits unsigned integer")
        sys.exit(0)

    flat_mask = np.zeros(mask_row_nb * mask_col_nb, dtype=np.bool_)
    for i in range(mask_row_nb * mask_col_nb):
        flat_mask[i] = 1 if bitmask & np.uint64(2 ** i) == 2 ** i else 0
    return flat_mask.reshape(mask_row_nb, mask_col_nb)


# We define here a function for obtaining an extended mask given a nominal-binary mask
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


# Let's try the code!
if __name__ == '__main__':
    x = np.array([100, 99, 90, 85, 60, 50])
    y = np.array([9, 1, 5, 50, 12, 8])
    sb = 45
    sd = 50.2
    sq = 7.2
    nsr_agg = nsr_AGG(x, y, sb, sd, sq)
    print(nsr_agg)
