import numpy as np

# Let's calculate  the critical SPR
def spr_crit(dback, nsr, td, ntr):
    term1 = 7.1 / dback
    term2 = min(nsr) / np.sqrt(td * ntr)
    return term1 * term2

# Let's calculate the NSR of the target
def NSRn(sb, sd, sq, ft, fc):
    n = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft
    return n

# Let's go Aussie mode, keegan, yo!
def nsr_AGG(x, y, sb, sd, sq):
    n = []
    for i in range(1, len(x)+1):
        n.append(np.sqrt(np.sum(x[:i] + y[:i] + sb + sd ** 2 + sq ** 2))/np.sum(x[:i]))
    return n

def aperture(ft, fc, sb, sd, sq):

    # First we compute the NSR of the system
    nsr = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft

    # Then we flatten that nsr and the fluxes
    nsr_1d = nsr.flatten()
    ft_1d = ft.flatten()
    fc_1d = fc.flatten()

    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting
    nsr_1d_index = np.argsort(nsr_1d)

    # Then we compute the intensity over those index for the target and the contaminant
    ft_1d = ft_1d[nsr_1d_index]
    fc_1d = fc_1d[nsr_1d_index]

    # Then we compute the aggregate noise-to-signal ratio
    n = []
    for i in range(1, len(ft_1d) + 1):
        n.append(np.sqrt(np.sum(ft_1d[:i] + fc_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i]))

    n = np.array(n)

    # We compute the nsr_agg over 1h
    nsr1h = ((10 ** 6) / (12 * np.sqrt(24))) * n

    # First we create a vector wiht only zeroes
    w = np.zeros(36)

    # Then we create a vector with just the amount of indexes of the mask
    mask = nsr_1d_index[:np.argmin(nsr1h) + 1]

    # Then we create our mask, we show the index where the mask vector has a value of 1
    w[mask] = 1

    return nsr1h, w

# Let's try the code!
if __name__ == '__main__':
    x = np.array([100, 99, 90, 85, 60, 50])
    y = np.array([9, 1, 5, 50, 12, 8])
    sb = 45
    sd = 50.2
    sq = 7.2
    nsr_agg = nsr_AGG(x, y, sb, sd, sq)
    print(nsr_agg)
