import numpy as np

# Let's calculate  the critical SPR
def spr_crit(dback, nsr, td, ntr):
    term1 = 7.1 / dback
    term2 = min(nsr) / np.sqrt(td * ntr)
    return term1 * term2

# Let's calculate the NSR of the system
def NSRn(sb, sd, sq, ft, fc):
    n = np.sqrt(ft + fc + sb ** 2 + sd ** 2 + sq ** 2) / ft
    return n

# Let's define a function that is just a multiplicatin for obtaining the NSR_1h, z is the factor of 10^6 / (12 * sqrt(number of cameras))
def nsr1h(z, nsr):
    return z*nsr


# Let's go Aussie mode, keegan, yo!
def nsr_AGG(x, y, sb, sd, sq):
    n = []
    for i in range(1, len(x)+1):
        n.append(np.sqrt(np.sum(x[:i] + y[:i] + sb ** 2 + sd ** 2 + sq ** 2))/np.sum(x[:i]))
    return n

# Let's try the code!
if __name__ == '__main__':
    x = np.array([100, 99, 90, 85, 60, 50])
    y = np.array([9, 1, 5, 50, 12, 8])
    sb = 45
    sd = 50.2
    sq = 7.2
    nsr_agg = nsr_AGG(x, y, sb, sd, sq)
    print(nsr_agg)
