import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

# Parameters for the plots
Pmin = 8
Pmax = 14
binsize = 0.5
ntr = 3

# We load the npy files with all the metrics of the nominal and secondary and extended masks
data = np.load('targets_P5.npy')
data_sec = np.load('targets_P5_contaminant.npy')
data_ext = np.load('targets_P5_extended.npy')
data_bray = np.load('targets_P5_bray.npy')

# We obtain the magnitude of all the targets
mag = data[:, 1]

# We obtain the number of contaminant stars that could create a false positive for each target (i.e. N_bad)
n_bad = data[:, 8]

# We create a useful mask for getting the magnitude range for the P5 sample only (P = 10.66 - 12.66)
mask_p5 = (mag >= 10) & (mag <= 13)

# Now we apply the mask for getting the magnitude range
mag_p5 = mag[mask_p5]

# Now we apply the mask again to estimate the number of N_bad in the P5 sample magnitude range given the nominal mask
n_bad_p5 = n_bad[mask_p5]

# Now we plot a percentage histogram like the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Number of potential N_bad')
plt.ylabel('Fraction of targets [%] ')
plt.show()


# Now we get all the etas and delta_obs
eta_t = data[:, 12]
delta_obs_t = data[:, 13]
eta_c = data_sec[:, 6]
delta_obs_c = data_sec[:, 7]
eta_ext = data_ext[:, 7]
delta_obs_ext = data_ext[:, 8]
nsr1h = data[:, 7]

# Now we get the etas and delta_obs for the P5 sample
eta_t_P5 = eta_t[mask_p5]
delta_obs_t_P5 = delta_obs_t[mask_p5]
eta_c_P5 = eta_c[mask_p5]
delta_obs_c_P5 = delta_obs_c[mask_p5]
eta_ext_P5 = eta_ext[mask_p5]
delta_obs_ext_P5 = delta_obs_ext[mask_p5]
nsr1h_p5 = nsr1h[mask_p5]

# Now we estimate the number of false positives detected by the secondary mask (Marchiori's approach)
false_positives_marchiori = ((eta_t_P5 > 7.1) & (delta_obs_c_P5 > delta_obs_t_P5) & (eta_c_P5 > 7.1)).sum()

# Now we can estimate the total the number of false positives
false_positives_p5 = (eta_t_P5 > 7.1).sum()

# Now we obtain the efficiency of Marchiori's method
print("The efficiency of Marchiori's method is:", "%.4f" % (false_positives_marchiori / false_positives_p5))

# Now we estimate the number of false positives detected by the Extended mask
false_positives_extended = ((eta_t_P5 > 7.1) & (delta_obs_ext_P5 > delta_obs_t_P5) & (eta_ext_P5 > 7.1)).sum()

# Now we obtain the efficiency of the extended mask
print("The efficiency of Extended mask's method is:", "%.4f" % (false_positives_extended / false_positives_p5))

# Now I will make a plot like Réza's to show the efficiency of both Marchiori and Extended masks method as a function
# of the magnitude of the target

"""
extended = []
secondary = []
magnitude = [10.5, 11, 11.5, 12, 12.5, 13]
#magnitude = []
pi = []
a = np.arange(5, 11)
for i in range(a):
    # First we define a mask
    Pi = Pmin + i*binsize
    print(i)
    print(Pi)
    pi.append(Pi)
    #m = (mag >= Pi - binsize/2.) & (mag < Pi + binsize/2.)
    m = (mag >= 10.5 - binsize/2.) & (mag < 10.5 + binsize/2.)
    #mag_eff = np.median(mag[m])
    fp_marchiori = ((eta_t[m] > 7.1) & (delta_obs_c[m] > delta_obs_t[m]) & (
                eta_c[m] > 7.1)).sum()
    fp_extended = ((eta_t[m] > 7.1) & (delta_obs_ext[m] > delta_obs_t[m]) & (
                eta_ext[m] > 7.1)).sum()
    fp = (eta_t[m] > 7.1).sum()
    e_marchiori = np.median(fp_marchiori / fp)
    e_extended = np.median(fp_extended / fp)
    secondary.append(e_marchiori)
    extended.append(e_extended)
    #magnitude.append(mag_eff)

extended = np.array(extended)
secondary = np.array(secondary)
magnitude = np.array(magnitude)
pi = np.array(pi)
plt.plot(magnitude, secondary * 100, 'o-', label='Secondary Mask')
plt.plot(magnitude, extended * 100, 'o-', label='Extended Mask')
plt.xlabel('P Magnitude')
plt.ylabel('Efficiency [%]')
plt.legend()
plt.show()
"""

"""
extended = []
secondary = []
#magnitude = [10.5, 11, 11.5, 12, 12.5, 13]
magnitude = []
for i in range(5, 11):
    mask_mag = (mag > (8 + i*binsize) - 0.25) & (mag < (8 + i*binsize) + 0.25)
    magnitude_eff = np.median(mag[mask_mag])
    fp_marchiori = ((eta_t[mask_mag] > 7.1) & (delta_obs_c[mask_mag] > delta_obs_t[mask_mag]) & (
                eta_c[mask_mag] > 7.1)).sum()
    fp_extended = ((eta_t[mask_mag] > 7.1) & (delta_obs_ext[mask_mag] > delta_obs_t[mask_mag]) & (
                eta_ext[mask_mag] > 7.1)).sum()
    fp = (eta_t[mask_mag] > 7.1).sum()
    e_marchiori = np.median(fp_marchiori / fp)
    e_extended = np.median(fp_extended / fp)
    secondary.append(e_marchiori)
    extended.append(e_extended)
    magnitude.append(magnitude_eff)

extended = np.array(extended)
secondary = np.array(secondary)
magnitude = np.array(magnitude)

plt.plot(magnitude, secondary * 100, 'o-', label='Secondary Mask')
plt.plot(magnitude, extended * 100, 'o-', label='Extended Mask')
plt.xlabel('P Magnitude')
plt.ylabel('Efficiency [%]')
plt.legend()
plt.show()
"""


# Now bray's assumption
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]
"""
Now we will obtain a plot showing the amount of false positives detected by Bray et al 2 x 2 mask in comparison with the
amount of false positives detected by our Nominal mask
"""

n_bad_bray_p5 = n_bad_bray[mask_p5]

# Now we plot a percentage histogram like the one presented by Marchiori
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8,
         label='Nominal Mask', alpha=0.5)
plt.hist(n_bad_bray_p5, bins=bins, weights=[1 / len(n_bad_bray_p5)] * len(n_bad_bray_p5), edgecolor='black', rwidth=0.8,
         label='Bray 2 x 2 Mask', alpha=0.5)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Number of potential N_bad')
plt.ylabel('Fraction of targets [%]')
plt.legend()
plt.show()

#nsr1h_bray_p5 = nsr1h_bray[mask_p5]

"""
Now we obtain the NSR for both Bray et al 2 x 2 and our nominal masks as a function of the Target magnitude
"""

for i in range(8, 14):
    mask_mag = (mag >= i - 0.25) & (mag <= i + 0.25)
    magnitude_bray = np.median(mag[mask_mag])
    nsr_nominal = np.median(nsr1h[mask_mag])
    nsr_bray = np.median(nsr1h_bray[mask_mag])
    plt.scatter([magnitude_bray], [nsr_nominal], color='b')
    plt.scatter([magnitude_bray], [nsr_bray], color='orange')

plt.show()

"""
Now we obtain several plots for showing the degeneracy of the Nominal, Secondary and Extended Masks
"""

# First we obtain the mask size of every mask
key_nom = data[:, 5]
key_sec = data[:, 2]
key_ext = data[:, 2]
size_nom = data[:, 6]
size_sec = data_sec[:, 3]
size_e = data_ext[:, 3]


for i in range(8, 14):
    mask_mag = (mag >= i - 0.25) & (mag <= i + 0.25)
    magnitude_size = np.median(mag[mask_mag])
    size_nominal = np.mean(size_nom[mask_mag])
    size_secondary = np.mean(size_sec[mask_mag])
    size_ext = np.mean(size_e[mask_mag])
    plt.scatter([magnitude_size], [size_nominal], color='b')
    plt.scatter([magnitude_size], [size_ext], color='orange')


plt.show()

# Now we obtain the degeneracy of the masks. For doing so we just need to know the number of unique mask keys

for i in range(8, 14):
    mask_mag = (mag >= i - 0.25) & (mag <= i + 0.25)
    magnitude_size = np.median(mag[mask_mag])
    key_nominal = len(np.unique(key_nom[mask_mag]))
    key_secondary = len(np.unique(key_sec[mask_mag]))
    key_extended = len(np.unique(key_ext[mask_mag]))
    plt.scatter([magnitude_size], [key_nominal * 0.333], color='b')
    plt.scatter([magnitude_size], [key_extended * 0.333], color='orange')
    plt.xlabel(" P Magnitude")
    plt.ylabel(r"Fraction of unique mask shapes [%]")

plt.show()
