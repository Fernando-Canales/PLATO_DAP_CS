import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

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
mask_p5 = (mag >= 10.66) & (mag <= 12.66)

# Now we apply the mask for getting the magnitude range
mag_p5 = mag[mask_p5]

# Now we apply the mask again to estimate the number of N_bad in the P5 sample magnitude range given the nominal mask
n_bad_p5 = n_bad[mask_p5]

# Now we make percentage histogram as the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()

# Now we compute the efficiency of Marchiori's method. For doing that, we compute the number of false positives detected
# by the secondary mask over the total number of false positives detected by the nominal mask. For doing so, we compute
# the number of cases where the statistical significance of the secondary mask (eta_c) is higher than the statistical
# significance of the nominal mask (eta_t) given that the statistical significance of the nominal mask is higher than
# 7.1 (eta_t > 7.1).

# We get the eta_t and eta_c
eta_t = data[:, 12]
eta_c = data_sec[:, 6]
nsr1h = data[:, 7]

# Now we get the eta_t and eta_c only for the stars of the P5 sample
eta_t_P5 = eta_t[mask_p5]
eta_c_P5 = eta_c[mask_p5]
nsr1h_p5 = nsr1h[mask_p5]

# Now we create a mask to estimate the number of false positives detected by Marchiori's secondary mask. This means, the
# number of cases where eta_c > eta_t given that eta_t > 7.1
mask_marchiori = (eta_c_P5 > eta_t_P5) & (eta_t_P5 > 7.1)

# Now we create a mask to estimate the total number of false positives of the P5 sample. This means, the number of cases
# where eta_t > 7.1
mask_total_false = (eta_t_P5 > 7.1)

# Now we can estimate the number of false positives detected by Marchiori's secondary mask
false_positives_marchiori = mask_marchiori.sum()

# Now we can estimate the total the number of false positives
false_positives_p5 = mask_total_false.sum()

# Now we obtain the efficiency of Marchiori's method
# print(((eta_c_P5 > eta_t_P5) & (eta_t_P5 > 7.1)).sum() / (eta_t_P5 > 7.1).sum())
print("The efficiency of Marchiori's method is:", "%.4f" % (false_positives_marchiori / false_positives_p5))

# Now we can plot a histogram that contains our estimations of the total number of false positives and the number of
# false positives detected by Marchiori's method for given magnitude bins
plt.hist(mag_p5[mask_total_false], bins=[10, 11, 12, 13], edgecolor='black', label='Nominal Mask', alpha=0.5)
plt.hist(mag_p5[mask_marchiori], bins=[10, 11, 12, 13], edgecolor='black', label='Secondary Mask', alpha=0.5)
plt.legend()
plt.show()

# Now we compute the efficiency of the extended mask method. For doing that, we compute the number of cases where the
# extended mask is able to detect a false positive over the total number of false positives detected by the nominal
# mask.

sprk_ext = data_ext[:, 5]
spr_crit_ex = data_ext[:, 6]
sprk = data[:, 10]
spr_crit = data[:, 9]

# The efficiency of the extended mask is this the ratio between the  number of contaminants stars for which
# (sprk_ext > spr_crit_ex) and (sprk > spr_crit) and (sprk_ext>sprk) and N_bad  (i.e. the number of contaminants for
# which (sprk > spr_crit)

# So, I have to count the  number of contaminants stars for which the condition above is verified and store this number

# I will use the mask again
sprk_ext_P5 = sprk_ext[mask_p5]
sprk_crit_ext_P5 = spr_crit_ex[mask_p5]
sprk_P5 = sprk[mask_p5]
spr_crit_P5 = spr_crit[mask_p5]

# I will obtain the efficiency
mask_extended = (sprk_ext_P5 > sprk_crit_ext_P5) & (sprk_P5 > spr_crit_P5) & (sprk_ext_P5 > sprk_P5)

mask_nominal = (sprk_P5 > spr_crit_P5)

# Now we compute the number of N_bad_ext
N_bad_extended = mask_extended.sum()

# Now we compute number of N_bad
N_bad = mask_nominal.sum()

# Now we obtain the efficiency of the extended mask
print("The efficiency of Extended mask's method is:", "%.4f" % (N_bad_extended / N_bad))

# Now we can plot a histogram that contains our estimations of the total number of false positives and the number of
# false positives detected by Extended mask's method for given magnitude bins
plt.hist(mag_p5[N_bad], bins=[10, 11, 12, 13], edgecolor='black', label='Nominal Mask', alpha=0.5)
plt.hist(mag_p5[N_bad_extended], bins=[10, 11, 12, 13], edgecolor='black', label='Secondary Mask', alpha=0.5)
plt.legend()
plt.show()

# Now bray's assumption
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]

n_bad_bray_p5 = n_bad_bray[mask_p5]
nsr1h_bray_p5 = nsr1h_bray[mask_p5]

# Now we make percentage histogram as the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8, label='Marchiori', alpha=0.5)
plt.hist(n_bad_bray_p5, bins=bins, weights=[1 / len(n_bad_bray_p5)] * len(n_bad_bray_p5), edgecolor='black', rwidth=0.8, label='Bray', alpha=0.5)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.legend()
plt.show()

plt.plot(nsr1h_p5)
plt.show()
