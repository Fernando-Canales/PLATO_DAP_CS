import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

# We load the npy files with all the metrics of the nominal and secondary masks
data = np.load('targets_P5.npy')
data_sec = np.load('targets_P5_contaminant.npy')
# We obtain the magnitude of all the targets
mag = data[:, 1]

# We obtain the number of contaminant stars that could create a false positive for each target
n_bad = data[:, 8]

# We create a useful mask for getting the magnitude range of only P5 sample (P = 10.66 - 12.66)
mask_p5 = (mag >= 10.66) & (mag <= 12.66)

# Now we apply the mask
mag_p5 = mag[mask_p5]

# Now we apply the mask again to estimate the number of contaminant stars that could create a false positive
n_bad_p5 = n_bad[mask_p5]

# Now we make percentage histogram as the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
plt.hist(n_bad_p5, bins=bins, weights=[1/len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8)
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

# Now we get the eta_t and eta_c only for the stars of the P5 sample
eta_t_P5 = eta_t[mask_p5]
eta_c_P5 = eta_c[mask_p5]

# Now we create a mask to estimate the number of false positives detected by Marchiori's secondary mask. This means, the
# the number of cases where eta_c > eta_t given that eta_t > 7.1
mask_marchiori = (eta_c_P5 > eta_t_P5) & (eta_t_P5 > 7.1)

# Now we create a mask to estimate the total number of false positives of the P5 sample. This means, the number of cases
# where eta_t > 7.1
mask_total_false = (eta_t_P5 > 7.1)

# Now we can estimate the number of false positives detected by Marchiori's secondary mask
false_positives_marchiori = mask_marchiori.sum()

# Now we can estimate the total the number of false positives
false_positives_p5 = mask_total_false.sum()

# Now we obtain the efficiency
#print(((eta_c_P5 > eta_t_P5) & (eta_t_P5 > 7.1)).sum() / (eta_t_P5 > 7.1).sum())
print(false_positives_marchiori / false_positives_p5)

# Now we can plot an histogram that contains our estimations of the total number of false positives and the number of
# false positives detected by Marchiori's method for given magnitude bins
plt.hist(mag_p5[mask_total_false], bins=[10, 11, 12, 13], edgecolor='black', label='Nominal Mask', alpha=0.5)
plt.hist(mag_p5[mask_marchiori], bins=[10, 11, 12, 13], edgecolor='black', label='Secondary Mask', alpha=0.5)
plt.legend()
plt.show()

# We load the npy file with all the metrics of the extended mask
data_ext = np.load('targets_P5_extended.npy')

# Now we compute the efficiency of the extended mask method. For doing that, we compute the number of cases where the
# extended mask is able to detect a false positive over the total number of false positives detected by the nominal
# mask.

#sprk_ext = data_ext[]
#spr_crit_ex =
#sprk =
#spr_crit =

#mask_extended =