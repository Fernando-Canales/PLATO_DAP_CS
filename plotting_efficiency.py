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
mask_p5 = (mag >= 10) & (mag <= 13)

# Now we apply the mask for getting the magnitude range
mag_p5 = mag[mask_p5]

# Now we apply the mask again to estimate the number of N_bad in the P5 sample magnitude range given the nominal mask
n_bad_p5 = n_bad[mask_p5]

# Now we plot a percentage histogram like the one presented by Marchiori
bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
plt.hist(n_bad_p5, bins=bins, weights=[1 / len(n_bad_p5)] * len(n_bad_p5), edgecolor='black', rwidth=0.8)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.title('N_bad')
plt.xlabel('N')
plt.show()

"""

The next thing to do is to compute the efficiency of both Marchiori (Secondary Mask) and Samadi (Extended mask) methods. 
For doing that is enough to compute the ratio between the number of false positives detected by each mask over the total
number of false positives detected by the Nominal mask. We summarize all this in the following two paragraphs:

We will call N_sec to the number of false positives detected by the Secondary Mask.
We will call N_ext to the number of false positives detected by the Extended mask.
We will call N_bad to the number of false positives detected by the Nominal mask (the TOTAL number of false positives).

The efficiency of the Secondary mask is therefore N_c / N_bad.
The efficiency of the Extended mask is therefore N_ext / N_bad.

We need to write the correct expressions to obtain these three quantities.

So far in our code, we've been assuming (and that is quite a strong assumption!) all contaminant stars for every target 
to be eclipsing binaries, every single one of them with the same transit depth. Furthermore, all the signals in our code
come from the contaminant stars. In such a case, the total number of detectable false positives for the Nominal mask 
(i.e N_bad) is just the number of contaminant stars that could generate a detectable signal on the Nominal mask 
(i.e eta_t > 7.1).

The number of detectable false positives by the Secondary mask is given by the number of stars that fulfils the 
following three requirements simultaneously (in no particular order):
 
1.- the signal of the transit is detectable on the Secondary mask (i.e. eta_c > 7.1) 
 
2.- the transit depth of the signal on the contaminant's photometry is greater than the produced transit depth on the 
target photometry (i.e. delta_obs_c > delta_obs_t). It is important to recall here that the transit depth of the 
target's photometry is produced by the signal coming from the contaminant star.
 
3.- the signal of the transit is detectable on the Nominal mask (i.e. eta_t > 7.1)
 
The number of detectable false positives by the Extended mask is computed in a complete analogous way to what was just
described for the Secondary mask. It is enough to change the subscripts "c" for "ext". 

"""

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
extended = []
secondary = []
magnitude = [10.5, 11, 11.5, 12, 12.5, 13]
for i in range(7, 13):
    mask_mag = (mag >= i + 0.5)
    mag_plot = mag[mask_mag]
    fp_marchiori = ((eta_t[mask_mag] > 7.1) & (delta_obs_c[mask_mag] > delta_obs_t[mask_mag]) & (
                eta_c[mask_mag] > 7.1)).sum()
    fp_extended = ((eta_t[mask_mag] > 7.1) & (delta_obs_ext[mask_mag] > delta_obs_t[mask_mag]) & (
                eta_ext[mask_mag] > 7.1)).sum()
    fp = (eta_t[mask_mag] > 7.1).sum()
    e_marchiori = np.median(fp_marchiori / fp)
    e_extended = np.median(fp_extended / fp)
    secondary.append(e_marchiori)
    extended.append(e_extended)
    # magnitude.append(min(mag_plot))
    # magnitude.append(np.mean(mag_plot))

extended = np.array(extended)
secondary = np.array(secondary)
magnitude = np.array(magnitude)

plt.plot(magnitude, secondary * 100, 'o', label='Secondary Mask')
plt.plot(magnitude, extended * 100, 'o', label='Extended Mask')
plt.xlabel('P Magnitude')
plt.ylabel('Efficiency [%]')
plt.legend()
plt.show()
# Now bray's assumption
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]

n_bad_bray_p5 = n_bad_bray[mask_p5]
nsr1h_bray_p5 = nsr1h_bray[mask_p5]
