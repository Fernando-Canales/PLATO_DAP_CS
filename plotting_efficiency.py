import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from pylab import *

dataDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/test_results/'

# Parameters for the plots
Pmin = 8
Pmax = 13
binsize = 0.5
nP = int((Pmax - Pmin) / binsize + 1)
fsize = 14

# We load the npy files with all the metrics of the nominal and secondary and extended masks
#data_mag = np.load('SFP_DR3_20220831.npy')
data = np.load(dataDIR + 'targets_P5.npy')
data_sec = np.load(dataDIR + 'targets_P5_contaminant.npy')
data_ext = np.load(dataDIR + 'targets_P5_extended.npy')
data_bray = np.load(dataDIR + 'targets_P5_bray.npy')
mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])

plt.plot(mag_value, star_count, 'o')
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel('Numb. of stars', fontsize=fsize)
plt.show()

# We obtain the magnitude of all the targets and the magnitude of the most problematic contaminants'
#mag_gaia = data_mag[:, 2]
mag = data[:, 1]
mag_bad = data[:, 3]

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
plt.xlabel('Number of potential N_bad', fontsize=fsize)
plt.ylabel('Fraction of targets [%] ', fontsize=fsize)
plt.show()

# Now we get all the etas and delta_obs
eta_t = data[:, 12]
delta_obs_t = data[:, 13]
eta_c = data_sec[:, 6]
delta_obs_c = data_sec[:, 7]
eta_ext = data_ext[:, 7]
delta_obs_ext = data_ext[:, 8]
nsr1h = data[:, 7]
nsr1h_sec = data_sec[:, 4]
nsr1h_ext = data_ext[:, 4]
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]


# We also obtain the mask size of every mask
key_nom = data[:, 5]
key_sec = data_sec[:, 2]
key_ext = data_ext[:, 2]
size_nom = data[:, 6]
size_sec = data_sec[:, 3]
size_e = data_ext[:, 3]


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
plt.xlabel('Number of potential N_bad', fontsize=fsize)
plt.ylabel('Fraction of targets [%]', fontsize=fsize)
plt.legend()
plt.show()


# Now I will make a plot like Réza's to show the efficiency of both Marchiori and Extended masks method as a function
# of the magnitude of the target
figure(0)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    sec = ((eta_t > 7.1) & (delta_obs_c > delta_obs_t) & (eta_c > 7.1))[m].sum() / (eta_t > 7.1)[m].sum()
    ext = ((eta_t > 7.1) & (delta_obs_ext > delta_obs_t) & (eta_ext > 7.1))[m].sum() / (eta_t > 7.1)[m].sum()
    scatter(Pi, sec*100, color='red')
    scatter(Pi, ext*100, color='blue')

xlabel("P Magnitude", fontsize=fsize)
ylabel("Efficiency[%]", fontsize=fsize)

#plt.show()

"""
Now we obtain the NSR for both Bray et al 2 x 2 and our nominal masks as a function of the Target magnitude
"""
figure(1)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    nsr_nominal = np.median(nsr1h[m])
    nsr_bray = np.median(nsr1h_bray[m])
    scatter(Pi, nsr_nominal, color='black')
    scatter(Pi, nsr_bray, color='orange')

xlabel(" P Magnitude", fontsize=fsize)
ylabel(r"NSR[ppm $hr^{\frac{1}{2}}$]", fontsize=fsize)

#plt.show()

"""
Now we obtain several plots for showing the average size of the Nominal, Secondary and Extended Masks as a function of 
the target magnitude
"""
figure(2)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    size_nominal = np.mean(size_nom[m])
    size_ext = np.mean(size_e[m])
    scatter(Pi, size_nominal, color='black')
    scatter(Pi, size_ext, color='blue')

xlabel(" P Magnitude", fontsize=fsize)
ylabel(r"Average mask size", fontsize=fsize)

"""
Now we plot the size of every mask as a function of the target P magnitude
"""
figure(3)
clf()
plot(mag, size_nom, 'k+', label='nominal mask')
plot(mag, size_sec, 'r+', label='secondary mask')
plot(mag, size_e, 'b+', label='extended mask')
legend()
xlabel('P')
ylabel(r'Mask size')


#plt.show()
figure(4)
clf()
for i in range(5, 30):
    Pi = 5 + i * binsize
    m_bad = (mag_bad >= Pi - binsize/2.) & (mag_bad <= Pi + binsize/2.)
    size_secondary = np.mean(size_sec[m_bad])
    scatter(Pi, size_secondary, color='red')

xlabel(" P Magnitude of the Contaminants", fontsize=fsize)
ylabel(r"Average sec. mask size", fontsize=fsize)

#plt.show()

"""
Now we obtain the degeneracy of the masks. For doing so we just need to know the number of unique mask keys. Let's
begin to plot the cumulative or total number of unique shapes of the secondary mask needed for all the most 
problematic contaminant stars
"""
figure(5)
clf()
for i in range(5, 30):
    Pi = 5 + i * binsize
    m_bad = (mag_bad <= Pi + binsize/2.)
    key_secondary = len(np.unique(key_sec[m_bad]))
    scatter(Pi, key_secondary, color='red')

xlabel("P Magnitude of the Contaminant", fontsize=fsize)
ylabel("Cum. count of mask shapes", fontsize=fsize)

#plt.show()

"""
Now we plot the cumulative or total number of nominal mask shapes to address the total number of target stars 
"""
figure(6)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag <= Pi + binsize/2.)
    pix_nominal = len(np.unique(key_nom[m]))
    pix_sec = len(np.unique(key_sec[m]))
    pix_ext = len(np.unique(key_ext[m]))
    scatter(Pi, pix_nominal, color='black')
    scatter(Pi, pix_sec, color='red')
    scatter(Pi, pix_ext, color='blue')

xlabel('P Magnitude', fontsize=fsize)
ylabel('Cum. count of mask shapes', fontsize=fsize)

#plt.show()

"""
Now let's plot the efficiency of the C.O.B. shift measurements
"""
delta_cob = data[:, 14]
eta_cob = data[:, 15]
sigma_cob = data[:, 16]
delta_cob_sec = data_sec[:, 8]
eta_cob_sec = data_sec[:, 9]
sigma_cob_sec = data_sec[:, 10]
delta_cob_ext = data_ext[:, 9]
eta_cob_ext = data_ext[:, 10]
sigma_cob_ext = data_ext[:, 11]

figure(7)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = ((eta_cob > 3) & (eta_t > 7.1))[m].sum() / (eta_t > 7.1)[m].sum()
    eff_cob_sec = ((eta_cob_sec > 3) & (eta_t > 7.1))[m].sum() / (eta_t > 7.1)[m].sum()
    eff_cob_ext = ((eta_cob_ext > 3) & (eta_t > 7.1))[m].sum() / (eta_t > 7.1)[m].sum()
    scatter(Pi, eff_cob*100, color='black')
    scatter(Pi, eff_cob_sec*100, color='red')
    scatter(Pi, eff_cob_ext*100, color='blue')

xlabel('P Magnitude', fontsize=fsize)
ylabel('Efficiency[%]', fontsize=fsize)

#plt.show()

"""
Now we plot the NSR over 1h  for every mask as a function of the target P magnitude
"""
figure(8)
clf()
plot(mag, nsr1h, 'k+', label='nominal mask')
plot(mag, nsr1h_sec, 'r+', label='secondary mask')
plot(mag, nsr1h_ext, 'b+', label='extended mask')
semilogy()
legend()
xlabel('P magnitude')
ylabel(r'$NSR_{1hr}$')

"""
Now we plot the COB shift as a function of the target P magnitude
"""
figure(9)
clf()
plot(mag, delta_cob, 'k+', label='nominal mask')
plot(mag, delta_cob_sec, 'r+', label='secondary mask')
plot(mag, delta_cob_ext, 'b+', label='extended mask')
xlabel('P magnitude')
ylabel(r'$\Delta_{COB}$')
semilogy()
legend()

figure(10)
clf()
plot(mag, sigma_cob, 'k+', label='nominal mask')
plot(mag, sigma_cob_sec, 'r+', label='secondary mask')
plot(mag, sigma_cob_ext, 'b+', label='extended mask')
xlabel('P magnitude')
ylabel(r'$\sigma_{COB} [pix]$')
semilogy()
legend()

figure(11)
clf()
plot(mag, sigma_cob, 'k+', label='noise')
plot(mag, delta_cob, 'b.', label='signal')
xlabel('P magnitude')
ylabel(r'$\Delta_{COB}$')
title('Nominal Mask')
semilogy()
legend()

show()