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
flux_trsh = 7.1
cob_trsh = 3
ntr = 3  # number of observed transit events
td = 4  # transit event duration in hours

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
eta_c_2 = data_sec[:, 15]
delta_obs_c_2 = data_sec[:, 16]
eta_ext = data_ext[:, 7]
delta_obs_ext = data_ext[:, 8]
eta_ext_2 = data_ext[:, 15]
delta_obs_ext_2 = data_ext[:, 16]
eta_ext_3 = data_ext[:, 23]
delta_obs_ext_3 = data_ext[:, 24]
nsr1h = data[:, 7]
nsr1h_sec = data_sec[:, 4]
nsr1h_sec_2 = data_sec[:, 13]
nsr1h_ext = data_ext[:, 4]
nsr1h_ext_2 = data_ext[:, 14]
nsr1h_ext_3 = data_ext[:, 22]
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]


# We also obtain the shape and size of every mask
key_nom = data[:, 5]
key_sec = data_sec[:, 2]
key_sec_2 = data_sec[:, 11]
key_ext = data_ext[:, 2]
key_ext_2 = data_ext[:, 12]
key_ext_3 = data_ext[:, 20]
size_nom = data[:, 6]
size_sec = data_sec[:, 3]
size_sec_2 = data_sec[:, 12]
size_e = data_ext[:, 3]
size_e_2 = data_ext[:, 13]
size_e_3 = data_ext[:, 21]


# Nozw all the metrics related to the COB
delta_cob = data[:, 14]
eta_cob = data[:, 15]
sigma_cob = data[:, 16]
delta_cob_sec = data_sec[:, 8]
eta_cob_sec = data_sec[:, 9]
sigma_cob_sec = data_sec[:, 10]
delta_cob_sec_2 = data_sec[:, 17]
eta_cob_sec_2 = data_sec[:, 18]
sigma_cob_sec_2 = data_sec[:, 19]
delta_cob_ext = data_ext[:, 9]
eta_cob_ext = data_ext[:, 10]
sigma_cob_ext = data_ext[:, 11]
delta_cob_ext_2 = data_ext[:, 16]
eta_cob_ext_2 = data_ext[:, 18]
sigma_cob_ext_2 = data_ext[:, 19]
delta_cob_ext_3 = data_ext[:, 25]
eta_cob_ext_3 = data_ext[:, 26]
sigma_cob_ext_3 = data_ext[:, 27]

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


"""
Now I will make a plot like Réza's to show the efficiency of both Marchiori and Extended masks method as a function
of the magnitude of the target
"""
# First the expressions for plotting the efficiency for detecting false positives for the different masks
fp = (eta_t > flux_trsh)  # false positive
sdr_flux = fp & (eta_c > flux_trsh) & (delta_obs_c > delta_obs_t)  # secondary mask false positive detection rate
s_2dr_flux = fp & (eta_c_2 > flux_trsh) & (delta_obs_c_2 > delta_obs_t)  # secondary mask false positive detection rate
edr_flux = fp & (eta_ext > flux_trsh) & (delta_obs_ext > delta_obs_t)  # extended mask false positive detection rate
e_2dr_flux = fp & (eta_ext_2 > flux_trsh) & (delta_obs_ext_2 > delta_obs_t)  # extended (2) mask false positive detection rate
e_3dr_flux = fp & (eta_ext_3 > flux_trsh) & (delta_obs_ext_3 > delta_obs_t)  # extended (3) mask false positive detection rate

figure(0)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    #scatter(Pi, sec, color='red')
    scatter(Pi, sec_2, color='green')
    #scatter(Pi, ext, color='blue')
    scatter(Pi, ext_2, color='cyan')
    #scatter(Pi, ext_3, color='magenta')

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
ylabel(r"$NSR_{1hr}[ppm \sqrt{hr}]$", fontsize=fsize)

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
    size_ext_2 = np.mean(size_e_2[m])
    size_ext_3 = np.mean(size_e_3[m])
    scatter(Pi, size_nominal, color='black')
    scatter(Pi, size_ext, color='blue')
    scatter(Pi, size_ext_2, color='cyan')
    scatter(Pi, size_ext_3, color='magenta')

xlabel(" P Magnitude", fontsize=fsize)
ylabel(r"Average mask size", fontsize=fsize)

"""
Now we plot the size of every mask as a function of the target P magnitude
"""
figure(3)
clf()
#plot(mag, size_nom, 'k+', label='Nominal Mask')
plot(mag, size_sec, 'r+', label='Secondary Mask')
plot(mag, size_sec_2, 'g+', label='Secondary Mask (1 pixel ring)')
#plot(mag, size_e, 'b+', label='Extended mask')
#plot(mag, size_e_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, size_e_3, 'm+', label='extended mask (3)')
legend()
xlabel('P Magnitude', fontsize=fsize)
ylabel(r'Mask size', fontsize=fsize)


#plt.show()
figure(4)
clf()
for i in range(5, 30):
    Pi = 5 + i * binsize
    m_bad = (mag_bad >= Pi - binsize/2.) & (mag_bad <= Pi + binsize/2.)
    size_secondary = np.mean(size_sec[m_bad])
    size_secondary_2 = np.mean(size_sec_2[m_bad])
    scatter(Pi, size_secondary, color='red')
    scatter(Pi, size_secondary_2, color='green')

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
    key_secondary_2 = len(np.unique(key_sec_2[m_bad]))
    scatter(Pi, key_secondary, color='red')
    scatter(Pi, key_secondary_2, color='red')

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
    pix_sec_2 = len(np.unique(key_sec_2[m]))
    pix_ext = len(np.unique(key_ext[m]))
    pix_ext_2 = len(np.unique(key_ext_2[m]))
    pix_ext_3 = len(np.unique(key_ext_3[m]))
    scatter(Pi, pix_nominal, color='black')
    scatter(Pi, pix_sec, color='red')
    scatter(Pi, pix_sec_2, color='green')
    scatter(Pi, pix_ext, color='blue')
    scatter(Pi, pix_ext_2, color='cyan')
    scatter(Pi, pix_ext_3, color='cyan')

xlabel('P Magnitude', fontsize=fsize)
ylabel('Cum. count of mask shapes', fontsize=fsize)

#plt.show()

"""
Now let's plot the efficiency of the C.O.B. shift measurements
"""
# First the expressions for computing the efficiency for detecting false positives of the cob shift for different masks
ndr_cob = fp & (eta_cob > cob_trsh)  # nominal mask false positive detection rate via cob shift
sdr_cob = fp & (eta_cob_sec > cob_trsh)  # secondary mask false positive detection rate via cob shift
s_2dr_cob = fp & (eta_cob_sec_2 > cob_trsh)  # secondary mask false positive detection rate via cob shift
edr_cob = fp & (eta_cob_ext > cob_trsh)  # extended mask false positive detection rate via cob shift
e_2dr_cob = fp & (eta_cob_ext_2 > cob_trsh)  # extended mask (2) false positive detection rate via cob shift
e_3dr_cob = fp & (eta_cob_ext_3 > cob_trsh)  # extended mask (3) false positive detection rate via cob shift

figure(7)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    #scatter(Pi, eff_cob_sec, color='red')
    scatter(Pi, eff_cob_sec_2, color='green')
    #scatter(Pi, eff_cob, color='black')
    #scatter(Pi, eff_cob_ext, color='blue')
    scatter(Pi, eff_cob_ext_2, color='cyan')
    #scatter(Pi, eff_cob_ext_2, color='magenta')

#legend()
xlabel('P Magnitude', fontsize=fsize)
ylabel('Efficiency[%]', fontsize=fsize)

"""
Now we plot the NSR over 1h  for every mask as a function of the target P magnitude
"""
figure(8)
clf()
#plot(mag, nsr1h, 'k+', label='Nominal mask')
#plot(mag, nsr1h_sec, 'r+', label='Secondary Mask')
plot(mag, nsr1h_sec_2, 'g+', label='Secondary Mask (1 pixel ring)')
plot(mag, nsr1h_ext, 'b+', label='Extended Mask')
plot(mag, nsr1h_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, nsr1h_ext_3, 'm+', label='extended mask (3)')
semilogy()
legend()
xlabel('P magnitude', fontsize=fsize)
ylabel(r'$NSR_{1hr} [ppm \sqrt{hr}]$', fontsize=fsize)

"""
Now we plot the statistical significance for every mask as a function of the target P magnitude
"""
figure(9)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    plot(Pi, np.median(eta_t[m]), 'ko')
    plot(Pi, np.median(eta_c[m]), 'ro')
    # plot(mag, eta_c_2, 'g+', label='Secondary Mask (1 pixel ring)')
    plot(Pi, np.median(eta_ext[m]), 'bo')
    # plot(mag, eta_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
    # plot(mag, eta_ext_3, 'm+', label='extended mask (3)')
    legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'], loc='best')
#semilogy()
#legend()
xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\eta$', fontsize=fsize)

"""
Now we plot the COB shift as a function of the target P magnitude
"""
figure(10)
clf()
#plot(mag, delta_cob, 'k+', label='nominal mask')
#plot(mag, delta_cob_sec, 'r+', label='secondary mask')
plot(mag, delta_cob_sec_2, 'g+', label='Secondary Mask (2 pixels ring)')
#plot(mag, delta_cob_ext, 'b+', label='extended mask')
plot(mag, delta_cob_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, delta_cob_ext_3, 'm+', label='extended mask(3)')
xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\Delta_{COB}$', fontsize=fsize)
semilogy()
legend()

"""
Now we plot the COB shift error as a function of the target P magnitude
"""
figure(11)
clf()
plot(mag, sigma_cob, 'k+', label='nominal mask')
plot(mag, sigma_cob_sec, 'r+', label='secondary mask')
#plot(mag, sigma_cob_sec_2, 'g+', label='secondary mask(2)')
plot(mag, sigma_cob_ext, 'b+', label='extended mask')
#plot(mag, sigma_cob_ext_2, 'c+', label='extended mask(2)')
#plot(mag, sigma_cob_ext_3, 'm+', label='extended mask(3)')
xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\sigma_{COB} [pix]$', fontsize=fsize)
semilogy()
legend()

"""
Now we plot the compqrison between the flux and COB shift methods as a function of the target P magnitude
"""
figure(12)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    plot(Pi, ext, 'b+')
    #plot(Pi, ext_2, 'c+')
    #plot(Pi, ext_3, 'm+')
    plot(Pi, eff_cob_sec, 'r^')
    #plot(Pi, eff_cob_sec_2, 'g^')
    #plot(Pi, eff_cob*100, 'k')
    plot(Pi, eff_cob_ext, 'b^')
    #plot(Pi, eff_cob_ext_2, 'c^')
    #plot(Pi, eff_cob_ext_3, 'm^')
    legend(['Sec. Mask Flux', 'Ext. Mask Flux',
            'Sec. Mask COB shift', 'Ext. Mask COB shift'], loc='best')

xlabel('P Magnitude', fontsize=fsize)
ylabel('Efficiency[%]', fontsize=fsize)


"""
Now we plot the delta obs
"""
figure(13)
clf()
ratio_delta_obs = delta_obs_c / delta_obs_c_2
plot(mag, delta_obs_c, 'r+', label='Secondary Mask')
plot(mag, delta_obs_c_2, 'g+', label='Secondary Mask (1 pixel ring)')
plot(mag, ratio_delta_obs, 'y+', label='$\delta_{obs_{sec}} / \delta_{obs_{sec_{1}}}$')
#plot(mag, delta_obs_ext, 'b+', label='Extended Mask')
#plot(mag, delta_obs_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, eta_ext_3, 'm+', label='extended mask (3)')
semilogy()
legend()
xlabel('P magnitude', fontsize=fsize)
#ylabel(r'$\delta_{obs_{sec}} / \delta_{obs_{sec_{1}}}$', fontsize=fsize)
ylabel(r'$\delta_{obs}[ppm]$', fontsize=fsize)

"""
Now we plot the eta ratios
"""
nom_eta = (eta_t / eta_cob)
sec_eta = (eta_c / eta_cob_sec)
ext_eta = (eta_ext / eta_cob_ext)

figure(14)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    plot(Pi, np.median(nom_eta[m]), 'ko')
    plot(Pi, np.median(sec_eta[m]), 'ro')
    plot(Pi, np.median(ext_eta[m]), 'bo')
    legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'],
           loc='best')
    #semilogy()
    #legend()

xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\eta_{flux} / \eta_{cob}$', fontsize=fsize)

figure(15)
clf()
plot(mag, eta_t/eta_cob, 'ko')
plot(mag, eta_c/eta_cob_sec, 'ro')
plot(mag, eta_ext/eta_cob_ext, 'bo')
legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'],
       loc='best')
semilogy()
xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\eta_{flux} / \eta_{cob}$', fontsize=fsize)


"""
Now we plot the delta_obs_ratio wit the median values
"""
figure(16)
clf()
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    plot(Pi, np.median(delta_obs_c_2[m]), 'yo')
    #plot(Pi, np.median(sec_eta[m]), 'ro')
    #plot(Pi, np.median(ext_eta[m]), 'bo')
    #legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'], loc='best')
    #semilogy()
    #legend()

xlabel('P magnitude', fontsize=fsize)
ylabel(r'$\delta_{obs_{sec}} / \delta_{obs_{sec_{1 pixel ring}}}$', fontsize=fsize)

show()