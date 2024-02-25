import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

#dataDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/test_results/'
dataDIR = '/home/fercho/double-aperture-photometry/test_results/'
# Parameters for the plots
Pmin = 10
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
data_sec = np.load(dataDIR + 'targets_P5_secondary.npy')
data_ext = np.load(dataDIR + 'targets_P5_extended.npy')
data_bray = np.load(dataDIR + 'targets_P5_bray.npy')
mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])
dback = 85000  # transit depth in ppm
td = 4         # transit duration in hours
ntr = 3        # number of transits in one hour

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
#eta_c_2 = data_sec[:, 15]
#delta_obs_c_2 = data_sec[:, 16]
eta_ext = data_ext[:, 7]
eta_ext_correct = data_ext[:, 8]
delta_obs_ext = data_ext[:, 9]
#eta_ext_2 = data_ext[:, 15]
#delta_obs_ext_2 = data_ext[:, 16]
#eta_ext_3 = data_ext[:, 23]
#delta_obs_ext_3 = data_ext[:, 24]
nsr1h = data[:, 7]
nsr1h_sec = data_sec[:, 4]
spr_crit = data[:, 8]
#nsr1h_sec_2 = data_sec[:, 13]
nsr1h_ext = data_ext[:, 4]
#nsr1h_ext_2 = data_ext[:, 14]
#nsr1h_ext_3 = data_ext[:, 22]
n_bad_bray = data_bray[:, 4]
nsr1h_bray = data_bray[:, 3]



n = data.shape[0]

# We also obtain the sprk and sprcrit
spr_tot_ext = data_ext[:, 13]
sprk_ext_10 = data_ext[:, 14:23]
eta_ext_10 = data_ext[:, 24:33]



# We also obtain the shape and size of every mask
key_nom = data[:, 5]
key_sec = data_sec[:, 2]
#key_sec_2 = data_sec[:, 11]
key_ext = data_ext[:, 2]
#key_ext_2 = data_ext[:, 12]
#key_ext_3 = data_ext[:, 20]
size_nom = data[:, 6]
size_sec = data_sec[:, 3]
#size_sec_2 = data_sec[:, 12]
size_e = data_ext[:, 3]
#size_e_2 = data_ext[:, 13]
#size_e_3 = data_ext[:, 21]


# Now all the metrics related to the COB
delta_cob = data[:, 14]
eta_cob = data[:, 15]
sigma_cob = data[:, 16]
#eta_cob_wrong = data[:, 17]
#sigma_cob_wrong = data[:, 18]
delta_cob_sec = data_sec[:, 8]
eta_cob_sec = data_sec[:, 9]
sigma_cob_sec = data_sec[:, 10]
#eta_cob_sec_wrong = data_sec[:, 20]
#sigma_cob_sec_wrong = data_sec[:, 21]
#delta_cob_sec_2 = data_sec[:, 17]
#eta_cob_sec_2 = data_sec[:, 18]
#sigma_cob_sec_2 = data_sec[:, 19]
#eta_cob_sec_2_wrong = data_sec[:, 22]
#igma_cob_sec_2_wrong = data_sec[:, 23]
delta_cob_ext = data_ext[:, 10]
eta_cob_ext = data_ext[:, 11]
sigma_cob_ext = data_ext[:, 12]
#eta_cob_ext_wrong = data_ext[:, 31]
#sigma_cob_ext_wrong = data_ext[:, 32]
#delta_cob_ext_2 = data_ext[:, 16]
#eta_cob_ext_2 = data_ext[:, 18]
#sigma_cob_ext_2 = data_ext[:, 19]
#eta_cob_ext_2_wrong = data_ext[:, 33]
#sigma_cob_ext_2_wrong = data_ext[:, 34]
#delta_cob_ext_3 = data_ext[:, 25]
#eta_cob_ext_3 = data_ext[:, 26]
#sigma_cob_ext_3 = data_ext[:, 27]
#eta_cob_ext_3_wrong = data_ext[:, 35]
#sigma_cob_ext_3_wrong = data_ext[:, 36]
SPRK10_first = data_ext[:,14:24]
eta_ext_bt = np.zeros((n, 10))
eta_bt = np.zeros((n, 10))

for i in range(n):
    eta_bt[i, :] = (SPRK10_first[i, :]/data[i, 8])*flux_trsh*(dback)*np.sqrt(td)
    eta_ext_bt[i, :] = dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 4] * (1 - data_ext[i, 13])) 

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
#s_2dr_flux = fp & (eta_c_2 > flux_trsh) & (delta_obs_c_2 > delta_obs_t)  # secondary mask false positive detection rate
edr_flux = fp & (eta_ext > flux_trsh) & (delta_obs_ext > delta_obs_t)  # extended mask false positive detection rate
edr_flux_correct = fp & (eta_ext_correct > flux_trsh) & (delta_obs_ext > delta_obs_t)  # extended mask false positive detection rate
#e_2dr_flux = fp & (eta_ext_2 > flux_trsh) & (delta_obs_ext_2 > delta_obs_t)  # extended (2) mask false positive detection rate
#e_3dr_flux = fp & (eta_ext_3 > flux_trsh) & (delta_obs_ext_3 > delta_obs_t)  # extended (3) mask false positive detection rate
#edr_flux_ext_overall = fp & (eta_ext_bt>flux_trsh)

plt.figure(0)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_correct = (edr_flux_correct[m].sum() / fp[m].sum()) * 100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    #scatter(Pi, sec, color='red')
    #plt.scatter(Pi, sec_2, color='green')
    #scatter(Pi, ext, color='blue')
    #plt.scatter(Pi, ext_2, color='cyan')
    #scatter(Pi, ext_3, color='magenta')

plt.xlabel("P Magnitude", fontsize=fsize)
plt.ylabel("Efficiency[%]", fontsize=fsize)
plt.title("Extended Mask", fontsize=fsize)
#plt.show()

"""
Now we obtain the NSR for both Bray et al 2 x 2 and our nominal masks as a function of the Target magnitude
"""
plt.figure(1)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    nsr_nominal = np.median(nsr1h[m])
    nsr_bray = np.median(nsr1h_bray[m])
    plt.scatter(Pi, nsr_nominal, color='black')
    plt.scatter(Pi, nsr_bray, color='orange')

plt.xlabel(" P Magnitude", fontsize=fsize)
plt.ylabel(r"$NSR_{1hr}[ppm \sqrt{hr}]$", fontsize=fsize)
plt.title("NSR for Bray and Marchiori nominal masks")

"""
Now we obtain several plots for showing the average size of the Nominal, Secondary and Extended Masks as a function of 
the target magnitude
"""
plt.figure(2)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    size_nominal = np.mean(size_nom[m])
    size_ext = np.mean(size_e[m])
    #size_sec = np.mean(size_sec[m])
    #size_ext_2 = np.mean(size_e_2[m])
    #size_ext_3 = np.mean(size_e_3[m])
    plt.scatter(Pi, size_nominal, color='black')
    plt.scatter(Pi, size_ext, color='blue')
    #plt.scatter(Pi, size_sec, color='red')
    #plt.scatter(Pi, size_ext_2, color='cyan')
    #plt.scatter(Pi, size_ext_3, color='magenta')

plt.xlabel(" P Magnitude", fontsize=fsize)
plt.ylabel(r"Average mask size", fontsize=fsize)
plt.title("Average mask size")

"""
Now we plot the size of every mask as a function of the target P magnitude
"""
plt.figure(3)
plt.plot(mag, size_nom, 'k+', label='Nominal Mask')
plt.plot(mag, size_sec, 'r+', label='Secondary Mask')
#plt.plot(mag, size_sec_2, 'g+', label='Secondary Mask (1 pixel ring)')
plt.plot(mag, size_e, 'b+', label='Extended mask')
#plot(mag, size_e_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, size_e_3, 'm+', label='extended mask (3)')
plt.legend()
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel(r'Mask size', fontsize=fsize)

plt.figure(4)
for i in range(5, 30):
    Pi = 5 + i * binsize
    m_bad = (mag_bad >= Pi - binsize/2.) & (mag_bad <= Pi + binsize/2.)
    size_secondary = np.mean(size_sec[m_bad])
    #size_secondary_2 = np.mean(size_sec_2[m_bad])
    plt.scatter(Pi, size_secondary, color='red')
    #plt.scatter(Pi, size_secondary_2, color='green')

plt.xlabel(" P Magnitude of the Contaminants", fontsize=fsize)
plt.ylabel(r"Average sec. mask size", fontsize=fsize)

"""
Now we obtain the degeneracy of the masks. For doing so we just need to know the number of unique mask keys. Let's
begin to plot the cumulative or total number of unique shapes of the secondary mask needed for all the most 
problematic contaminant stars
"""
plt.figure(5)
for i in range(5, 30):
    Pi = 5 + i * binsize
    m_bad = (mag_bad <= Pi + binsize/2.)
    key_secondary = len(np.unique(key_sec[m_bad]))
    #key_secondary_2 = len(np.unique(key_sec_2[m_bad]))
    plt.scatter(Pi, key_secondary, color='red')
    #plt.scatter(Pi, key_secondary_2, color='red')

plt.xlabel("P Magnitude of the Contaminant", fontsize=fsize)
plt.ylabel("Cum. count of mask shapes", fontsize=fsize)

#plt.show()

"""
Now we plot the cumulative or total number of nominal mask shapes to address the total number of target stars 
"""
plt.figure(6)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag <= Pi + binsize/2.)
    pix_nominal = len(np.unique(key_nom[m]))
    pix_sec = len(np.unique(key_sec[m]))
    #pix_sec_2 = len(np.unique(key_sec_2[m]))
    pix_ext = len(np.unique(key_ext[m]))
    #pix_ext_2 = len(np.unique(key_ext_2[m]))
    #pix_ext_3 = len(np.unique(key_ext_3[m]))
    plt.scatter(Pi, pix_nominal, color='black')
    plt.scatter(Pi, pix_sec, color='red')
    #plt.scatter(Pi, pix_sec_2, color='green')
    plt.scatter(Pi, pix_ext, color='blue')
    #plt.scatter(Pi, pix_ext_2, color='cyan')
    #plt.scatter(Pi, pix_ext_3, color='cyan')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Cum. count of mask shapes', fontsize=fsize)

"""
Now let's plot the efficiency of the C.O.B. shift measurements
"""
# First the expressions for computing the efficiency for detecting false positives of the cob shift for different masks
ndr_cob = fp & (eta_cob > cob_trsh)  # nominal mask false positive detection rate via cob shift
sdr_cob = fp & (eta_cob_sec > cob_trsh)  # secondary mask false positive detection rate via cob shift
#s_2dr_cob = fp & (eta_cob_sec_2 > cob_trsh)  # secondary mask false positive detection rate via cob shift
edr_cob = fp & (eta_cob_ext > cob_trsh)  # extended mask false positive detection rate via cob shift
#e_2dr_cob = fp & (eta_cob_ext_2 > cob_trsh)  # extended mask (2) false positive detection rate via cob shift
#e_3dr_cob = fp & (eta_cob_ext_3 > cob_trsh)  # extended mask (3) false positive detection rate via cob shift

plt.figure(7)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    plt.scatter(Pi, eff_cob_sec, color='red')
    #plt.scatter(Pi, eff_cob_sec_2, color='green')
    plt.scatter(Pi, eff_cob, color='black')
    plt.scatter(Pi, eff_cob_ext, color='blue')
    #plt.scatter(Pi, eff_cob_ext_2, color='cyan')
    #scatter(Pi, eff_cob_ext_2, color='magenta')

#legend()
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('COB shift Efficiency for every mask', fontsize=fsize)

"""
Now we plot the NSR over 1h  for every mask as a function of the target P magnitude
"""
plt.figure(8)
plt.plot(mag, nsr1h, 'k+', label='Nominal mask')
plt.plot(mag, nsr1h_sec, 'r+', label='Secondary Mask')
#plt.plot(mag, nsr1h_sec_2, 'g+', label='Secondary Mask (1 pixel ring)')
plt.plot(mag, nsr1h_ext, 'b+', label='Extended Mask')
#plt.plot(mag, nsr1h_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, nsr1h_ext_3, 'm+', label='extended mask (3)')
plt.semilogy()
plt.legend()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$NSR_{1hr} [ppm \sqrt{hr}]$', fontsize=fsize)
plt.title(r'$NSR_{1hr}$ for every mask')

"""
Now we plot the statistical significance for every mask as a function of the target P magnitude
"""
plt.figure(9)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    plt.plot(Pi, np.median(eta_t[m]), 'ko')
    plt.plot(Pi, np.median(eta_c[m]), 'ro')
    # plot(mag, eta_c_2, 'g+', label='Secondary Mask (1 pixel ring)')
    plt.plot(Pi, np.median(eta_ext[m]), 'bo')
    plt.plot(Pi, np.median(eta_ext_correct[m]), 'go')
    # plot(mag, eta_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
    # plot(mag, eta_ext_3, 'm+', label='extended mask (3)')
    plt.legend(['Nominal Mask', 'Secondary Mask', 'Extended mask', 'Extended corr'], loc='best')

plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\eta$', fontsize=fsize)

"""
Now we plot the COB shift as a function of the target P magnitude
"""
plt.figure(10)
plt.plot(mag, delta_cob, 'k+', label='nominal mask')
plt.plot(mag, delta_cob_sec, 'r+', label='secondary mask')
#plt.plot(mag, delta_cob_sec_2, 'g+', label='Secondary Mask (2 pixels ring)')
plt.plot(mag, delta_cob_ext, 'b+', label='extended mask')
#plt.plot(mag, delta_cob_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, delta_cob_ext_3, 'm+', label='extended mask(3)')
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\Delta_{COB}$', fontsize=fsize)
plt.semilogy()
plt.legend()

"""
Now we plot the COB shift error as a function of the target P magnitude
"""
plt.figure(11)
plt.plot(mag, sigma_cob, 'k+', label='nominal mask')
plt.plot(mag, sigma_cob_sec, 'r+', label='secondary mask')
#plot(mag, sigma_cob_sec_2, 'g+', label='secondary mask(2)')
plt.plot(mag, sigma_cob_ext, 'b+', label='extended mask')
#plot(mag, sigma_cob_ext_2, 'c+', label='extended mask(2)')
#plot(mag, sigma_cob_ext_3, 'm+', label='extended mask(3)')
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma_{COB} [pix]$', fontsize=fsize)
plt.semilogy()
plt.legend()

"""
Now we plot the comparison between the flux and COB shift methods as a function of the target P magnitude
"""
plt.figure(12)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_corr = (edr_flux_correct[m].sum() / fp[m].sum()) * 100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    plt.plot(Pi, ext, 'b+')
    plt.plot(Pi, ext_corr, 'g+')
    #plot(Pi, ext_2, 'c+')
    #plot(Pi, ext_3, 'm+')
    plt.plot(Pi, eff_cob, 'k^')
    plt.plot(Pi, eff_cob_sec, 'r^')
    #plot(Pi, eff_cob_sec_2, 'g^')
    #plot(Pi, eff_cob*100, 'k')
    plt.plot(Pi, eff_cob_ext, 'b^')
    #plot(Pi, eff_cob_ext_2, 'c^')
    #plot(Pi, eff_cob_ext_3, 'm^')
    plt.legend(['Sec. Mask Flux', 'Ext. Mask Flux', 'Ext. corr. flux', 'Nom. Mask COB shift',
            'Sec. Mask COB shift', 'Ext. Mask COB shift'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison of Every Method for Every Mask', fontsize=fsize)


"""
Now we plot the delta obs
"""
plt.figure(13)
#ratio_delta_obs = delta_obs_c / delta_obs_c_2
plt.plot(mag, delta_obs_c, 'r+', label='Secondary Mask')
#plt.plot(mag, delta_obs_c_2, 'g+', label='Secondary Mask (1 pixel ring)')
#plt.plot(mag, ratio_delta_obs, 'y+', label='$\delta_{obs_{sec}} / \delta_{obs_{sec_{1}}}$')
#plot(mag, delta_obs_ext, 'b+', label='Extended Mask')
#plot(mag, delta_obs_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
#plot(mag, eta_ext_3, 'm+', label='extended mask (3)')
plt.semilogy()
plt.legend()
plt.xlabel('P magnitude', fontsize=fsize)
#ylabel(r'$\delta_{obs_{sec}} / \delta_{obs_{sec_{1}}}$', fontsize=fsize)
plt.ylabel(r'$\delta_{obs}[ppm]$', fontsize=fsize)

"""
Now we plot the eta ratios
"""
nom_eta = (eta_t / eta_cob)
sec_eta = (eta_c / eta_cob_sec)
ext_eta = (eta_ext / eta_cob_ext)
ext_eta_correct = (eta_ext_correct / eta_cob_ext)

plt.figure(14)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    plt.plot(Pi, np.median(nom_eta[m]), 'ko', label='Nominal Mask')
    plt.plot(Pi, np.median(sec_eta[m]), 'ro', label='Secondary Mask')
    plt.plot(Pi, np.median(ext_eta[m]), 'bo', label='Extended Mask')
    plt.plot(Pi, np.median(ext_eta[m]), 'go', label='Extended Mask corrected')

plt.legend()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\eta_{flux} / \eta_{cob}$', fontsize=fsize)

plt.figure(15)
plt.plot(mag, eta_t/eta_cob, 'ko')
plt.plot(mag, eta_c/eta_cob_sec, 'ro')
plt.plot(mag, eta_ext/eta_cob_ext, 'bo')
plt.legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'],
       loc='best')
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\eta_{flux} / \eta_{cob}$', fontsize=fsize)


"""
Now we plot the delta_obs_ratio wit the median values
"""
plt.figure(16)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #plt.plot(Pi, np.median(delta_obs_c_2[m]), 'yo')
    #plot(Pi, np.median(sec_eta[m]), 'ro')
    #plot(Pi, np.median(ext_eta[m]), 'bo')
    #legend(['Nominal Mask', 'Secondary Mask', 'Extended Mask'], loc='best')
    #semilogy()
    #legend()

plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\delta_{obs_{sec}} / \delta_{obs_{sec_{1 pixel ring}}}$', fontsize=fsize)


"""
Now we plot the two expressions for the error of the COB (the ones )
"""
plt.figure(17)
plt.plot(mag, sigma_cob, 'ko')
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma_{\Delta C} [pixel]$', fontsize=fsize)
plt.title(r'COB error expression that does not dependend on $\delta_{back}$')

plt.figure(18)
#plt.plot(mag, sigma_cob_wrong, 'bo')
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma_{\Delta C} [pixel] $', fontsize=fsize)
plt.title(r'COB error expression that does depend on $\delta_{back}$')

"""
Trying to reproduce Réza's histogram
"""
plt.figure(19)
mask_eta = (eta_cob_ext > 3) & (eta_ext > 7.1)
r = (eta_ext / eta_cob_ext)[mask_eta]
plt.hist(r, bins=50, range=[0, 10])



""" 
Comparing double-aperture photometry with (nominal) COB shift
"""

plt.figure(20)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    plt.plot(Pi, ext, 'b+')
    #plot(Pi, ext_2, 'c+')
    #plot(Pi, ext_3, 'm+')
    #plt.plot(Pi, eff_cob, 'k^')
    #plt.plot(Pi, eff_cob_sec, 'r^')
    #plot(Pi, eff_cob_sec_2, 'g^')
    #plot(Pi, eff_cob*100, 'k')
    #plt.plot(Pi, eff_cob_ext, 'b^')
    #plot(Pi, eff_cob_ext_2, 'c^')
    #plot(Pi, eff_cob_ext_3, 'm^')
    plt.legend(['Sec. Mask Flux', 'Ext. Mask Flux'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison for the two types of DAP', fontsize=fsize)



"""
Now we plot the comparison between the flux and COB shift methods as a function of the target P magnitude
"""
plt.figure(21)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_corr = (edr_flux_correct[m].sum() / fp[m].sum()) * 100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    #plt.plot(Pi, ext, 'b+')
    plt.plot(Pi, ext_corr, 'b+')
    #plot(Pi, ext_2, 'c+')
    #plot(Pi, ext_3, 'm+')
    plt.plot(Pi, eff_cob, 'k^')
    plt.plot(Pi, eff_cob_sec, 'r^')
    #plot(Pi, eff_cob_sec_2, 'g^')
    #plot(Pi, eff_cob*100, 'k')
    plt.plot(Pi, eff_cob_ext, 'b^')
    #plot(Pi, eff_cob_ext_2, 'c^')
    #plot(Pi, eff_cob_ext_3, 'm^')
    plt.legend(['Sec. Mask Flux', 'Ext. Mask Flux', 'Nom. Mask COB shift',
            'Sec. Mask COB shift', 'Ext. Mask COB shift'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison of Every Method for Every Mask', fontsize=fsize)



"""
Now we plot the comparison between extended mask and the correct version of it as a function of the target P magnitude
"""
plt.figure(22)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    s = (eta_bt>flux_trsh)[m,:].sum()
    eff_cob = (ndr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_sec = (sdr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (edr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (sdr_flux[m].sum() / fp[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (edr_flux[m].sum() / fp[m].sum()) * 100
    ext_corr = (edr_flux_correct[m].sum() / fp[m].sum()) * 100
    ext_overall =((eta_ext_bt > flux_trsh) & (eta_bt>flux_trsh))[m, :].sum() / s*100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    #plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    plt.plot(Pi, ext, 'b+')
    plt.plot(Pi, ext_corr, 'g+')
    #plot(Pi, ext_2, 'c+')
    #plot(Pi, ext_3, 'm+')
    #plt.plot(Pi, eff_cob, 'k^')
    #plt.plot(Pi, eff_cob_sec, 'r^')
    #plot(Pi, eff_cob_sec_2, 'g^')
    #plot(Pi, eff_cob*100, 'k')
    #plt.plot(Pi, eff_cob_ext, 'b^')
    #plot(Pi, eff_cob_ext_2, 'c^')
    #plot(Pi, eff_cob_ext_3, 'm^')
    plt.legend(['Ext. Mask highest SPR', 'Ext. Mask all contaminants'], loc='best')

plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
plt.title('Efficieny Comparison for the Extended mask', fontsize=fsize)

plt.show()