import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from imagette import ran_unique_int

#dataDIR = '/home/fgutierrez/biruni3/Sep17_real_MC_T1413/test_results/'
dataDIR = '/home/fercho/double-aperture-photometry/test_results/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_000/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_001/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_002/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_003/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_004/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_005/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_006/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_007/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_008/'
#dataDIR = '/home/fercho/double-aperture-photometry/test_results/output_009/'
cataDIR = '/home/fercho/double-aperture-photometry/catalogues_stars/'
# Parameters for the plots
Pmin = 10
Pmax = 13
binsize = 0.5
nP = int((Pmax - Pmin) / binsize + 1)
fsize = 14
flux_trsh = 7.1
cob_trsh = 3
ntr = 3  # number of observed transit events
#td = 4  # transit event duration in hours
n_tar = 7000
# We load the npy files with all the metrics of the nominal and secondary and extended masks
#data_mag = np.load('SFP_DR3_20220831.npy')

data = np.load(dataDIR + 'targets_P5.npy')
#0: ID_t
#1: P_t
#2: n_c
#3: magnitude of the most sinificant contaminant
#4: distance from the target to the most significant contaminant
#5: w_t_key
#6: w_t_size 
#7: NSR1h
#8: n_bad
#9: SPR_crit
#10: sprk[ind_sprk]
#11: SPR_tot
#12: eta_t
#13: delta_obs_t
#14: abs_cob
#15: eta_cob
#16: sigma_1_24
#17-26: SPRK_10first
#27-36: eta_10first
#37: eta_true_positive_corrected_24_cameras
#38: eta_true_positive_corrected_6_cameras
#39: eta_true_positive_corrected_24_cameras_super-earth
#40: eta_true_positive_corrected_6_cameras_super-earth
#41: eta_cob_6_cameras
#42: abs_cob_6_cameras
#43: sigma_1_6_cameras
#44: SPR_crit_6_cameras
#45: eta_t_6_cameras
#46-55: eta_cob_10first
#56-65: sigma_cob_10first
#66-75: abs_cob_shift_10first
#76-85: eta_cob_10first_6_cameras
#86-95: sigma_cob_10first_6_cameras
#96-105: abs_cob_shift_10first_6_cameras
 
data_sec = np.load(dataDIR + 'targets_P5_secondary.npy')
#0: ID_t
#1: P_t
#2: secondary_mask_key
#3: secondary_mask_size
#4: nsr_1h_24_cameras_secondary_mask
#5: spr_tot_secondary_mask
#6: eta_c
#7: delta_obs_secondary_mask
#8: abs_cob_c
#9: eta_cob_c
#10: sigma_1_24_c
#11: eta_c_6_cameras
#12: eta_cob_c_6_cameras
#13: abs_cob_c_6_cameras
#14: sigma_1_6_cameras_c
#15: SPR_tot_sec_6_cameras
 
data_ext = np.load(dataDIR + 'targets_P5_extended.npy')
#0: ID_t
#1: P_t
#2: extended mask_key
#3: extended_mask_size
#4: NSR_ext_1h_24_cameras
#5: sprk_ext[ind_sprk]
#6: SPR_crit_ext
#7: eta_ext
#8: delta_obs_ext
#9: abs_cob_ext
#10: eta_cob_ext
#11: sigma_1_24_ext
#12: n_bad_ext
#13: SPR_tot_ext
#14-23: SPRK_ext_10first
#24-33: eta_ext_10first
#34-43: eta_ext_10first_6_cameras
#44: NSR_ext_1h_6_cameras
#45-54: eta_cob_ext_10first
#55-64: sigma_cob_ext_10first
#65-74: abs_cob_shift_ext_10first
#75-84: eta_cob_ext_10first_6_cameras
#85-94: sigma_cob_ext_10first_6_cameras
#95-104: abs_cob_shift_ext_10first_6_cameras

data_bray = np.load(dataDIR + 'targets_P5_bray.npy')
#0: ID_t
#1: P_t
#2: n_c
#3: NSR_bray_1h
#4: n_bad_bray
#5: SPR_crit_bray
#6: sprk_bray[ind_sprk]
#7: SPR_tot_bray

mag_value, star_count = np.loadtxt(dataDIR + 'star_count.txt', unpack=True, usecols=[0, 1])
#dback = 85000  # transit depth in ppm
#td = 4         # transit duration in hours
ntr = 3        # number of transits in one hour
n = data.shape[0]

dback_set = np.loadtxt(cataDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt')
dback_n = dback_set.shape[0]
seed = 123434434

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
eta_t_24_cameras = data[:, 37]
eta_t_6_cameras_earth_like = data[:, 38]
eta_t_24_cameras_superearth = data[:, 39]
eta_t_6_cameras_superearth = data[:, 40]
eta_t_6_cameras = data[:, 45]
delta_obs_t = data[:, 13]
eta_c = data_sec[:, 6]
eta_c_6_cameras = data_sec[:, 11]
delta_obs_c = data_sec[:, 7]
#delta_obs_c_6_cameras = data_sec[:, 13]
#eta_c_2 = data_sec[:, 15]
#delta_obs_c_2 = data_sec[:, 16]
eta_ext = data_ext[:, 7]
#eta_ext_correct = data_ext[:, 8]
delta_obs_ext = data_ext[:, 8]
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
eta_cob_6_cameras = data[:, 41]
delta_cob_6_cameras = data[:, 42]
sigma_cob_6_cameras = data[:, 43]
#eta_cob_wrong = data[:, 17]
#sigma_cob_wrong = data[:, 18]
delta_cob_sec = data_sec[:, 8]
eta_cob_sec = data_sec[:, 9]
sigma_cob_sec = data_sec[:, 10]
eta_cob_sec_6_cameras = data_sec[:, 12]
delta_cob_sec_6_cameras = data_sec[:, 13]
sigma_cob_sec_6_cameras = data_sec[:, 14]
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

# We get now the 10 first values for each param. of the nominal cob shift
eta_cob_10first = data[:, 46:56]
sigma_cob_10first = data[:, 56:66]
abs_cob_shift_10first = data[:, 66:76]
eta_cob_10first_6_cameras = data[:, 76:86]
sigma_cob_10first_6_cameras = data[:, 86:96]
abs_cob_shift_10first = data[:, 96:106]

# We get now the 10 first values for each param. of the extended cob shift
SPRK10_first = data[:,17:27]
SPRK10_first_ext = data_ext[:,14:24]
eta_cob_ext_10first_24_cameras = data_ext[:, 45:55]
sigma_cob_ext_10first_24_cameras = data_ext[:, 55:65]
delta_cob_ext_10first_24_cameras = data_ext[:, 65:75]
eta_cob_ext_10first_6_cameras = data_ext[:, 75:85]
sigma_cob_ext_10first_6_cameras = data_ext[:, 85:95]
delta_cob_ext_10first_6_cameras = data_ext[:, 95:105]

nbad_sp = np.zeros(n) # small planet 2<R_Ee -> 4*84ppm
eta_ext_bt = np.zeros((n, 10))
eta_ext_bt_6_cameras = np.zeros((n,10))
eta_bt = np.zeros((n, 10))
eta_bt_6_cameras = np.zeros((n, 10))
delta_obs = np.zeros((n,10))
delta_obs_ext = np.zeros((n,10))
delta_obs_ext_6_cameras = np.zeros((n, 10))

for i in range(n):
    j = ran_unique_int(10,interval=[0,dback_n-1]) # random sort of a BT (background transit)
    dback = dback_set[j,0] # transit depth
    td = dback_set[j,1] # transit duration
    #dback = np.ones(10)*85000
    #td = np.ones(10)*4.
    eta_bt[i, :] = (SPRK10_first[i, :]/data[i, 9])*flux_trsh*(dback/85000)*np.sqrt(td/4)
    eta_bt_6_cameras[i, :] = (SPRK10_first[i, :]/data[i, 44])*flux_trsh*(dback/85000)*np.sqrt(td/4)
    eta_ext_bt[i, :] = dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 4] * (1 - data_ext[i, 13]))
    eta_ext_bt_6_cameras[i, :] = dback*data_ext[i, 14:24]*np.sqrt(td*ntr)/(data_ext[i, 44] * (1 - data_ext[i, 13]))
    delta_obs[i,:] = dback*SPRK10_first[i,:] # observed transit depth
    #delta_int = delta_obs[i,:]/(1. -data_nommask[i,9] ) # inferred intrinsic transit depth
    delta_obs_ext[i,:] = dback*data_ext[i,14:24] # observed transit depth
     ##dback[:] = 85000
     ##td[:] = 4.
    delta_obs_ext_6_cameras[i, :] = dback*data_ext[i,14:24] # observed transit depth with 6 cameras
    delta_int = delta_obs_t[i]/ (1 - data[i, 11])
    nbad_sp[i] = np.sum( (eta_bt[i,:]>7.1) & (delta_int<4*84. ))

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

# First, the expressions for the efficiency false positives detections for the different masks
fp_24_cameras = (eta_t > flux_trsh)  # eta_t > 7.1: false positive
fp_6_cameras = (eta_t_6_cameras > flux_trsh)
secondary_mask_conditions_24_cameras = (eta_c > flux_trsh) & (delta_obs_c > delta_obs_t) & fp_24_cameras  # secondary mask efficiency condition for 24 cameras
secondary_mask_conditions_6_cameras = (eta_c_6_cameras > flux_trsh) & (delta_obs_c > delta_obs_t) & fp_6_cameras # secondary mask efficiency condition for 6 cameras
#s_2dr_flux = fp & (eta_c_2 > flux_trsh) & (delta_obs_c_2 > delta_obs_t)  # secondary mask false positive detection rate
efd = fp_24_cameras & (eta_ext > flux_trsh) & (data_ext[:,8] > data[:,13])  # extended mask false positive detection rate
#efd_correct = fp & (eta_ext_correct > flux_trsh) & (delta_obs_ext > delta_obs_t)  # extended mask false positive detection rate
#e_2dr_flux = fp & (eta_ext_2 > flux_trsh) & (delta_obs_ext_2 > delta_obs_t)  # extended (2) mask false positive detection rate
#e_3dr_flux = fp & (eta_ext_3 > flux_trsh) & (delta_obs_ext_3 > delta_obs_t)  # extended (3) mask false positive detection rate
#efd_ext_overall = fp & (eta_ext_bt>flux_trsh)
# Second, the expressions for computing the efficiency for detecting false positives of the cob shift for different masks
cd = fp_24_cameras & (eta_cob > cob_trsh)  # nominal mask false positive detection rate via cob shift
cd_6_cameras = fp_6_cameras & (eta_cob_6_cameras > cob_trsh)
secondary_mask_conditions_cob_24_cameras = (eta_cob_sec > cob_trsh) & fp_24_cameras  # secondary mask efficiency condition for 24 cameras and cob shift
secondary_mask_conditions_cob_6_cameras = (eta_cob_sec_6_cameras > cob_trsh) & fp_6_cameras # secondary mask efficiency condition for 6 cameras and cob shift
#s_2dr_cob = fp & (eta_cob_sec_2 > cob_trsh)  # secondary mask false positive detection rate via cob shift
ecd = fp_24_cameras & (eta_cob_ext > cob_trsh)  # extended mask false positive detection rate via cob shift
#e_2dr_cob = fp & (eta_cob_ext_2 > cob_trsh)  # extended mask (2) false positive detection rate via cob shift
#e_3dr_cob = fp & (eta_cob_ext_3 > cob_trsh)  # extended mask (3) false positive detection rate via cob shift


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

plt.figure(7)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (cd[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_24_cameras[m].sum()) * 100
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
    #plt.plot(Pi, np.median(eta_t[m]), 'ko')
    #plt.plot(Pi, np.median(eta_c[m]), 'ro')
    plt.plot(Pi, np.median(eta_t_24_cameras),  'ko')
    plt.plot(Pi, np.median(eta_t_6_cameras_earth_like),  'ro')
    # plot(mag, eta_c_2, 'g+', label='Secondary Mask (1 pixel ring)')
    #plt.plot(Pi, np.median(eta_ext[m]), 'bo')
    #plt.plot(Pi, np.median(eta_ext_correct[m]), 'go')
    # plot(mag, eta_ext_2, 'c+', label='Extended Mask (2 pixels ring)')
    # plot(mag, eta_ext_3, 'm+', label='extended mask (3)')
    plt.legend(['24 cameras', '6 cameras'], loc='best')

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

def create_plot_11():
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
    eff_cob = (cd[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (efd[m].sum() / fp_24_cameras[m].sum()) * 100
    #ext_corr = (efd_correct[m].sum() / fp[m].sum()) * 100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    plt.plot(Pi, ext, 'b+')
    #plt.plot(Pi, ext_corr, 'g+')
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
#ext_eta_correct = (eta_ext_correct / eta_cob_ext)

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
Now we plot the two expressions for the error of the COB (the ones )
"""
plt.figure(17)
plt.plot(mag, sigma_cob, 'ko')
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$\sigma_{\Delta C} [pixel]$', fontsize=fsize)
plt.title(r'COB error expression that does not dependend on $\delta_{back}$')



""" 
Comparing double-aperture photometry with (nominal) COB shift
"""

plt.figure(20)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    eff_cob = (cd[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_24_cameras[m].sum()) * 100
    sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    ext = (efd[m].sum() / fp_24_cameras[m].sum()) * 100
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
    eff_cob = (cd[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_sec_2 = (s_2dr_cob[m].sum() / fp[m].sum()) * 100
    eff_cob_ext = (ecd[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_ext_2 = (e_2dr_cob[m].sum() / fp[m].sum()) * 100
    #eff_cob_ext_3 = (e_3dr_cob[m].sum() / fp[m].sum()) * 100
    sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    sec_6 = (secondary_mask_conditions_6_cameras[m].sum() / fp_6_cameras[m].sum()) * 100
    #sec_2 = (s_2dr_flux[m].sum() / fp[m].sum()) * 100
    ext = (efd[m].sum() / fp_24_cameras[m].sum()) * 100
    #ext_corr = (efd_correct[m].sum() / fp[m].sum()) * 100
    #ext_2 = (e_2dr_flux[m].sum() / fp[m].sum()) * 100
    #ext_3 = (e_3dr_flux[m].sum() / fp[m].sum()) * 100
    plt.plot(Pi, sec, 'r+')
    #plot(Pi, sec_2, 'g+')
    #plt.plot(Pi, ext, 'b+')
    #plt.plot(Pi, ext_corr, 'b+')
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
    fp_ext_overall_24_cameras = (eta_bt>flux_trsh)[m,:].sum()
    fp_ext_overall_6_cameras = (eta_bt_6_cameras>flux_trsh)[m,:].sum()
    eff_ext_overall = ((eta_ext_bt > flux_trsh) & (delta_obs_ext>delta_obs) & (eta_bt>flux_trsh))[m,:].sum() / fp_ext_overall_24_cameras*100.
    eff_ext_overall_6_cameras = ((eta_ext_bt_6_cameras > flux_trsh) & (delta_obs_ext_6_cameras>delta_obs) & (eta_bt_6_cameras>flux_trsh))[m,:].sum() / fp_ext_overall_6_cameras*100.
    eff_sec_6_cameras = (secondary_mask_conditions_6_cameras[m].sum() / fp_6_cameras[m].sum()) * 100
    eff_sec = (secondary_mask_conditions_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100

    error = np.sqrt(n_tar * eff_ext_overall * (100 - eff_ext_overall)) / n_tar
    error_6_cameras = np.sqrt(n_tar * eff_ext_overall_6_cameras * (100 - eff_ext_overall_6_cameras)) / n_tar
    error_sec_6_cameras = np.sqrt(n_tar * eff_sec_6_cameras * (100 - eff_sec_6_cameras)) / n_tar
    error_sec = np.sqrt(n_tar * eff_sec * (100 - eff_sec)) / n_tar
    
    
    # Plotting with linestyle='-'
    plt.errorbar(Pi, eff_sec, yerr=error_sec, fmt='o', color='purple', ecolor='purple', capsize=5, label='Sec. Mask (24 cameras)' if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_sec_6_cameras, yerr=error_sec_6_cameras, fmt='o',  color='green', ecolor='green', capsize=5, label='Sec. Mask (6 cameras)' if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall, yerr=error, fmt='s',  color='blue', ecolor='blue', capsize=5, label='Ext. Mask (24 cameras)'  if i == 0 else "", markersize=4)
    plt.errorbar(Pi, eff_ext_overall_6_cameras, yerr=error_6_cameras, fmt='s',  color='red', ecolor='red', capsize=5, label='Ext. Mask (6 cameras)'  if i == 0 else "", markersize=4)
    plt.fill_between([9, 11.7], [20, 20], [100, 100], color='aqua', alpha=0.1)
    plt.fill_between([11, 13.4], [20,20], [100, 100], color='plum', alpha=0.1)

    # Plot lines connecting points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_sec, eff_sec], color='purple', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_sec_6_cameras, eff_sec_6_cameras], color='green', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall, eff_ext_overall], color='blue', linestyle='-', markersize=0)
        plt.plot([prev_Pi, Pi], [prev_eff_ext_overall_6_cameras, eff_ext_overall_6_cameras], color='red', linestyle='-', markersize=0)
    
    prev_Pi, prev_eff_sec, prev_eff_sec_6_cameras, prev_eff_ext_overall, prev_eff_ext_overall_6_cameras = Pi, eff_sec, eff_sec_6_cameras, eff_ext_overall, eff_ext_overall_6_cameras

    #plt.legend(['Ext. Mask (10 contaminants)'], loc='best')
    plt.vlines(11.7, ymin=20, ymax = 100, linestyles='dashed', colors='green')
    plt.vlines(11, ymin=20, ymax=100, linestyles='dashdot', colors='red')
    plt.ylim(20, 100)
    plt.xlim(9.9, 13.4)
    plt.text(10, 70.1,'Earth-like planet detection\nregion (24 cameras)', color='green', weight='bold')
    plt.text(11, 54, 'On-board light curve processing region', color='red',  weight='bold')

# Display legend
plt.legend()
plt.xlabel('P Magnitude', fontsize=fsize)
plt.ylabel('Efficiency[%]', fontsize=fsize)
#plt.title('Double-Aperture Photometry Comparison', fontsize=fsize)

"""
Now we plot the statistical significance for every mask as a function of the target P magnitude
"""
plt.figure(23)

plt.plot(mag, eta_t_24_cameras, 'o', markersize=4,  label='Earth-like planets (24 cameras)', )
#plt.plot(mag, eta_t_24_cameras_superearth, 'o', markersize=4,  label='Super-Earths (24 cameras)')
plt.plot(mag, eta_t_6_cameras_earth_like, 'o', markersize=4, label='Jovian planets (6 cameras)')
plt.plot(mag, eta_t_6_cameras_superearth, 'o', markersize=4, label='Super-Earth (6 cameras)')
plt.hlines(7.1, xmin=10, xmax=13, linestyles='dashed', colors='red')
plt.legend()
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel(r'$ \eta $', fontsize=fsize)
plt.xlim(10, 13)


"""
Now we plot the efficiency of the COB shift (all contaminants)
    
"""
plt.figure(24)
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    s_24_cameras = (eta_bt>flux_trsh)[m,:].sum()
    s_6_cameras = (eta_bt_6_cameras>flux_trsh)[m,:].sum()
    #eff_ext = (efd[m].sum() / fp[m].sum()) * 100
    eff_ext_cob_overall = ((eta_cob_ext_10first_24_cameras > cob_trsh) & (eta_bt>flux_trsh))[m,:].sum() / s_24_cameras*100.
    eff_ext_cob_overall_6_cameras = ((eta_cob_ext_10first_6_cameras > cob_trsh) &  (eta_bt_6_cameras>flux_trsh))[m,:].sum() / s_6_cameras*100.
    eff_cob = ((eta_cob_10first > cob_trsh) & (eta_bt>flux_trsh))[m,:].sum() / s_24_cameras*100
    eff_cob_6_cameras = ((eta_cob_10first_6_cameras > cob_trsh) & (eta_bt>flux_trsh))[m,:].sum() / s_6_cameras*100
    #eff_cob = (cd[m].sum() / fp_24_cameras[m].sum()) * 100
    #eff_cob_6_cameras = (cd_6_cameras[m].sum() / fp_6_cameras[m].sum()) * 100
    eff_cob_sec = (secondary_mask_conditions_cob_24_cameras[m].sum() / fp_24_cameras[m].sum()) * 100
    eff_cob_sec_6_cameras = (secondary_mask_conditions_cob_6_cameras[m].sum() / fp_6_cameras[m].sum()) * 100
    
    #Computing the errors
    error_ext_cob = np.sqrt(n_tar * eff_ext_cob_overall * (100 - eff_ext_cob_overall)) / n_tar
    error_ext_cob_6_cameras = np.sqrt(n_tar * eff_ext_cob_overall_6_cameras * (100 - eff_ext_cob_overall_6_cameras)) / n_tar
    error_cob = np.sqrt(n_tar * eff_cob * (100 - eff_cob)) / n_tar
    error_cob_6_cameras = np.sqrt(n_tar * eff_cob_6_cameras * (100 - eff_cob_6_cameras)) / n_tar
    error_cob_sec = np.sqrt(n_tar * eff_cob_sec * (100 - eff_cob_sec)) / n_tar
    error_cob_sec_6_cameras = np.sqrt(n_tar * eff_cob_sec_6_cameras * (100 - eff_cob_sec_6_cameras)) / n_tar
    
    plt.errorbar(Pi, eff_ext_cob_overall, fmt='s', yerr=error_ext_cob, label='Ext. Mask (24 cameras)' if i == 0 else "", color='blue', ecolor='blue', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_ext_cob_overall_6_cameras, fmt='s', yerr=error_ext_cob_6_cameras, label='Ext. Mask (6 cameras)' if i == 0 else "", color='red', ecolor='red', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob, fmt='*', yerr=error_cob, label='Nom. Mask (24 cameras)' if i == 0 else "", color='orange', ecolor='orange', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_6_cameras, fmt='*', yerr=error_cob_6_cameras, label='Nom. Mask (6 cameras)' if i == 0 else "", color='olive', ecolor='olive', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec, fmt='o', yerr=error_cob_sec, label='Sec. Mask (24 cameras)' if i == 0 else "", color='purple', ecolor='purple', capsize=5, markersize=4)
    plt.errorbar(Pi, eff_cob_sec_6_cameras, fmt='o', yerr=error_cob_sec_6_cameras, label='Sec. Mask (6 cameras)' if i == 0 else "", color='green', ecolor='green', capsize=5, markersize=4)
    plt.fill_between([9, 11.7], [20, 20], [100, 100], color='aqua', alpha=0.1)
    plt.fill_between([11, 13.4], [20,20], [100, 100], color='plum', alpha=0.1)
    
    # Plot lines connecting the points
    if i > 0:
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall, eff_ext_cob_overall], color='blue', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob, eff_cob], color='orange', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec, eff_cob_sec], color='purple', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_ext_cob_overall_6_cameras, eff_ext_cob_overall_6_cameras], color='red', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_6_cameras, eff_cob_6_cameras], color='olive', linestyle='-')
        plt.plot([prev_Pi, Pi], [prev_eff_cob_sec_6_cameras, eff_cob_sec_6_cameras], color='green', linestyle='-')   
    # Update previous values
    prev_Pi, prev_eff_ext_cob_overall, prev_eff_cob, prev_eff_cob_sec, prev_eff_ext_cob_overall_6_cameras, prev_eff_cob_6_cameras, prev_eff_cob_sec_6_cameras = Pi, eff_ext_cob_overall, eff_cob, eff_cob_sec, eff_ext_cob_overall_6_cameras, eff_cob_6_cameras, eff_cob_sec_6_cameras


    #plt.legend(['Ext. Mask (10 contaminants)'], loc='best')
    plt.vlines(11.7, ymin=20, ymax = 100, linestyles='dashed', colors='green')
    plt.vlines(11, ymin=20, ymax=100, linestyles='dashdot', colors='red')
    plt.ylim(20, 100)
    plt.xlim(9.9, 13.4)
    plt.text(10, 91.1,'Earth-like planet detection \nregion (24 cameras)', color='green', weight='bold')
    plt.text(11, 65, 'On-board light curve processing region', color='red', weight='bold')
    
    
    
plt.legend()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel('Efficiency [%]', fontsize=fsize)


"""" 
Now we check the n_bad variable with the new nbad_sp parameter
"""
plt.figure(25)

plt.plot(mag, nbad_sp, 'o', label='Nbad_sp')
#plt.plot(mag, n_bad, 'o', label='Nbad')

plt.legend()
plt.semilogy()
plt.xlabel('P magnitude', fontsize=fsize)
plt.ylabel('counts')


""" 
Now we make the histograms for the conditions
"""
plt.figure(26)
mag_value = 12.5

# Lists to store data for plotting histograms
count_list = []
count_list_6_cameras = []

for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #if Pi == mag_value:
    condition_eta_ext = (eta_ext_bt > flux_trsh) & (delta_obs_ext>delta_obs)
    filtered_eta_ext = condition_eta_ext[m, :]
    count = filtered_eta_ext.sum()
    count_list.append(count) 
    condition_eta_ext_6_cameras = (eta_ext_bt_6_cameras > flux_trsh) & (delta_obs_ext>delta_obs)
    filtered_eta_ext_6_cameras = condition_eta_ext_6_cameras[m, :]
    count_6_cameras = filtered_eta_ext_6_cameras.sum()
    count_list_6_cameras.append(count_6_cameras)
        
  

plt.hist(count_list, bins='auto', alpha=0.5, label='eta_ext')
plt.hist(count_list_6_cameras, bins='auto', alpha=0.5, label='eta_ext_6_cameras')
#plt.hist(condition_eta_ext_list, label='eta_ext', alpha=0.5)
#plt.hist(condition_eta_ext_6_cameras_list, label='eta_ext_6_cameras', alpha=0.5)  
plt.xlabel('Number of Targets Fulfilling Condition')
plt.ylabel('Frequency')
plt.legend()    
         

"""
The condition of eta_ext
"""
plt.figure(28)
mag_value = 12.5

# Lists to store data for plotting histograms
count_eta_ext = []
count_eta_ext_6_cameras = []
magnitude_ranges = []
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #if Pi == mag_value:
    condition_eta_ext = (delta_obs_ext>delta_obs) & (eta_ext_bt>flux_trsh)
    filtered_eta_ext = condition_eta_ext[m, :]
    count_ext = filtered_eta_ext.sum()
    count_eta_ext.append(count_ext) 
    condition_eta_ext_6_cameras = (delta_obs_ext>delta_obs) & (eta_ext_bt_6_cameras>flux_trsh)
    filtered_eta_ext_6_cameras = condition_eta_ext_6_cameras[m, :]
    count_ext_6_cameras = filtered_eta_ext_6_cameras.sum()
    count_eta_ext_6_cameras.append(count_ext_6_cameras)
    magnitude_ranges.append(str(Pi))

plt.bar(magnitude_ranges, count_eta_ext, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{ext} \geq \eta_{min})$ (24 cameras)')
plt.bar(magnitude_ranges, count_eta_ext_6_cameras, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{ext} \geq \eta_{min})$ (6 cameras)')
#plt.hist(condition_eta_ext_list, label='eta_ext', alpha=0.5)
#plt.hist(condition_eta_ext_6_cameras_list, label='eta_ext_6_cameras', alpha=0.5)
plt.legend()  
plt.xlabel('Magnitude Bin')
plt.ylabel('Counts')
plt.legend()



plt.figure(29)
mag_value = 12.5

# Lists to store data for plotting histograms
count_eta_nom = []
count_eta_nom_6_cameras = []
magnitude_ranges = []
for i in range(nP):
    Pi = Pmin + i * binsize
    m = (mag >= Pi - binsize/2.) & (mag <= Pi + binsize/2.)
    #if Pi == mag_value:
    condition_eta_nom = (delta_obs_ext>delta_obs) & (eta_bt>flux_trsh)
    filtered_eta_nom = condition_eta_nom[m, :]
    count_nom = filtered_eta_nom.sum()
    count_eta_nom.append(count_nom) 
    condition_eta_nom_6_cameras = (delta_obs_ext>delta_obs) & (eta_bt_6_cameras>flux_trsh)
    filtered_eta_nom_6_cameras = condition_eta_nom_6_cameras[m, :]
    count_nom_6_cameras = filtered_eta_nom_6_cameras.sum()
    count_eta_nom_6_cameras.append(count_nom_6_cameras)
    magnitude_ranges.append(str(Pi))
        
  

plt.bar(magnitude_ranges, count_eta_nom, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$ (24 cameras)')
plt.bar(magnitude_ranges, count_eta_nom_6_cameras, alpha=0.5, label=r'$(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$ (6 cameras)')
#plt.hist(condition_eta_ext_list, label='eta_ext', alpha=0.5)
#plt.hist(condition_eta_ext_6_cameras_list, label='eta_ext_6_cameras', alpha=0.5)
plt.legend()  
plt.xlabel('Magnitude Bin')
plt.ylabel('Counts')
plt.legend()



"""
The histogram of eta_ext
"""
plt.figure(30)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

eta_and_delta_mask = (eta_bt[mag_mask] >= flux_trsh) & (delta_obs_ext[mag_mask] >= delta_obs[mag_mask])
eta_and_delta_mask_6_cameras = (eta_bt_6_cameras[mag_mask] >= flux_trsh) & (delta_obs_ext[mag_mask] >= delta_obs[mag_mask])

final_eta_ext = eta_ext[eta_and_delta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[eta_and_delta_mask_6_cameras]

number_below_threshold = len(np.where(final_eta_ext < flux_trsh)[0])
number_above_threshold = len(np.where(final_eta_ext > flux_trsh)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=150, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 100, number_below_threshold, color='blue', weight='bold')
plt.text(50, 100, number_above_threshold, color='blue', weight='bold')
plt.text(0, 70, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 70, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given $(\delta_{back}^{ext} \geq \delta_{back}^{nom}) & (\eta_{nom} \geq \eta_{min})$')
plt.ylabel('Counts')
plt.ylim(0, 150)
plt.legend()



plt.figure(31)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_nom = eta_bt[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]
eta_nom_6_cameras = eta_bt_6_cameras[mag_mask]
eta_ext_flat = eta_ext.flatten()
eta_nom_flat = eta_nom.flatten()
eta_ext_flat_6_cameras = eta_ext_6_cameras.flatten()
eta_nom_flat_6_cameras = eta_nom_6_cameras.flatten()
plt.plot(eta_ext_flat/eta_nom_flat, 'o', alpha=0.5)
plt.plot(eta_ext_flat_6_cameras/eta_nom_flat_6_cameras, 'o', alpha=0.5)
#plt.ylim(0,20)
plt.grid(True)
#plt.plot(eta_nom, eta_ext)


plt.figure(32)

mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_nom_24_cameras = eta_bt[mag_mask]
eta_nom_6_cameras = eta_bt_6_cameras[mag_mask]

eta_nom_24_cameras_flat = eta_nom_24_cameras.flatten()

eta_nom_6_cameras_flat = eta_nom_6_cameras.flatten()


number_below_threshold = len(np.where(eta_nom_24_cameras_flat < 7.1)[0])
number_above_threshold = len(np.where(eta_nom_24_cameras_flat > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(eta_nom_6_cameras_flat < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(eta_nom_6_cameras_flat > 7.1)[0])

plt.hist(eta_nom_24_cameras_flat, bins=15, range=(0, 50), alpha=0.5, label='24 cameras')
plt.hist(eta_nom_6_cameras_flat, bins=15, range=(0,50), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=19800, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
plt.text(0, 17500, number_below_threshold, color='blue', weight='bold')
plt.text(50, 17500, number_above_threshold, color='blue', weight='bold')
plt.text(0, 12500, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 12500, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{k}^{nom}$')
plt.legend()


plt.figure(33)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_c = eta_c[mag_mask]
plt.scatter(eta_c, eta_ext[:, 0], alpha=0.5)
plt.ylim(0, 20)
plt.xlim(0, 20)
plt.xlabel(r'$\eta_{kmax}^{sec}$')
plt.ylabel(r'$\eta_{ext}$')


plt.figure(34)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
delta_obs_ext = delta_obs_ext[mag_mask]
delta_obs_c = delta_obs_c[mag_mask]
plt.scatter(delta_obs_c, delta_obs_ext[:, 0], alpha=0.5)
plt.xlabel(r'$\delta_{kmax}^{sec}$')
plt.ylabel(r'$\delta_{ext}$')


"""
The histogram of eta_ext
"""
plt.figure(35)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

eta_mask = (eta_bt[mag_mask] >= flux_trsh)
eta_mask_6_cameras = (eta_bt_6_cameras[mag_mask] >= flux_trsh)

final_eta_ext = eta_ext[eta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[eta_mask_6_cameras]

number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=150, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 100, number_below_threshold, color='blue', weight='bold')
plt.text(50, 100, number_above_threshold, color='blue', weight='bold')
plt.text(0, 70, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 70, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given $(\eta_{nom} \geq \eta_{min})$')
plt.ylabel('Counts')
plt.ylim(0, 150)
plt.legend()


"""
The histogram of eta_ext
"""
plt.figure(36)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

delta_mask = (delta_obs_ext >= delta_obs[mag_mask])

final_eta_ext = eta_ext[delta_mask]
final_eta_ext_6_cameras = eta_ext_6_cameras[delta_mask]

number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=22000, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
#plt.plot(eta_nom, eta_ext)
plt.text(0, 1000, number_below_threshold, color='blue', weight='bold')
plt.text(50, 1000, number_above_threshold, color='blue', weight='bold')
plt.text(0, 7000, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 7000, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$ given  $(\delta_{back}^{ext} \geq \delta_{back}^{nom})$')
plt.ylabel('Counts')
#plt.ylim(0, 200)
plt.legend()



plt.figure(37)
mag_mask = (mag >= 12.25) & (mag <= 12.75)
eta_ext = eta_ext_bt[mag_mask]
eta_ext_6_cameras = eta_ext_bt_6_cameras[mag_mask]

final_eta_ext = eta_ext.flatten()
final_eta_ext_6_cameras = eta_ext_6_cameras.flatten()


number_below_threshold = len(np.where(final_eta_ext < 7.1)[0])
number_above_threshold = len(np.where(final_eta_ext > 7.1)[0])

number_below_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras < 7.1)[0])
number_above_threshold_6_cameras = len(np.where(final_eta_ext_6_cameras > 7.1)[0])

plt.hist(final_eta_ext, bins=51, range=(0, 100), alpha=0.5, label='24 cameras')
plt.hist(final_eta_ext_6_cameras, bins=51, range=(0,100), alpha=0.5, label='6 cameras')
plt.vlines(7.1, ymin=0, ymax=20210, linestyles='dashdot', colors='red', label=r'$\eta_{min} = 7.1$')
plt.text(0, 17000, number_below_threshold, color='blue', weight='bold')
plt.text(50, 17000, number_above_threshold, color='blue', weight='bold')
plt.text(0, 16000, number_below_threshold_6_cameras, color='orange', weight='bold')
plt.text(50, 16000, number_above_threshold_6_cameras, color='orange', weight='bold')
plt.xlabel(r'$\eta_{ext}$')
plt.ylabel('Counts')
#plt.ylim(0, 200)
plt.legend()


l1 = (delta_obs_ext >= delta_obs[mag_mask]) & np.logical_and(flux_trsh <= eta_ext, eta_ext <= 2*flux_trsh) & (eta_bt[mag_mask] >= 2 *flux_trsh)
l2 = (delta_obs_ext >= delta_obs[mag_mask]) & np.logical_and(flux_trsh <= eta_bt[mag_mask], eta_bt[mag_mask] <= 2*flux_trsh) & (eta_ext >= 2*flux_trsh)
l3 = (delta_obs_ext >= delta_obs[mag_mask]) & np.logical_and(flux_trsh <= eta_bt[mag_mask], eta_bt[mag_mask] <= 2*flux_trsh) & np.logical_and(flux_trsh <= eta_ext, eta_ext <= (2*flux_trsh))
l4 = np.logical_and(flux_trsh < eta_bt[mag_mask], eta_bt[mag_mask] < 2*flux_trsh)


print('L1 is:', np.sum(l1))
print('L2 is:', np.sum(l2))
print('L3 is:',np.sum(l3))
print('L4 is:',np.sum(l4))
print('L1 + L2 + L3 is:', np.sum(l1 + l2 + l3))





plt.show()
 