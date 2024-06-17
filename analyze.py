from pylab import *
from pylab import *
from imagette import ran_unique_int

# Dir = '10000/'
Dir = 'PSF_Focus_0mu_0.2pxdif/'
Pmin = 8. # 10.5
Pmax = 13
binsize = 0.5
cob_thr= 3.
ntr = 3  # number of transits
flx_trh = 7.1

data_nommask = np.load(Dir+'targets_P5.npy')
# 0: ID_t
# 1: P_t
# 2: PSF index
# 3: n_c
# 4: w_t_key
# 5: w_t_size
# 6: NSR1h
# 7: n_bad
# 8: SPR_crit
# 9: SPR_tot
# 10: ID_c
# 11: P_c
# 12: spr max
# 13: eta_t
# 14: delta_obs_t
# 15: delta_COB
# 16: delta_COB_sig_1h_24c
# 17: eta_cob
# 18: d_c
# 19: Gammax
# 20: Gammay
# 21: n_bad_wrong (Victor's formula)
# 22-31: 10 first SPRk values
# 32-41: 10 first Gamma values
# 42-51: 10 first delta_COB_sig_1h_24c

data_2ndmask = np.load(Dir+'targets_P5_2ndmask.npy')
# 0: ID_t
# 1: P_t
# 2: ID_c
# 3: P_c
# 4: NSR1h_c
# 5: w_c_key
# 6: w_c_size
# 7: eta_c
# 8: delta_obs_c
# 9: delta_COB_c
# 10: delta_COB_sig_1h_24c_c
# 11: eta_cob_c
# 12: spr_tot_c
# 13: Gammax_c
# 14: Gammay_c



data_extmask = np.load(Dir+'targets_P5_extended.npy')
# 0: ID_t
# 1: P_t
# 2: w_ext_key
# 3: w_ext_size
# 4: NSR_ext_1h
# 5: n_bad_ext
# 6: SPR_crit_ext
# 7: SPR_tot_ext
# 8: ID_c
# 9: P_c
# 10: eta_ext
# 11: delta_obs_ext
# 12: delta_COB_ext
# 13: delta_COB_sig_1h_24c_ext
# 14: eta_cob_ext
# 15: delta_obs_ext
# 16: eta_ext_max
# 17: delta_obs_ext_max
# 18: n_det_ext
# 19: Gammax_ext
# 20: Gammay_ext
# 21-30: 10 first SPRk values
# 31-40: 10 first Gamma values
# 41-50: 10 first delta_COB_sig_1h_24c

nP = int(round ( (Pmax-Pmin)/binsize + 1 ))

j = np.argsort(data_nommask[:,1])
data_nommask = data_nommask[j]
data_2ndmask = data_2ndmask[j]
data_extmask = data_extmask[j]

P = data_nommask[:,1]
m = (P>=Pmin-binsize/2.) & (P<Pmax+binsize/2.)

data_nommask = data_nommask[m]
data_2ndmask = data_2ndmask[m]
data_extmask = data_extmask[m]
P = P[m]

fp = (data_nommask[:,13] > flx_trh) # eta_t>7: false positive (not yet identified as such)
sfd = fp & (data_2ndmask[:,7] > flx_trh)  &  (data_2ndmask[:,8] > data_nommask[:,14])  # eta_c>7.1 and delta_obt_c > delta_obs_t : secondary mask detection
efd = fp & (data_extmask[:,10] > flx_trh)   & (data_extmask[:,15] > data_nommask[:,14])  # eta_ext>7.1 and delta_obt_ext > delta_obs_t : extended flux detection

cd = fp & (data_nommask[:,17] > cob_thr) # eta_cob>3: COB detection
scd = fp & (data_2ndmask[:,11] > cob_thr)  # eta_cob_c>3: s-COB detection
ecd = fp & (data_extmask[:,14] > cob_thr)  # eta_cob_ext>3: e-COB detection



n = data_nommask.shape[0]

SPRk = data_nommask[:,22:32]
# SPRk = data_extmask[:,22:32]
dback_set = np.loadtxt('KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt')
dback_n = dback_set.shape[0]
seed = 123434434
np.random.seed(seed)
nbad = np.zeros(n)
nbad_cob = np.zeros(n)

nbad_sp = np.zeros(n) # small planet 2<R_Ee -> 4*84ppm

nbad6c = np.zeros(n)

delta_obs = np.zeros((n,10))
delta_obs_ext = np.zeros((n,10))
eta_bt = np.zeros((n,10))
eta_cob_bt = np.zeros((n,10))
eta_cob_ext_bt = np.zeros((n,10))
eta_ext_bt = np.zeros((n,10))
ntr = 3  # number of transits
nbad_ext = np.zeros(n)
nbad_ext_cob = np.zeros(n)

# compute n_bad  for a real distribution in delta_back: nb of contaminant stars eta>eta_min=7.1
for i in range(n):
    j = ran_unique_int(10,interval=[0,dback_n-1]) # random sort of a BT (background transit)
    dback = dback_set[j,0] # transit depth
    td = dback_set[j,1] # transit duration
    # dback = np.ones(10)*85000
    # td = np.ones(10)*4.
    eta_bt[i,:] = (SPRk[i,:]/data_nommask[i,8])*flx_trh *(dback/85000)*sqrt(td/4.) # significance of the BT in the nominal flux
    eta_ext_bt[i,:] = dback*data_extmask[i,21:31]*np.sqrt(td*ntr)/(1-data_extmask[i,7])/data_extmask[i,4] # significance in the extended mask
    nbad[i] = np.sum(eta_bt[i,:]>flx_trh) # number of false detection in the nominal mask
    nbad_ext[i] = np.sum(eta_ext_bt[i,:]>flx_trh)  # number of false detection in the extended mask
    delta_obs[i,:] = dback*SPRk[i,:] # observed transit depth
    delta_int = delta_obs[i,:]/(1. -data_nommask[i,9] ) # inferred intrinsic transit depth
    delta_obs_ext[i,:] = dback*data_extmask[i,21:31] # observed transit depth
    nbad_sp[i] =     np.sum( (eta_bt[i,:]>flx_trh) & (delta_int<4*84. ))
    eta6c = eta_bt[i,:]*sqrt(6./24) #  significance for 6 cameras
    nbad6c[i] = np.sum(eta6c>flx_trh) # number of BT detected with 6 cameras
##    dback[:] = 85000
    Lambda =  dback*1e-6/(1.-dback*1e-6*data_nommask[i,22:32])
##    td[:] = 4.
    eta_cob_bt[i,:] = Lambda*data_nommask[i,32:42]*sqrt(td*ntr)/data_nommask[i,42:52] # significance of centroid shift in the nominal mask
    nbad_cob[i] = np.sum(eta_cob_bt[i,:]>cob_thr)
    Lambda =  dback*1e-6/(1.-dback*1e-6*data_extmask[i,21:31])
    eta_cob_ext_bt[i,:] = Lambda*data_extmask[i,31:41]*sqrt(td*ntr)/data_extmask[i,41:51] # significance of centroid shift in the extended mask
    nbad_ext_cob[i] = np.sum(eta_cob_ext_bt[i,:]>cob_thr)
#STOP
# shapes number
figure(0)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m =  (P < Pi + binsize/2.) & fp
    scatter([Pi],[len(np.unique(np.array(data_nommask[m,4],dtype=np.int64)))],color='k')
    scatter([Pi],[len(np.unique(np.array(data_2ndmask[m,5],dtype=np.int64)))],color='r')
    scatter([Pi],[len(np.unique(np.array(data_extmask[m,2],dtype=np.int64)))],color='b')
    for k in range(1,4):
        p =  m.sum()//k
        j = np.array(np.arange(p)*m.sum()/p,dtype=np.int32)
        t =     len(np.unique(np.array(data_2ndmask[m,5][j],dtype=np.int64)))
        print(Pi,k,p,m.sum(),t)

    #      scatter([Pi],[len(np.unique(np.array(data_2ndmask[m,5][j],dtype=np.int64)))],color='r')




# ,label='fraction: %f' % (1./k)

# plot(P,nunique_nomask,'k+')
# plot(P,nunique_2ndmask,'r+')
# plot(P,nunique_extmask,'b+')
xlabel('P')
ylabel(r'Cumul. number of Mask shapes')
title(Dir)


# Mask size
figure(1)
clf()
plot(P,data_nommask[:,5],'k+')
plot(P,data_2ndmask[:,6],'r+')
plot(P,data_extmask[:,3],'b+')
xlabel('P')
ylabel(r'Mask size')
title(Dir)

# efficiency: flux
figure(2)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s = fp[m].sum()
    eff_s = sfd[m].sum()/s * 100.
    eff_ext = efd[m].sum()/s * 100.
    scatter([Pi],[eff_ext],color='b')
    scatter([Pi],[eff_s],color='r')
xlabel('P')
ylabel(r'Efficieny [%%]')
title(Dir)


# efficiency: COB
figure(3)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s = fp[m].sum()
    eff_s = scd[m].sum()/s * 100.
    eff_ext = ecd[m].sum()/s * 100.
    eff_n = cd[m].sum()/s * 100.
    scatter([Pi],[eff_n],color='k')
    scatter([Pi],[eff_s],color='r')
    scatter([Pi],[eff_ext],color='b')

xlabel('P')
ylabel(r'Efficieny [%%]')
title(Dir)


# n_bad
figure(4)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    scatter([Pi],[data_nommask[m,7].sum()],color='k')
    scatter([Pi],[nbad[m].sum()],color='m')
    scatter([Pi],[data_nommask[m,21].sum()],color='g')
    scatter([Pi],[data_extmask[m,5].sum()],color='b')
    scatter([Pi],[nbad_sp[m].sum()],color='c')

xlabel('P')
ylabel(r'$n_{bad}$')
title(Dir)

# SPR_crit
figure(5)
clf()
plot(P[fp],data_nommask[fp,14],'k+')
plot(P[fp],data_2ndmask[fp,8],'r+')
plot(P[fp],data_extmask[fp,15],'b+')
xlabel('P')
ylabel(r'SPR$_{crit}$')
semilogy()
title(Dir)

# eta ratio
figure(6)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.) & fp
    scatter([Pi],[np.median(data_2ndmask[m,7] / data_nommask[m,13])],color='r')
    scatter([Pi],[np.median(data_extmask[m,10] / data_nommask[m,13])],color='b')
# plot(P[fp],data_2ndmask[fp,7] / data_nommask[fp,13],'r+')
# plot(P[fp],data_extmask[fp,10] / data_nommask[fp,13],'b+')
xlabel('P')
ylabel(r'$\eta/\eta_{T}$')
title(Dir)

# NSR ratio
figure(7)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.) & fp
    scatter([Pi],[np.median(data_2ndmask[m,4] / data_nommask[m,6])],color='r')
    scatter([Pi],[np.median(data_extmask[m,4] / data_nommask[m,6])],color='b')
# plot(P[fp],data_2ndmask[fp,4] / data_nommask[fp,6],'r+')
# plot(P[fp],data_extmask[fp,4] / data_nommask[fp,6],'b+')
xlabel('P')
ylabel(r'NSR$_{1h}$/NSR$_{1h,T}$')
title(Dir)

# NSR values
figure(8)
clf()
plot(P[fp],data_nommask[fp,6],'k+')
plot(P[fp],data_2ndmask[fp,4] ,'r+')
plot(P[fp],data_extmask[fp,4],'b+')
xlabel('P')
ylabel(r'NSR$_{1h}$')
semilogy()
title(Dir)

# eta: wrong formula
figure(9)
clf()
plot(P[fp],data_nommask[fp,13]*(1.-data_nommask[fp,9]),'k+')
plot(P[fp],data_2ndmask[fp,7]*(1.-data_2ndmask[fp,12]),'r+')
plot(P[fp],data_extmask[fp,10]*(1.-data_extmask[fp,7]) ,'b+')
xlabel('P')
ylabel(r'$\eta$')
semilogy()
title(r'$\eta$: wrong formula')
title(Dir)




# eta: correct formulas
figure(10)
clf()
plot(P[fp],data_nommask[fp,13] ,'k+')
plot(P[fp],data_2ndmask[fp,7] ,'r+')
plot(P[fp],data_extmask[fp,10]  ,'b+')
xlabel('P')
ylabel(r'$\eta$')
title(r'$\eta$: correct formula')
semilogy()
title(Dir)

# delta_obs
figure(11)
clf()
plot(P[fp],data_nommask[fp,14],'k+')
plot(P[fp],data_2ndmask[fp,8] ,'r+')
plot(P[fp],data_extmask[fp,11],'b+')
xlabel('P')
ylabel(r'$\delta_{obs}$')
semilogy()
title(Dir)

# delta_COB
figure(12)
clf()
plot(P[fp],data_nommask[fp,15],'k+')
plot(P[fp],data_2ndmask[fp,9] ,'r+')
plot(P[fp],data_extmask[fp,12],'b+')
xlabel('P')
ylabel(r'$\delta_{COB}$')
semilogy()
title(Dir)

# eta_COB
figure(13)
clf()
plot(P[fp],data_nommask[fp,17],'k+')
plot(P[fp],data_2ndmask[fp,11] ,'r+')
plot(P[fp],data_extmask[fp,14],'b+')
xlabel('P')
ylabel(r'$\eta_{COB}$')
semilogy()
title(Dir)

# data_nommask[:,7]: histogram of n_bad : nb of contaminant stars > SPR_crit (fixed delta_back value)
figure(14)
clf()
hist((data_nommask[:,7],nbad,nbad6c),range=[0,10],color=('m','k','g'),density=True,histtype='bar')

# histogram of n_bad : nb of contaminant stars eta>eta_min (for a real distribution in delta_back)

xlabel(r'$n_{bad}$')
title(Dir)


# histogram of n_bad_ext : nb of contaminant stars eta>eta_min (for a real distribution in delta_back)
figure(17)
clf()
hist(nbad_ext,range=[0,10],color='k',density=True)
title(Dir)

xlabel(r'$n_{bad}$')

# data_nommask[:,7]: histogram of n_bad : nb of contaminant stars > SPR_crit (fixed delta_back value)
figure(15)
clf()
hist((nbad_cob,nbad_ext_cob),range=[0,10],color=('k','b'),density=True,histtype='bar')

# histogram of n_bad : nb of contaminant stars eta>eta_min (for a real distribution in delta_back)

xlabel(r'$n_{bad,COB}$')
title(Dir)



# sigma_centroid
figure(16)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    scatter([Pi],[np.median(data_nommask[m,16])],color='k')


xlabel('P')
ylabel('COB error [pixel]')
title(Dir)

semilogy()


# efficiency: flux (all contaminants)
figure(17)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s =(eta_bt>flx_trh)[m,:].sum()
    eff_ext =  ( (eta_ext_bt>flx_trh) & (delta_obs_ext>delta_obs) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
    scatter([Pi],[eff_ext],color='g',marker='s')

    s = fp[m].sum()
    eff_s = sfd[m].sum()/s * 100.
    eff_ext = efd[m].sum()/s * 100.
    scatter([Pi],[eff_ext],color='b')
    scatter([Pi],[eff_s],color='r')

ylabel(r'Efficieny [%%]')
title(Dir)



# efficiency: COB (all contaminants)
figure(18)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s =(eta_bt>flx_trh)[m,:].sum()
    eff_s = ( (eta_cob_bt>cob_thr) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
    eff_ext =  ( (eta_cob_ext_bt>cob_thr) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
    scatter([Pi],[eff_s],color='k',marker='s')
    scatter([Pi],[eff_ext],color='b',marker='s')


ylabel(r'Efficieny [%%]')
title(Dir)
show()


# c=np.loadtxt('KeplerEclipsinBinaryCatalog_DR3_2019.csv',usecols=(3,4,1,5,6),delimiter=',')
# m=(c[:,0]>0) & (c[:,1]>0) & (c[:,3]>0) & (c[:,4]>0)
# p = m.sum()
#
# depth = np.append(c[m,0],c[m,1])*1e6
# td =  np.append(c[m,2]*c[m,3],c[m,2]*c[m,4])*86400./3600.
# data = np.zeros((p*2,2))
# data[:,0] = depth
# data[:,1] = td
# np.savetxt('KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt',data,header='Transit depth (pdeph and sdepth)  in ppm, transit duration in hours',fmt=("%8.3f %8.3f"))
