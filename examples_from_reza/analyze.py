from pylab import * # type: ignore
from pylab import * # type: ignore
from pylab import * # type: ignore
from imagette import ran_unique_int
import math


# Dir = '40000/'
# Dir = 'PSF_Focus_0mu_0.2pxdif/'
Dir = 'test40/'
# test4: new dback and td such as to have gamma=0.5
# test5:  extended mask extension by 2 pixels
# test6: optimal extended mask built by joining together (union) the secondary mask of each contaminant star that can create a FP
# test7: v1 optimal extended mask: for each contaminant, take the extra pixels that maximize the significance
# CatalogueDIR='/home/fercho/double-aperture-photometry/catalogues_stars/'
CatalogueDIR= './'
Pmin =  10.
Pmax =  13
binsize = 0.5
cob_thr= 3.
ntr = 3  # number of transits
flx_trh =  7.1
ext_flx_trh = 3. # 7.1
sec_flx_trh = 3. # 7.1
dback_ref =  132000.
td_ref = 6.72*0.46**2
depth_sig_scaling = 3. #

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
# 52-61: 10 first eta_COB
# 62-71: 10 first eta_10first
# 72-81: IDs of the 10 first contaminants
# 82-91: delta x
# 92-101: delta y
# 102: x target position in the imagette
# 103: Y target position in the imagette
# 104:  PSF index for the contaminant stars
# 105: theoretical n-mask NSR

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
# 13: Gamma_c
# 14-23: 10 first SPRtot_c values
# 24-33: 10 first Gamma_c values
# 34-43: 10 first delta_COB_sig_1h_24c
# 44-53: 10 first NSR1h_c
# 54-63: 10 first w_c_key
# 64-73: 10 first eta_COB_c
# 74-83: 10 first eta_c_10first
# 84: theoretical s-mask NSR


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
# 19: Gamma_ext
# 20: NSR_chi
# 21-30: 10 first SPRk values
# 31-40: 10 first Gamma values
# 41-50: 10 first delta_COB_sig_1h_24c
# 51-60: 10 first eta_COB_ext
# 61-70: 10 first eta_ext_10first
# 71-80: 10 first eta_dtd  (significance of the differential transit depth)

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


u = (data_2ndmask[:,6]<=1)
data_2ndmask[u,9] = 0. # delta_COB_c
data_2ndmask[u,11] = 0. # eta_cob_c

#
# data_nommask[:,13] *= (dback_ref/85000.)*sqrt(td_ref/4.)
# data_2ndmask[:,7]  *= (dback_ref/85000.)*sqrt(td_ref/4.)
# data_2ndmask[:,11]  *= (dback_ref/85000.)*sqrt(td_ref/4.)

sig_depth_s_sec = data_2ndmask[:,4]*(1.-data_2ndmask[:,12])/sqrt(td_ref*ntr) # NSR*(1-SPRtot)
sig_depth_nom =   data_nommask[:,6]*(1. - data_nommask[:,9]) /sqrt(td_ref*ntr)  # NSR*(1-SPRtot)
sig_depth_s = np.sqrt(sig_depth_nom**2+ sig_depth_s_sec**2)

sig_depth_ext = data_extmask[:, 4] * (1. - data_extmask[:, 7]) /sqrt(td_ref*ntr)  # NSR*(1-SPRtot)


sig_depth = np.sqrt(sig_depth_nom**2+ sig_depth_ext**2)
fp = (data_nommask[:,13] > flx_trh) # eta_t>7: false positive (not yet identified as such)
sfd = fp & (data_2ndmask[:,7] > sec_flx_trh)  &  (data_2ndmask[:,8] > data_nommask[:,14] + depth_sig_scaling*sig_depth_s)  # eta_c>7.1 and delta_obt_c > delta_obs_t : secondary mask detection
efd = fp & (data_extmask[:,10] > ext_flx_trh)   & (data_extmask[:,15] > data_nommask[:,14]+ depth_sig_scaling*sig_depth)  # eta_ext>7.1 and delta_obt_ext > delta_obs_t : extended flux detection

cd = fp & (data_nommask[:,17] > cob_thr) # eta_cob>3: COB detection
scd = fp & (data_2ndmask[:,11] > cob_thr)  & (data_2ndmask[:,9]/10.  > data_2ndmask[:,10] )  # eta_cob_c>3: s-COB detection  delta_COB_c>delta_COB_sig_1h_24c_c/10.
ecd = fp & (data_extmask[:,14] > cob_thr)  # eta_cob_ext>3: e-COB detection


n = data_nommask.shape[0]

SPRk = data_nommask[:,22:32]
# SPRk = data_extmask[:,22:32]
dback_set = np.loadtxt(CatalogueDIR + 'KeplerEclipsinBinaryCatalog_DR3_2019_depth.txt')
dback_n = dback_set.shape[0]
seed = 23434434
np.random.seed(seed)
nbad = np.zeros(n)
nbad_cob = np.zeros(n)

nbad_sp = np.zeros(n) # small planet 2<R_Ee -> 4*84ppm

nbad6c = np.zeros(n)

eta_ext_cob_prueba = data_extmask[:, 51:61]
delta_obs = np.zeros((n,10))
delta_obs_ext = np.zeros((n,10))
eta_bt = np.zeros((n,10))
eta_cob_bt = np.zeros((n,10))
eta_cob_ext_bt = np.zeros((n,10))
eta_ext_bt = np.zeros((n,10))
ntr = 3  # number of transits
nbad_ext = np.zeros(n)
nbad_ext_cob = np.zeros(n)


spr_t_ext = data_extmask[:,4]*(1. - data_extmask[:,7])  /(dback_ref*math.sqrt(td_ref*ntr))
spr_t_ext = np.repeat(np.reshape(spr_t_ext,(n,1)),10 ,axis=1)

spr_t_nom = data_nommask[:,6]*(1. - data_nommask[:,9])  /(dback_ref*math.sqrt(td_ref*ntr))
spr_t_nom = np.repeat(np.reshape(spr_t_nom,(n,1)),10 ,axis=1)

spr_t = np.sqrt(spr_t_nom**2 + spr_t_ext**2)



sig_depth_ext = np.reshape(sig_depth_ext,(n,1))
sig_depth_ext = np.repeat(sig_depth_ext,10 ,axis=1)

sig_depth_nom = np.reshape(sig_depth_nom,(n,1))
sig_depth_nom = np.repeat(sig_depth_nom,10 ,axis=1)

sig_depth = np.sqrt(sig_depth_nom**2+ sig_depth_ext**2)


# sig_depth_ext = np.zeros((n,10))
# sig_depth_nom = np.zeros((n,10))
# compute n_bad  for a real distribution in delta_back: nb of contaminant stars eta>eta_min=7.1
for i in range(n):
    # j = ran_unique_int(10,interval=[0,dback_n-1]) # random sort of a BT (background transit)
    # dback = dback_set[j,0] # transit depth
    # td = dback_set[j,1] # transit duration
    dback = np.ones(10)*dback_ref
    td = np.ones(10)*td_ref
    eta_bt[i,:] = (SPRk[i,:]/data_nommask[i,8])*flx_trh *(dback/dback_ref)*sqrt(td/td_ref) # significance of the BT in the nominal flux
    eta_ext_bt[i,:] = dback*data_extmask[i,21:31]*np.sqrt(td*ntr)/(1-data_extmask[i,7])/data_extmask[i,4] # significance in the extended mask

    nbad[i] = np.sum(eta_bt[i,:]>flx_trh) # number of false detection in the nominal mask
    nbad_ext[i] = np.sum(eta_ext_bt[i,:]>ext_flx_trh)  # number of false detection in the extended mask
    delta_obs[i,:] = dback*SPRk[i,:] # observed transit depth

    delta_int = delta_obs[i,:]/(1. -data_nommask[i,9] ) # inferred intrinsic transit depth
    delta_obs_ext[i,:] = dback*data_extmask[i,21:31] # observed transit depth
    nbad_sp[i] = np.sum( (eta_bt[i,:]>flx_trh) & (delta_int<4*84. )) # for small planets
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

    all = np.append(data_nommask[m,4],data_2ndmask[m,5])
    all = np.append(all,np.array(data_extmask[m,2]))
    scatter([Pi],[len(np.unique(np.array(all,dtype=np.int64)))],color='g')

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
figure(15)
clf()
hist(nbad_ext,range=[0,10],color='k',density=True)
title(Dir)

xlabel(r'$n_{bad}$')

# data_nommask[:,7]: histogram of n_bad : nb of contaminant stars > SPR_crit (fixed delta_back value)
figure(16)
clf()
hist((nbad_cob,nbad_ext_cob),range=[0,10],color=('k','b'),density=True,histtype='bar')

# histogram of n_bad : nb of contaminant stars eta>eta_min (for a real distribution in delta_back)

xlabel(r'$n_{bad,COB}$')
title(Dir)



# sigma_centroid
figure(17)
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
figure(18)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s =(eta_bt>flx_trh)[m,:].sum()
    eff_ext =  ( (eta_ext_bt>ext_flx_trh) & (delta_obs_ext>delta_obs+depth_sig_scaling*sig_depth) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
    scatter([Pi],[eff_ext],color='b',marker='s')



    s = fp[m].sum()
    eff_s = sfd[m].sum()/s * 100.
    eff_ext = efd[m].sum()/s * 100. # most prominent contaminant star
    scatter([Pi],[eff_ext],color='g')
    scatter([Pi],[eff_s],color='r')

ylabel(r'Efficieny [%%]')
title(Dir)




figure(19)
clf()
for i in range(nP):
    Pi = Pmin + i*binsize
    m = (P >= Pi -binsize/2.) & (P < Pi + binsize/2.)
    s =(eta_bt>flx_trh)[m,:].sum()
    eff_s = scd[m].sum()/s * 100.
    eff_nom = ( (eta_cob_bt>cob_thr) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
    eff_ext =  ( (eta_cob_ext_bt>cob_thr) & (eta_bt>flx_trh))[m,:].sum()/s * 100.
#    eff_ext_prueba = ((eta_ext_cob_prueba>cob_thr) & (eta_bt > flx_trh))[m,:].sum()/s * 100
    scatter([Pi], [eff_s], color='r', marker='s',label='SCOB')
    scatter([Pi], [eff_nom], color='k',marker='s',label='NCOB')
    scatter([Pi], [eff_ext], color='b',marker='s',label='ECOB')
#    scatter([Pi], [eff_ext_prueba], color='g', marker='s')
ylabel(r'Efficiency [%%]')
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


s = (eta_bt>flx_trh)
a = (eta_ext_bt>ext_flx_trh) & (data_extmask[:,21:31]>data_nommask[:,22:32]+depth_sig_scaling*spr_t)
eff_ext = ( a  & s).sum()/s.sum()*100.
print('extended flux efficiency: %f' % eff_ext)


r =(1.-data_extmask[:,7]) / (1. - data_nommask[:,9] )# (1-SPR_to^EXT)/(1-SPR_tot^nom)
r = np.reshape(r,(r.shape[0],1))
r = np.repeat(r,10,axis=1)
a = (eta_ext_bt>ext_flx_trh) & (data_extmask[:,21:31]>r*data_nommask[:,22:32]+depth_sig_scaling*spr_t)
eff_ext = ( a  & s).sum()/s.sum()*100.
print('extended flux efficiency - planet hypothesis: %f' % eff_ext)


a_noisefree  = (eta_ext_bt>ext_flx_trh) & (delta_obs_ext>delta_obs)
eff_ext_noisefree = ( a_noisefree   & s).sum()/s.sum()*100.
print('extended flux efficiency with noise free condition on delta: %f' % eff_ext_noisefree)


s_max = (eta_bt[:,0]>flx_trh)
a_max = (eta_ext_bt[:,0]>ext_flx_trh) & (delta_obs_ext[:,0]>delta_obs[:,0]+depth_sig_scaling*sig_depth[:,0])
eff_ext_max = ( a_max  & s_max).sum()/s_max.sum()*100.
print('extended flux efficiency most significant cont: %f' % eff_ext_max)


n_ext_det = np.sum(a & s,axis=1)


eff_s = sfd.sum() / fp.sum() * 100.
print('secondary mask efficiency: %f' % eff_s)


b = (eta_cob_bt>cob_thr)
# efficiency: COB (all contaminants)
eff_cob_nom = ( b & s).sum()/s.sum()*100.
print('nominal cob efficiency: %f' % eff_cob_nom)

n_cob_det = np.sum(b & s,axis=1)


c = (eta_cob_ext_bt>cob_thr)
eff_cob_ext = ( c & s).sum()/s.sum()*100.
print('extended cob efficiency: %f' % eff_cob_ext)

n_cob_ext_det = np.sum(c & s,axis=1)


eff_cob_s = ( scd & fp).sum()/fp.sum()*100.
print('secondary mask cob efficiency: %f' % eff_cob_s)


N = s.sum()
f = ( (a==False)  & c & s).sum()/N
print('fraction only detected by the ECOB but not by the EFX %f' % f)
f = ( (a)  & (c==False) & s).sum()/N
print('fraction only detected by the EFX but not by the ECOB: %f' % f)

f = ( (a==False)  & b & s).sum()/N
print('fraction only detected by the NCOB but not the EFX %f' % f)
f = ( (a)  & (b==False) & s).sum()/N
print('fraction only detected by the  EFX but not by NCOBB: %f' % f)


f = ( (c)  & (b==False) & s).sum()/N
print('fraction only detected by the  ECOB but not by the NCOB: %f' % f)

# metric assigment algo

metric_priority = -np.ones(n,dtype=np.int32)
nfp_detected = np.zeros(n,dtype=np.int32)
for i in range(n):
    nfp = s[i].sum()
    if( nfp ==0): # no FP
        metric_priority[i] = 0
    elif (nfp ==1): # only one FP, SFX or NCOB
        if(data_2ndmask[i,7] > sec_flx_trh):
            metric_priority[i] = 1 # SFX
            nfp_detected[i] = 1
        else:
            metric_priority[i] = 3 # NCOB
            ncob = np.sum((eta_cob_bt[i] > cob_thr) & s[i])
            nfp_detected[i] = ncob
    else: # 2 or more FP
        nefx = np.sum((eta_ext_bt[i] > ext_flx_trh) & (data_extmask[i,21:31] >data_nommask[i,22:32] + depth_sig_scaling * spr_t[i]) & s[i])
        ncob = np.sum((eta_cob_bt[i]>cob_thr) & s[i])
        necob = np.sum((eta_cob_ext_bt[i] > cob_thr)  & s[i])
        if(nefx >= ncob):
            metric_priority[i] = 2  # EFX
            nfp_detected[i] = nefx
        elif(ncob>=necob):
            metric_priority[i] = 3  # NCOB
            nfp_detected[i] = ncob
        else:
            metric_priority[i] = 4  # ECOB
            nfp_detected[i] = necob

figure(20)
clf()
h,_,_ = hist(metric_priority,density=True,range=(-0.5,4.5),bins=5)
print(h)
