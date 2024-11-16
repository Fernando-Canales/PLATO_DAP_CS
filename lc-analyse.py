from pylab import *

# fname= 'sim-transit'
## fname= 'sim-transit-oem'  # optimal extended mask
## fname= 'sim-transit-true-oem'
## fname= 'sim-transit-4.5px'
# fname= 'sim-transit-1px' 
# fname= 'sim-transit-oem-m13'
fname = 'test'

bm = np.load(fname+'_bm.npy')
em = np.load(fname+'_em.npy')
IhrC = np.load(fname+'_IhrC.npy')
Ihr = np.load(fname+'_IhrT.npy')

param = np.load(fname+'_param.npz')

i0,j0 = int(param['i0']),int(param['j0'])
bg = float(param['bg']) 
x0C,y0C = param['x0C'],param['y0C']
x0,y0 = param['x0'],param['y0']

data = np.load(fname+'.npy')

# data[:,i]
# 0 : time
# 1 : dx 
# 2 : dy 
#
# 3 : Nominal flux
# 4 : Extended flux
# 5 : Nominal cob X
# 6 : Nominal cob Y
# 7 : Extented cob x
# 8 : Extented cob y
#
# 9-14: same for target only
# 15-20: same for contaminant star only
# 21-26:  same for contaminant star only, transit free
# 27-32: same for target only, transit free
# 33-39: Ftot, no transit

Dfx = data[:,4]-data[:,3] # flux diffrence Eflux - Nflux
#  transit free flux:
NFx0 = data[:,33] #
EFx0 = data[:,34] #
Dfx0 =  EFx0 - NFx0


# transit free COB
sf1 = data[:,27]/NFx0
sf2 = data[:,21]/NFx0
NCOBx0 = (data[:,29]*sf1 + data[:,23]*sf2)
NCOBy0 = (data[:,30]*sf2 + data[:,24]*sf2)

sf1 = data[:,28]/EFx0
sf2 = data[:,22]/EFx0
ECOBx0 = (data[:,31]*sf1 + data[:,25]*sf2)
ECOBy0 = (data[:,32]*sf2 + data[:,26]*sf2)

close('all')
set_cmap('binary')

figure(0)
clf()
imshow(bm,origin='lower',interpolation='none',extent=(0,6,0,6))
plot([x0-i0],[y0-j0],'bo',ms=10)
plot([x0C-i0],[y0C-j0],'ro',ms=10)
savefig(fname+'_0.png')

figure(1)
clf()
imshow(em,origin='lower',interpolation='none',extent=(0,6,0,6))
plot([x0-i0],[y0-j0],'bo',ms=10)
plot([x0C-i0],[y0C-j0],'ro',ms=10)
savefig(fname+'_1.png')

figure(2)
clf()
imshow(em-bm,origin='lower',interpolation='none',extent=(0,6,0,6))
plot([x0-i0],[y0-j0],'bo',ms=10)
plot([x0C-i0],[y0C-j0],'ro',ms=10)
savefig(fname+'_2.png')

set_cmap('hot')

figure(3)
clf()
imshow(Ihr+ IhrC,origin='lower',interpolation='none',extent=(0,6,0,6))


# FLUX
figure(10) # Normalised flux vs. time
clf()

plot((data[:,3]/data[0,3]-1.)*100.,'k') # N-flux
plot(( data[:,4]/data[0,4]-1.)*100.,'g') # E-flux
plot((Dfx/Dfx[0]-1.)*100.,'r')  # differential flux (e-flux minus n-flux)
xlabel('Time [s]')
ylabel(r'$\Delta Flux$ [%]')
title('non-detrended light curves')
savefig(fname+'_10.png')

figure(11)  # detrended flux vs. time
clf()
plot((data[:,3]/NFx0-1.)*1e2,'k') # N-flux
plot((data[:,4]/EFx0-1.)*1e2,'g') # E-flux
plot((Dfx/Dfx0-1.)*1E2,'r') # differential flux (e-flux minus n-flux)
title('detrended light curves')
xlabel('Time [s]')
ylabel(r'$\Delta Flux$ [%]')
savefig(fname+'_11.png')


# COB
figure(20) # COB variations, along X
clf()
plot((data[:,5]-data[0,5])*1e3,'k') # NCOB
plot((data[:,7]-data[0,7])*1e3,'g') # ECOB
#plot(data[:,6]-data[0,6],'k:')
#plot(data[:,8]-data[0,6],'r:')
xlabel('Time [s]')
ylabel(r'$\Delta COB$ [mpx]')
title('non-detrended COB time-series')
savefig(fname+'_20.png')

figure(21) # COB variations along X, drift removed
clf()
plot( (data[:,5]-data[0,5]-data[:,1])*1e3,'k') # NCOB
plot( (data[:,7]-data[0,7]-data[:,1])*1e3,'g') # ECOB
title('detrended COB time-series')
ylabel(r'$\Delta COB$ [mpx]')
xlabel('Time [s]')
savefig(fname+'_21.png')


figure(23)
clf()
plot(data[:,5]/NCOBx0,'k')
plot(data[:,7]/ECOBx0,'g')
savefig(fname+'_23.png')


# COB versus flux
figure(30)
clf()
plot( (data[:,3]/data[0,3]-1.)*100.,(data[:,5]-data[0,5]-data[:,1])*1e3,'k+') # Ncob vs nflux
plot( (data[:,4]/data[0,4]-1.)*100.,(data[:,7]-data[0,7]-data[:,1])*1e3,'g+') # ecob vs eflux
xlabel(r'$\Delta Flux$ [%]')
ylabel(r'$\Delta COB$ [mpx]')
savefig(fname+'_30.png')

show()
