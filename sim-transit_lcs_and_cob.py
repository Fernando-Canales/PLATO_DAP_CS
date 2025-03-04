"""
This script plots the light curves
and COB of the extended and nominal masks

Fernando, 25.Feb.2025
"""
from pylab import * # type: ignore
import h5py # type: ignore
import spline2dbase # type: ignore
import math
import numpy as np #type:ignore
import matplotlib.pyplot as plt # type:ignore

fpsf = 'PSF_Focus_0mu_0.2pxdif.npz'
psfidx = 100 # index of the PSF to be used
subres = 64 # psf sub-pixel resolution
bsres = 20 # b-spline resolution
PSFsizex = 8 # PSF size in pixels
PSFsizey = 8
colormap = 'Pastel2'
resultsDIR = '/home/fercho/double-aperture-photometry/plots_pdfs/'

contaminants = [
    (4.5, 4.5),
    (1.5, 4.5),
    (0.25, 0.75),
    (1.0, 1.5),
    (2.5, 1.5)
]

nexp = 1500 #  number of exposures
driftamplitude = 0 # amplitude of the drif in pix/90 days

ron = 52.  # readout noise [e-]
zero_point = 20.62 # camera zero point
integration_time = 21. # integration time [s]
fl = 247.5 # focal length [mm]
gain = 25. # full gain [e-/ADU]
Nc = 24 # number of cameras
addnoise = True # add random noise or not (photon noise + detector noise)
bg = 2500. # background level [e-/pix]
bgerr = 100. # background residual error [e-/pix]

Mpx = 6 # size of the imagette
Mag = 11 # PLATO magnitude , target
MagC = Mag + 3 # magnitude contaminant
DistC = 1.6 # distance of the contaminant from the target
x0,y0 = 3, 3 # target position in the imagette


# stellar variability
ffluxvar = 'fluxtransitcorot7b_90d.hdf5' # Corot 7b depth=3.958e-4
fluxvarsf = 0.5/(1.958e-4) #  scaling factor for the false transit (transit in the contaminant)
fluxvarsftf = 0. # scaling factor for the true transit (transit in the target)


# fname= 'sim-transit-oem' #  distance 1pix MagC=16, OEM
# fname= 'sim-transit' #  distance 1pix MagC=16
# fname= 'sim-transit-1px' # distance 1pix MagC=13
# fname= 'sim-transit-4.5px' # distance 4.5pix MagC=13
# fname= 'sim-transit-oem-m13' # distance 1pix MagC=13, OEM
fname= 'test' # name of the output file


def add_poisson_noise(image):
    '''
    image_new =  add_poisson_noise(image)
    add poisson noise to an input image

    :param image: the input image
    :return: an image with poisson noise added
    '''
    image_new = np.random.poisson(lam=image.flatten(), size=(image.size))
    image_new = np.array(image_new.reshape(image.shape),dtype=image.dtype)
    return image_new

def add_gauss_noise(image,sigma):
    '''
    image_new = add_gauss_noise(image,sigma)
    add Gaussian noise to an input image

    :param image: the input image
    :param sigma: standard deviation associated with the Gaussian distribution
    :return: an image with Gaussian noise added
    '''
    noise = np.random.normal(size=image.size,scale=sigma)
    return image + noise.reshape(image.shape)


def barycenter(array,mask=None,x=None,y=None,subres=1):
    if(isinstance(x, (np.ndarray, np.generic) )== False):
        x=(np.arange(0,array.shape[1])+0.5)/float(subres)
    if(isinstance(y, (np.ndarray, np.generic) )== False):
        y=(np.arange(0,array.shape[0])+0.5)/float(subres)

    if(x.ndim ==1 & y.ndim ==1):
        x,y = np.meshgrid(x,y)

    if(mask is not None):
        weight=np.sum(array*mask)
        bx= np.sum(array*x*mask)/weight
        by= np.sum(array*y*mask)/weight
    else:
        tmp=np.sum(array)
        bx= np.sum(array*x)/tmp
        by= np.sum(array*y)/tmp
    return bx,by

def NSR(M,sb, sd, sq, It, Ic):
    nsr = np.sqrt(np.sum((It + Ic + sb + sd ** 2 + sq ** 2)*M))/ np.sum(It*M)
    return nsr


def aperture_computation(ft, fc, sb, sd, sq):
    """ We are following the procedure described
    in subsection 4.6.3. of Marchiori paper
    (Binary mask)

    Args:
        ft (_float_): target flux
        fc (_float_): contaminant flux
        sb (_float_): background noise
        sd (_float_): detector noise
        sq (_float_): quantization noise
    """
    # First we compute the NSR of the system. Eq. (36) of Marchiori's paper
    nsr = np.sqrt(ft + fc + sb + sd ** 2 + sq ** 2) / ft

    # Then we flatten the nsr and also the fluxes
    nsr_1d = nsr.flatten()
    ft_1d = ft.flatten()
    fc_1d = fc.flatten()

    # Then we sort the 1-D nsr in increasing order and obtain the index of the elements of the array before sorting them
    nsr_1d_index = np.argsort(nsr_1d)

    # Then we obtain the target and contaminant flux for such indexes
    ft_1d = ft_1d[nsr_1d_index]
    fc_1d = fc_1d[nsr_1d_index]

    # Then we compute the aggregat noise-to-signal ratio. Eq. (37) in Marchiori's paper
    nsr_agg = np.zeros(len(ft_1d))
    for i in range(1, len(ft_1d) + 1):
        nsr_agg[i - 1] = np.sqrt(np.sum(ft_1d[:i] + fc_1d[:i] + sb + sd ** 2 + sq ** 2)) / np.sum(ft_1d[:i])

    # We create now a vector full with zeros
    aperture = np.zeros(36)

    # Then we create a boolean mask
    boolean_mask_for_the_aperture = nsr_1d_index[:np.argmin(nsr_agg) + 1]

    # Then we obtain the nominal mask
    aperture[boolean_mask_for_the_aperture] = 1

    # Then we reshape the aperture
    aperture = aperture.reshape((6, 6))

    return aperture, nsr_agg


# We define here a function for obtaining an extended mask given a nominal-binary mask
def extended_binary_mask(mask, W):
    ny, nx = mask.shape
    maske = np.zeros((ny, nx))
    maske[:, :] = mask[:, :]
    for j in range(ny):
        for i in range(nx):
            if mask[j, i] > 1.0 - 1e-5:
                for k in range(-W, W + 1):
                    if (j + k >= 0) & (j + k < ny):
                        for m in range(-W, W + 1):
                            if (i + m >= 0) & (i + m < nx):
                                maske[j + k, i + m] = 1
    return maske

psfdata = np.load(fpsf)

psfbs = np.array(psfdata['psfbs'][psfidx])
pxc = float(psfdata['pxc'][psfidx])
pyc = float(psfdata['pyc'][psfidx])
xpsf_pix = float(psfdata['xpsf_pix'][psfidx])
ypsf_pix = float(psfdata['ypsf_pix'][psfidx])

xpsf = xpsf_pix*18e-3 # in mm
ypsf = ypsf_pix*18e-3 # in mm

r = math.sqrt(xpsf**2+ypsf**2)
angradius = math.atan2(r,fl) # angular radius [rad]
print('PSF angular radius = %f [deg]' % (angradius*180./math.pi))


# stellar variability (here a transità
fluxvar = h5py.File(ffluxvar)['Variation']

# postion of the contaminant in the imagette
x0C, y0C = x0+DistC/np.sqrt(2.),y0+DistC/np.sqrt(2.)

# reference flux (at P=0) in e-/exposure
flux_m0_ref = 10.**(0.4*zero_point)*integration_time*(math.cos(angradius))**2

fxT = flux_m0_ref*10.**( -(Mag/2.5) ) # target flux [e-/exposure]
fxC = flux_m0_ref*10.**( -(MagC/2.5) ) # contaminant flux [e-/exposure]

print('target flux = %f [ke-] = %f [ke-/s] (single camera)' % (fxT*1e-3,fxT*1e-3/integration_time))

sizex = Mpx
sizey = Mpx

# imagette obtained by integrating the b-spline decomposition of the PSF
IT = spline2dbase.Spline2Imagette(psfbs,bsres,sizex,sizey,offx=x0-pxc,offy=y0-pyc)
IC = spline2dbase.Spline2Imagette(psfbs,bsres,sizex,sizey,offx=x0C-pxc,offy=y0C-pyc)
IT *= fxT
IC *= fxC

IhrT = spline2dbase.Spline2Imagette(psfbs,bsres,sizex,sizey,offx=x0-pxc,offy=y0-pyc,subres=subres)
IhrC = spline2dbase.Spline2Imagette(psfbs,bsres,sizex,sizey,offx=x0C-pxc,offy=y0C-pyc,subres=subres)
IhrT *= fxT
IhrC *= fxC


#data = np.zeros((nexp,3+6*6))
data = np.zeros((nexp, 40))
i0 = int(round(x0-Mpx/2.))
j0 = int(round(y0-Mpx/2.))

DxC = np.array([x0C-x0])
DyC = np.array([y0C-y0])
MagC = np.array([MagC])

# calculation of the nominal mask (optimal binary mask)
sq = gain/math.sqrt(12.) # quantification noise
bm = aperture_computation(IT, IC, bg, ron,sq)[0]
NSRn = NSR(bm,bg,ron, sq, IT, IC)
print('SNRn = %f' % (1./NSRn))

# calculation of the corresponding extended mask
em = extended_binary_mask(bm,1)
NSRe = NSR(em,bg,ron, sq, IT, IC)
print('SNRe = %f' % (1./NSRe))

# calculation of the corresponding secondary mask
sm = aperture_computation(IC, IT, bg, ron, sq)[0]
NSRs = NSR(sm, bg, ron, sq, IC, IT)
print('SNRs = %f' % (1./NSRs))


MagTC = Mag + 2.5*(math.log10(fxT) - math.log10(fxC+fxT) )

print('Target coordinates: %f' % x0, y0)
print('Contaminant coordinates: %f' % x0C, y0C)


for t in range(nexp):
    dx = t*25./(86400.*90)* driftamplitude /math.sqrt(2.)
    dy = dx

    IT0 = spline2dbase.Spline2Imagette(psfbs, bsres, sizex, sizey, offx=x0+dx - pxc, offy=y0+dy - pyc)
    IC0 = spline2dbase.Spline2Imagette(psfbs, bsres, sizex, sizey, offx=x0C+dx - pxc, offy=y0C+dy - pyc)

    IT0 *= fxT
    IC0 *= fxC

    Itot0 = IT0 + IC0 + bgerr # Transit free

    # adding the stellar variability
    IT = IT0 * (1. + (-1. + fluxvar[t, 1]) * fluxvarsftf)  # with a transit in the target
    IC = IC0* (1. + (-1. + fluxvar[t, 1]) * fluxvarsf) # with a transit in the contaminant

    Itot = IT + IC + bg # with Transit
    if(addnoise):
        Itot = add_poisson_noise(Itot*Nc)/Nc
        Itot = add_gauss_noise(Itot,math.sqrt( (ron**2 + sq**2)/Nc))
    Itot = Itot -bg + bgerr

    data[t,0] = t*25.
    data[t,1] = dx
    data[t,2] = dy
    s = 3
    data[t,s] = np.sum(Itot*bm)
    data[t,1+s] = np.sum(Itot*em)
    data[t,2+s:4+s] = barycenter(Itot,mask=bm)
    data[t,4+s:6+s] = barycenter(Itot,mask=em)
    #data[t,4+s:6+s] = barycenter(Itot,mask=sm)
    #data[t, 7 + s] = np.sum(Itot * sm)  # Secondary flux for target
    
    s = 9
    data[t,s] = np.sum(IT*bm)
    data[t,1+s] = np.sum(IC*em)
    data[t,2+s:4+s] = barycenter(IT,mask=bm)
    data[t,4+s:6+s] = barycenter(IT,mask=em)
    #data[t,4+s:6+s] = barycenter(IC,mask=sm)
    #data[t, 7 + s + 6] = np.sum(IC * sm)  # Secondary flux for contaminant

    s = 15
    data[t,s] = np.sum(IC*bm)
    data[t,1+s] = np.sum(IT*em)
    data[t,2+s:4+s] = barycenter(IC,mask=bm)
    data[t,4+s:6+s] = barycenter(IC,mask=em)
    #data[t,4+s:6+s] = barycenter(IT,mask=sm)
    #data[t, 7 + s + 12] = np.sum(IT * sm)  # Secondary flux for contaminant

    s = 21
    data[t,s] = np.sum(IC0*bm)
    data[t,1+s] = np.sum(IT0*em)
    data[t,2+s:4+s] = barycenter(IC0,mask=bm)
    data[t,4+s:6+s] = barycenter(IC0,mask=em)
    #data[t,4+s:6+s] = barycenter(IT0,mask=sm)
    #data[t,7+s] = np.sum(IT0*sm)

    s = 27
    data[t,s] = np.sum(IT0*bm)
    data[t,1+s] = np.sum(IC0*em)
    data[t,2+s:4+s] = barycenter(IT0,mask=bm)
    data[t,4+s:6+s] = barycenter(IT0,mask=em)
    #data[t,4+s:6+s] = barycenter(IC0,mask=sm)
    #data[t,7+s] = np.sum(IC0*sm)

    s = 33
    data[t,s] = np.sum(Itot0*bm)
    data[t,1+s] = np.sum(Itot0*em)
    data[t,2+s:4+s] = barycenter(Itot0,mask=bm)
    data[t,4+s:6+s] = barycenter(Itot0,mask=em)
    #data[t,4+s:6+s] = barycenter(Itot0,mask=sm)
    #data[t,7+s] = np.sum(Itot0*sm)

# saving data
np.save(fname+'.npy',data)
np.save(fname+'_bm.npy',bm)
np.save(fname+'_em.npy',em)
np.save(fname+'_sm.npy',sm)
np.save(fname+'_IhrC.npy',IhrC)
np.save(fname+'_IhrT.npy',IhrT)
np.savez(fname+'_param',x0=x0,x0C=x0C,y0=y0,y0C=y0C,i0=0,j0=0,bg=bg)

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
# 21-26: same for contaminant star only, no transit
# 27-32: same for target only, no transit
# 33-40: Ftot, no transit

centroid_x_nom = data[:, 5]
centroid_y_nom = data[:, 6]
centroid_x_ext = data[:, 7]
centroid_y_ext = data[:, 8]

# Compute radial shifts (or adapt this formula to however you define the shift).
centroid_shift_nom = np.sqrt((centroid_x_nom/data[0,5])**2 + (centroid_y_nom/data[0,6])**2)
centroid_shift_ext = np.sqrt((centroid_x_ext/data[0, 7])**2 + (centroid_y_ext/data[0,8])**2)

max_cob = max(centroid_shift_nom.max(), centroid_shift_ext.max())
# Function to add dashed lines for the nominal mask
def plot_nominal_mask_contour(nominal_mask, plot_color):
    for j in range(6):
        for i in range(6):
            if nominal_mask[j, i] == 1:
                plt.plot([i, i+1], [j, j], color=plot_color, linestyle='--', lw=3)  # Bottom edge
                plt.plot([i, i+1], [j+1, j+1], color=plot_color, linestyle='--', lw=3)  # Top edge
                plt.plot([i, i], [j, j+1], color=plot_color, linestyle='--', lw=3)  # Left edge
                plt.plot([i+1, i+1], [j, j+1], color=plot_color, linestyle='--', lw=3)  # Right edge

plt.figure(1, figsize=(7, 6), dpi=300)
#plt.scatter(data[:, 0]/ 3600., data[:, 3]/np.mean(data[0:99, 3]), label= 'Nominal Flux', color='black',  marker='o', s=3)
#plt.scatter(data[:, 0]/ 3600., data[:, 4]/np.mean(data[0:99, 4]), label='Extended Flux', color='brown',  marker='P', s=3)
plt.scatter(data[:, 0]/3600., centroid_shift_nom - max_cob, label='Nominal COB', color='black', marker='o', s=3)
plt.scatter(data[:, 0]/3600., centroid_shift_ext - max_cob, label='Extended COB', color='blue', marker='o', s=3) 
plt.xlabel('Time [hours]', fontsize=14)
plt.ylabel('Centroid shift [pixel]', fontsize=14)
plt.legend()
# Adjust the layout manually to increase the space at the bottom for the x-axis label
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)  # Increase bottom margin
plt.savefig(resultsDIR+'Secondary_and_Nominal_COBs.png', dpi=300, format='png')
plt.show()


plt.figure(2, figsize=(7,6), dpi=300)
NSR_agg = aperture_computation(IT, IC, bg, ron,sq)[1]
plt.plot(13333*NSR_agg, 'o-')
plt.axvline(x=4, color='orange', linestyle='--')
plt.xlabel('Pixels')
plt.ylabel(r'$\mathrm{NSR}_{\mathrm{agg}}$ over 1 hour and 24 cam. [ppm $\mathrm{hr}^{\frac{1}{2}}$]')
plt.savefig(resultsDIR+'NSR_agg_sim_transit.png', dpi=300, format='png')
plt.show()

# Plotting the masks
plt.imshow(bm, origin='lower', extent=(0,6,0,6), cmap=colormap)
plt.grid(True, linewidth=2)
plt.scatter(x0, y0, s=80, c='green', zorder=5)  # 'c' sets the face color of the marker
#plt.scatter(x0C, y0C, s=80, c='cyan', marker='^', zorder=5)  # 'c' sets the face color of the marker
x_c = [c[0] for c in contaminants]
y_c = [c[1] for c in contaminants]
#plt.scatter(x_c, y_c, s=80, c='cyan', marker='^', zorder=5)
plt.savefig(resultsDIR+'nominal_mask_example_target_and_contaminant.png', dpi=300, format='png')
plt.show()
plt.imshow(em, origin='lower', extent=(0,6,0,6), cmap=colormap)
plt.grid(True, linewidth=2)
plt.scatter(x0, y0, s=80, c='green', zorder=5)  # 'c' sets the face color of the marker
plt.scatter(x0C, y0C, s=80, c='cyan', marker='^',  zorder=5)  # 'c' sets the face color of the marker
plot_nominal_mask_contour(nominal_mask=bm, plot_color='black')
plt.savefig(resultsDIR+'extended_mask_example_target_and_contaminant.png',  dpi=300, format='png')
plt.show()
plt.imshow(sm, origin='lower', extent=(0,6,0,6), cmap=colormap)
plt.grid(True, linewidth=2)
plt.scatter(x0, y0, s=80, c='green', zorder=5)  # 'c' sets the face color of the marker
plt.scatter(x0C, y0C, s=80, c='cyan', marker='^', zorder=5)  # 'c' sets the face color of the marker
plot_nominal_mask_contour(nominal_mask=bm, plot_color='white')
plt.savefig(resultsDIR+'secondary_mask_example_target_and_contaminant.png',  dpi=300, format='png')
plt.show()