from pylab import *
import h5py
import spline2dbase
import math

fpsf = 'PSF_Focus_0mu_0.2pxdif.npz'
psfidx = 100 # index of the PSF to be used
subres = 64 # psf sub-pixel resolution
bsres = 20 # b-spline resolution
PSFsizex = 8 # PSF size in pixels
PSFsizey = 8

nexp = 1500 #  number of exposures
driftamplitude = 1. # amplitude of the drif in pix/90 days

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
MagC = Mag + 2 # magnitude contaminant
DistC = 1. # distance of the contaminant from the target
x0,y0 = 3.1, 2.9 # target position in the imagette

# stellar variability
ffluxvar = 'fluxtransitcorot7b_90d.hdf5' # Corot 7b depth=3.958e-4
fluxvarsf = 0.5/(3.958e-4) #  scaling factor for the false transit (transit in the contaminant)
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

    return aperture


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
x0C, y0C = x0+DistC/sqrt(2.),y0+DistC/sqrt(2.)

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

data = np.zeros((nexp,3+6*6))

i0 = int(round(x0-Mpx/2.))
j0 = int(round(y0-Mpx/2.))

DxC = np.array([x0C-x0])
DyC = np.array([y0C-y0])
MagC = np.array([MagC])

# calculation of the nominal mask (optimal binary mask)
sq = gain/math.sqrt(12.) # quantification noise
bm = aperture_computation(IT, IC, bg, ron,sq)
NSRn = NSR(bm,bg,ron, sq, IT, IC)
print('SNRn = %f' % (1./NSRn))

# calculation of the corresponding extended mask
em = extended_binary_mask(bm,1)
NSRe = NSR(em,bg,ron, sq, IT, IC)
print('SNRe = %f' % (1./NSRe))

MagTC = Mag + 2.5*(math.log10(fxT) - math.log10(fxC+fxT) )

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
    
    s = 9
    data[t,s] = np.sum(IT*bm)
    data[t,1+s] = np.sum(IT*em)
    data[t,2+s:4+s] = barycenter(IT,mask=bm)
    data[t,4+s:6+s] = barycenter(IT,mask=em)

    s = 15
    data[t,s] = np.sum(IC*bm)
    data[t,1+s] = np.sum(IC*em)
    data[t,2+s:4+s] = barycenter(IC,mask=bm)
    data[t,4+s:6+s] = barycenter(IC,mask=em)

    s = 21
    data[t,s] = np.sum(IC0*bm)
    data[t,1+s] = np.sum(IC0*em)
    data[t,2+s:4+s] = barycenter(IC0,mask=bm)
    data[t,4+s:6+s] = barycenter(IC0,mask=em)

    s = 27
    data[t,s] = np.sum(IT0*bm)
    data[t,1+s] = np.sum(IT0*em)
    data[t,2+s:4+s] = barycenter(IT0,mask=bm)
    data[t,4+s:6+s] = barycenter(IT0,mask=em)

    s = 33
    data[t,s] = np.sum(Itot0*bm)
    data[t,1+s] = np.sum(Itot0*em)
    data[t,2+s:4+s] = barycenter(Itot0,mask=bm)
    data[t,4+s:6+s] = barycenter(Itot0,mask=em)

    # print((("%i %f %f") %  (t,dx,dy)))
    
    

# saving data
np.save(fname+'.npy',data)

figure(0)
clf()
plot(data[:,0],data[:,33],'k')
plot(data[:,0],data[:,3],'r')

show()

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
