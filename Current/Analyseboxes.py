import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo

###############################Taken from Main#################################

psdwidth = 1

fname = 'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'#'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'
box_info = boximport.parse_filename(fname)

#Define size of view and resolution
# defines our universe
CosmoUnits=Cosmo.CosmoUnits()
z = box_info['z']
theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files
dtheta = float(theta/size)

box=boximport.readbox(fname)
twenty1 = box.box_data  #so we have it seperate - this is 3D!

###############################################################################

image3D = np.load('image3Darray%s.npy'% size)
sigma3D = np.load('sigma3Darray%s.npy'% size)


###########################POWERSPECTRUM########################################

#This calculates 3D powerspectrum, after all slices are done
imagepowerspectrum = func.powerspectrum3D(image3D,psdwidth,size)
print 'done imagepowerspectrum'
sigmapowerspectrum = func.powerspectrum3D(sigma3D,psdwidth,size)
print 'done sigmapowerspectrum'
twenty1powerspectrum = func.powerspectrum3D(twenty1,psdwidth,size)
print 'done twenty1powerspectrum'

spatialfreq=np.fft.fftfreq(int(size/psdwidth), dtheta)
spatialfreq=spatialfreq[:int(size/(psdwidth*2))]    #this is used to give axis for power spectrum plots
plt.loglog(spatialfreq,imagepowerspectrum)
plt.loglog(spatialfreq,sigmapowerspectrum)
plt.loglog(spatialfreq,twenty1powerspectrum)
plt.xlabel('k')
plt.ylabel('P(k)')
plt.show()

plt.loglog(spatialfreq,spatialfreq**3*imagepowerspectrum/(2*np.pi**2))
plt.loglog(spatialfreq,spatialfreq**3*sigmapowerspectrum/(2*np.pi**2))
plt.loglog(spatialfreq,spatialfreq**3*twenty1powerspectrum/(2*np.pi**2))
plt.xlabel('k')
plt.ylabel('k$^3$ P(k)/2$\pi^2$')
plt.show()

#Compute rms between image and inputed 21cm - only do this if willing to wait
print 'rms between 21cmbox and image is', func.rmscalc(twenty1,image3D,size)

################################################################################