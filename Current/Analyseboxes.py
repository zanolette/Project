import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo

###############################Taken from Main#################################

psdwidth = 2

fname = 'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'#'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'#
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

dx=float(box_info['dim'])/float(box_info['BoxSize'])

###############################################################################

#image3D = np.load('image3Darray%s.npy'% size)
#sigma3D = np.load('sigma3Darray%s.npy'% size)

#image3D = np.load('image3Darraydim%s,%sMpc,z%s,time%s.npy'%(size,box_info['BoxSize'],z,120))
#sigma3D = np.load('sigma3Darraydim%s,%sMpc,z%s,time%s.npy'%(size,box_info['BoxSize'],z,120))

image3D = np.load('image3Darraydim%s,%sMpc,z%s.npy'%(size,box_info['BoxSize'],z))
sigma3D = np.load('sigma3Darraydim%s,%sMpc,z%s.npy'%(size,box_info['BoxSize'],z))

image3Dinv=np.fft.fftn(image3D)
image3Dinv=np.fft.fftshift(image3Dinv)

sigma3Dinv=np.fft.fftn(sigma3D)
sigma3Dinv=np.fft.fftshift(sigma3Dinv)

twenty1inv = np.fft.fftn(twenty1)   #gives 3D FFT of 21cm box!
twenty1inv = np.fft.fftshift(twenty1inv)


realps = np.loadtxt('ps_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_100_200Mpc_v2.txt', delimiter='\t')
print realps[:,0]

###########################2Drep###############################################

#func.visualizereionizationslicebyslice(image3D,twenty1, size, z, theta)
#func.visualizereionization(twenty1, size, z, box_info['BoxSize'])

#################################################################################

#Gives evolution of 2D power spectrum
#func.powerspectrumevolution(image3D,psdwidth,size,dtheta)

###########################POWERSPECTRUM########################################


#This calculates 3D powerspectrum, after all slices are done
imagek, imagepowerspectrum = func.powerspectrum3D(image3Dinv,psdwidth,size,dtheta,dx, z)
print 'done imagepowerspectrum'
print len(imagek), len(imagepowerspectrum)
sigmak, sigmapowerspectrum = func.powerspectrum3D(sigma3Dinv,psdwidth,size,dtheta, dx, z)
print 'done sigmapowerspectrum'
twenty1k, twenty1powerspectrum = func.powerspectrum3D(twenty1inv,psdwidth,size,dtheta,dx, z)
print 'done twenty1powerspectrum'



#spatialfreq=np.fft.fftfreq(int(size/psdwidth), dtheta)
#spatialfreq=spatialfreq[:int(size/(psdwidth*2))]    #this is used to give axis for power spectrum plots
plt.loglog(imagek,imagepowerspectrum)
#plt.loglog(sigmak,sigmapowerspectrum)
plt.loglog(twenty1k,twenty1powerspectrum)
plt.loglog(realps[:,0],realps[:,1]*(2*np.pi**2)/(realps[:,0]**3))
plt.xlabel('k')
plt.ylabel('P(k)')
plt.show()

plt.loglog(imagek,(imagek)**3*imagepowerspectrum/(2*np.pi**2))
#plt.loglog(sigmak,sigmak**3*sigmapowerspectrum/(2*np.pi**2))
plt.loglog(twenty1k,(twenty1k)**3*twenty1powerspectrum/(2*np.pi**2))
plt.loglog(realps[:,0],realps[:,1])

plt.xlabel('k')
plt.ylabel('k$^3$ P(k)/2$\pi^2$')
plt.show()

#Compute rms between image and inputed 21cm - only do this if willing to wait
#print 'rms between 21cmbox and image is', func.rmscalc(twenty1,image3D,size)

################################################################################



print twenty1powerspectrum[0], realps[0,1]
