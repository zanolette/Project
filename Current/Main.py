

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo

# defines our universe
CosmoUnits=Cosmo.CosmoUnits()

#remember generic print function - func.printgraph (image, xrange, yrange, xlabel, ylabel,scalemin,scalemax)

#getting 21cm box information
fname = 'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'
box_info = boximport.parse_filename(fname)


#Define size of view and resolution
z = box_info['z']
theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
print theta

size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files

ci = int(size/2)

dtheta = float(theta/size)
print dtheta

eps = 0.5    #this is instrument efficiency

#DEPENDS ON Z! FIND THIS
B = 8000000        #Bandwidth (Hz) - taking this estimate for now, from MG Santos



#imports 21cm box and takes single z slice
zhat,twenty1 = func.twentyonecmmatrix(fname,theta/2)    #z will be used later to compute rms of image and measured sky


#gets 2D psd of image then gets 1D radial psd - this is INPUT power spectrum
psdwidth = 3    #can change this!
abszhat=np.abs(zhat)**2
radialpsd = func.radial_data(abszhat,psdwidth)
del abszhat
spatialfreq=np.fft.fftfreq(int(size/psdwidth), dtheta)
spatialfreq=spatialfreq[:int(size/(psdwidth*2))]

fig=plt.loglog(spatialfreq,radialpsd)
plt.xlabel('k')
plt.ylabel('P(k)')
plt.show()


#Define Wavelength - find this out from z!!
lda=.21106*(1+z)


# Now we import array positions
(dx,dy,dz)=func.importarray('MWAhalved.txt',lda) #'MWAcoordinate.txt''vla.a.cfg'


#DEFINE EXPERIMENT PARAMETERS
H = 0.
tint = 60.      #interval in seconds
dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
totalintegrationtime = 10    #total time in hours
timestepsneeded= 1#int(totalintegrationtime * 60 * 24 / tint) # unitlessmeasurement of number of steps needed
delta = 90./180. * np.pi    #declination angle
scaling = 1./(size*dtheta)


#CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
#spatialfreq=np.fft.fftfreq(size, dtheta)
#print spatialfreq[4] - spatialfreq[3]
#print scaling


#Apply rotation matrix onto baseline vector and maps onto fourier plane.
(image, UVcount) = func.rotationmatrix(dx, dy, dz, scaling, H, dH, timestepsneeded, delta, size, zhat)


####################################MEASURE##################################################


# for loop that goes through the fourier space matrix and adds noise according to the number of times the signal gets sampled
# need to make sure this is correct according to pritchard - noise on complex thing
tsyst = 50000 + 60000*((1+z)/4.73)**2.55  #(mK) this is from "Probing . . . with the SKA" MG Santos


for i in range (size):
    for j in range (size):
        if image[i][j] != 0:
            sigma = tsyst/(eps*np.sqrt(UVcount[i][j]*tint*B))       #error eqn according to NRAO course + Pritchard
            real=np.random.normal(np.real(image[i][j]), sigma, 1)
            imaginary = np.random.normal(np.imag(image[i][j]), sigma, 1)
            image[i][j]=real[0] + imaginary[0]*1j


#gets 2D psd of image then gets 1D radial psd - this is OUTPUT power spectrum
absimage=np.abs(image)**2
radialpsd2 = func.radial_data(absimage,psdwidth)

fig=plt.loglog(spatialfreq,radialpsd2)
plt.xlabel('k')
plt.ylabel('P(k)')
plt.show()

#THIS IS TO FIND THE PSF
func.psfcrosssection(dtheta, image)


#shows sample fourier image
func.invert(image, dtheta)

#Compute rms between image and inputed 21cm
print 'rms between 21cmbox and image is', func.rmscalc(twenty1,image)


