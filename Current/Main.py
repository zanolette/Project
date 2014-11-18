

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
fname = 'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'
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
dl = box_info['BoxSize']/size
psdwidth = 3    #can change this!


#DEFINE EXPERIMENT PARAMETERS
H = 0.
tint = 60.      #interval in seconds
dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
totalintegrationtime = 1    #total time in hours
timestepsneeded= 1#int(totalintegrationtime * 60 * 24 / tint) # unitlessmeasurement of number of steps needed
delta = 90./180. * np.pi    #declination angle
scaling = 1./(size*dtheta)


#Define Wavelength - find this out from z!!
lda=.21106*(1+z)


# Now we import array positions
(dx,dy,dz)=func.importarray('MWAhalved.txt',lda) #this assumes lda doesn't change too much over slices!!! 'MWAcoordinate.txt''vla.a.cfg' 'MWAhalved.txt'


#Apply rotation matrix onto baseline vector and maps onto fourier plane.
UVcount = func.rotationmatrix(dx, dy, dz, scaling, H, dH, timestepsneeded, delta, size)


# need to make sure this is correct according to pritchard - noise on complex thing
tsyst = 50000 + 60000*((1+z)/4.73)**2.55  #(mK) this is from "Probing . . . with the SKA" MG Santos
B = 1420.41e6*dl/((1+z)*CosmoUnits.Dcomovingrad(z))        #Bandwidth (Hz) - this is given by frequency inteval. DeltaF = f1-f2, seperated by deltaz = z* deltaL/Dcomovingrad

sigma3D = np.zeros((size,size,size))    #initialised so we can save all slices here
image3D = np.zeros((size,size,size))

#imports 21cm box and takes single z slice
#zhat,twenty1 = func.twentyonecmmatrix(fname,theta/2)    #z will be used later to compute rms of image and measured sky
box=boximport.readbox(fname)
twenty1 = box.box_data  #so we have it seperate - this is 3D!

##########################This is where we introduce a slice################################

for slice in range(size):   #iterates over all slices

    print slice

    zhat=np.fft.fft2(twenty1[slice])     #then FFT of the specific slice
    zhat=np.fft.fftshift(zhat)  #this gives us the FFT of twenty1

    '''
    #gets 2D psd of image then gets 1D radial psd of this slice - this is INPUT power spectrum
    abszhat=np.abs(zhat)**2
    radialpsd = func.powerspectrum2D(abszhat,psdwidth)
    del abszhat
    spatialfreq=np.fft.fftfreq(int(size/psdwidth), dtheta)
    spatialfreq=spatialfreq[:int(size/(psdwidth*2))]    #this is used to give axis for power spectrum plots


    #this is power spectrum of 21cmbox
    fig=plt.loglog(spatialfreq,radialpsd)
    plt.xlabel('k')
    plt.ylabel('P(k)')
    plt.show()
    '''


    #CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
    #spatialfreq=np.fft.fftfreq(size, dtheta)
    #print spatialfreq[4] - spatialfreq[3]
    #print scaling


    ####################################MEASURE##################################################


    # for loop that goes through the fourier space matrix and adds noise according to the number of times the signal gets sampled


    #using UV count - this now merges the UVcoverage and the Image

    image = np.zeros((size,size),'complex')
    sigma = np.zeros((size,size))

    for i in range (size):
        for j in range (size):
            if UVcount[i][j] != 0:
                sigma[i][j] = tsyst/(eps*np.sqrt(UVcount[i][j]*tint*B))       #saved seperately to calculate power spectrum seperately, error eqn according to NRAO course + Pritchard
                real=np.random.normal(np.real(zhat[i][j]), sigma[i][j]/np.sqrt(2), 1)     #sqrt(2) here as real and imag components share it
                imaginary = np.random.normal(np.imag(zhat[i][j]), sigma[i][j]/np.sqrt(2), 1)
                image[i][j]=real[0] + imaginary[0]*1j

    '''
    #gets 2D psd of image then gets 1D radial psd - this is OUTPUT power spectrum
    absimage=np.abs(image)**2
    radialpsd2 = func.powerspectrum2D(absimage,psdwidth)
    noiseradialpsd2 = func.powerspectrum2D(sigma,psdwidth)
    fig1=plt.loglog(spatialfreq,noiseradialpsd2)
    plt.loglog(spatialfreq,radialpsd2)
    plt.loglog(spatialfreq,radialpsd)
    plt.xlabel('k')
    plt.ylabel('P(k)')
    plt.show()

    #THIS IS TO FIND THE PSF
    func.psfcrosssection(dtheta, image)
    '''

    sigma = func.invert(sigma, dtheta)  #at the moment this overwrites inverse space image after its inverted
    sigma3D[slice] = sigma  #saves this image slice in real space

    image = func.invert(image, dtheta)  #at the moment this overwrites inverse space image after its inverted
    image3D[slice] = image  #saves this image slice in real space

np.save('image3Darraydim%s,%sMpc,z%s'%(size,box_info['BoxSize'],z),image3D)
np.save('sigma3Darraydim%s,%sMpc,z%s'%(size,box_info['BoxSize'],z),sigma3D)

##############################IS DONE IN ANALYSE BOXES NOW##########################################

'''
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
'''