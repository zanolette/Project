

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
fname = 'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'#'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'#
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
timestepsneeded= 1 #int(totalintegrationtime * 60 * 24 / tint) # unitlessmeasurement of number of steps needed
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

#imports 21cm box and takes single z slice
#zhat,twenty1 = func.twentyonecmmatrix(fname,theta/2)    #z will be used later to compute rms of image and measured sky
box=boximport.readbox(fname)
twenty1 = box.box_data  #so we have it seperate - this is 3D!


twenty1inverse = np.fft.fftn(twenty1)   #gives 3D FFT of 21cm box!
twenty1inverse = np.fft.fftshift(twenty1inverse)

image3Dinverse = np.zeros((size,size,size),'complex')   #so we can save into a 3D box in fourier space
sigma3Dinverse = np.zeros((size,size,size))

##########################This is where we introduce a slice################################

for slice in range(size):   #iterates over all slices

    print slice

    zhat = twenty1inverse[slice]    #takes a 2D slice of 21cm in fourier space


    #CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
    #spatialfreq=np.fft.fftfreq(size, dtheta)
    #print spatialfreq[4] - spatialfreq[3]
    #print scaling


    ####################################MEASURE##################################################


    # for loop that goes through the fourier space matrix and adds noise according to the number of times the signal gets sampled


    #using UV count - this now merges the UVcoverage and the Image

    for i in range (size):
        for j in range (size):
            if UVcount[i][j] != 0:
                sigma3Dinverse[slice][i][j] = tsyst/(eps*np.sqrt(UVcount[i][j]*tint*B))       #saved seperately to calculate power spectrum seperately, error eqn according to NRAO course + Pritchard
                real=np.random.normal(np.real(zhat[i][j]), sigma3Dinverse[slice][i][j]/np.sqrt(2), 1)     #sqrt(2) here as real and imag components share it
                imaginary = np.random.normal(np.imag(zhat[i][j]), sigma3Dinverse[slice][i][j]/np.sqrt(2), 1)
                image3Dinverse[slice][i][j]=real[0] + imaginary[0]*1j


#func.phasecomparison(twenty1inverse, image3Dinverse, size)
'''
image3D = np.fft.ifftn(image3Dinverse)
image3D = np.abs(image3D)

#func.visualizereionizationagainstz(image3D, size, z, theta)

distimagex, distimagey= func.bubblesizedistribution(image3D, size,0.5)
dist21x, dist21y= func.bubblesizedistribution(twenty1, size,0.5)

ionizedsites = np.sum(dist21x*dist21y)
ionizedproportion=ionizedsites/(size**3)

print 'Ionized Proportion is:'
print ionizedproportion

plt.loglog(distimagex,distimagey)
plt.loglog(dist21x, dist21y)
plt.xlabel('Bubble Size')
plt.ylabel('Number of Bubbles')
plt.show()

'''
#THIS IS TO FIND THE PSF
#func.psfcrosssection(dtheta, image3Dinverse[int(size/2.)],size)
'''
(kaxis,Powerspectrum) = func.powerspectrum3D(twenty1inverse,psdwidth,size,dtheta)
plt.loglog(kaxis,Powerspectrum)#*(kaxis**3)/(2*np.pi**2))
print 1
(kaxis,Powerspectrum)=func.powerspectrum3D(image3Dinverse,psdwidth,size,dtheta)
plt.loglog(kaxis,Powerspectrum)#*(kaxis**3)/(2*np.pi**2))
print 2
(kaxis,Powerspectrum)=func.powerspectrum3D(sigma3Dinverse,psdwidth,size,dtheta)
plt.loglog(kaxis,Powerspectrum)#*(kaxis**3)/(2*np.pi**2))
print 3
plt.xlabel('k')
plt.ylabel('P(k)')

plt.show()
'''
    #THIS IS TO FIND THE PSF
    #func.psfcrosssection(dtheta, image,size)




#image3D = np.fft.ifftn(image3Dinverse)
#dont think we need to shift -
#image3D = np.abs(image3D)

#sigma3D = np.fft.ifftn(sigma3Dinverse)
#######do we need a shift here?###########
#sigma3D = np.abs(sigma3D)

#func.visualizereionizationslicebyslice(image3D,twenty1, size, z, theta)



#np.save('image3Darraydim%s,%sMpc,z%s,test'%(size,box_info['BoxSize'],z),image3D)
#np.save('sigma3Darraydim%s,%sMpc,z%s,test'%(size,box_info['BoxSize'],z),sigma3D)

#np.save('image3Darraydim%s,%sMpc,z%s,time%s'%(size,box_info['BoxSize'],z,timestepsneeded),image3D)
#np.save('sigma3Darraydim%s,%sMpc,z%s,time&s'%(size,box_info['BoxSize'],z,timestepsneeded),sigma3D)


##############################POWER SPECTRUM##########################################
#download real power spectrum
realps = np.loadtxt('ps_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_100_200Mpc_v2.txt', delimiter='\t')

imagek, imagepowerspectrum , imagedeldel= func.powerspectrum3D(image3Dinverse,psdwidth,size,dtheta,float(box_info['dim'])/float(box_info['BoxSize']), z) # this is the size of steps in real space dx=float(box_info['dim'])/float(box_info['BoxSize'])
print 'done imagepowerspectrum'
print len(imagek), len(imagepowerspectrum)
sigmak, sigmapowerspectrum , sigmadeldel= func.powerspectrum3D(sigma3Dinverse,psdwidth,size,dtheta, float(box_info['dim'])/float(box_info['BoxSize']), z)
print 'done sigmapowerspectrum'
twenty1k, twenty1powerspectrum, twenty1deldel= func.powerspectrum3D(twenty1inverse,psdwidth,size,dtheta,float(box_info['dim'])/float(box_info['BoxSize']), z)
print 'done twenty1powerspectrum'



#spatialfreq=np.fft.fftfreq(int(size/psdwidth), dtheta)
#spatialfreq=spatialfreq[:int(size/(psdwidth*2))]    #this is used to give axis for power spectrum plots
plt.loglog(imagek,imagedeldel)
plt.loglog(sigmak,sigmadeldel)
plt.loglog(twenty1k,twenty1deldel)
plt.loglog(realps[:,0],realps[:,1])
plt.xlim(0.02,3)

plt.xlabel('k (MPc$^{-1}$)')
plt.ylabel('k$^3$ P(k)/2$\pi^2$')
plt.show()



plt.loglog(imagek,imagepowerspectrum)
plt.loglog(sigmak,sigmapowerspectrum)
plt.loglog(twenty1k,twenty1powerspectrum)
plt.loglog(realps[:,0],realps[:,1]/(realps[:,0]**3))    #Important: here, as with other plot, we have no 2Pi**2 factor
plt.ylim(0,10000)
plt.xlim(0.02,3)

plt.xlabel('k (MPc$^{-1}$)')
plt.ylabel('P(k)')
plt.show()
