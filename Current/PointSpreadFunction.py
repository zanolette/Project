#Code To print off the point spread function of an array

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo

CosmoUnits=Cosmo.CosmoUnits()

z = 20
theta = CosmoUnits.thetaboxsize(z,400)
size = 200

ci = int(size/2)

dtheta = float(theta/size)
print dtheta

arrayname = 'vla.b.cfg'
Array = 'VLA 24 Hours Exposure'

#Define Wavelength - find this out from z!!
lda=0.21106*(1+z)   #in meters

if (arrayname=='LOFAR66CoreHB.txt'):
    Aeff = min(lda**2/3.,1.5625)
    N = 16*24*66 #number of tiles in a station * number of stations in an array
else:   #hence assuming its MWA
    Aeff = min(lda**2/3.,14.5/16)
    N = 16*128

Atotal = Aeff*N  #instrument efficiency: here Aeff is effective area per dipole, N is no. dipoles per station

#DEPENDS ON Z! FIND THIS
dl = 400/size   #so this is how many Mpc per index of box
psdwidth = 1    #can change this!

#DEFINE EXPERIMENT PARAMETERS
H = 0.
tint = 300.      #interval in seconds
dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
totalintegrationtime = 24    #total time in hours

daysofscanning = 1# keep as 1 unless you want more than one day.
timestepsneeded= int(totalintegrationtime * 60 * 60 / tint) # unitlessmeasurement of number of steps needed

delta = 90./180. * np.pi    #declination angle
scaling = 1./(size*dtheta) #the length of one interval in inverse space
print scaling
#CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
#spatialfreq=np.fft.fftfreq(size, dtheta)
#print spatialfreq[1] - spatialfreq[0]
#print scaling


# Now we import array positions
(dx,dy,dz)=func.importarray(arrayname,lda) #this assumes lda doesn't change too much over slices!!! 'MWAcoordinate.txt''vla.a.cfg' 'MWAhalved.txt''MWA128.txt'

#this calculated the longest baseline
baselinelengths=np.sqrt(np.multiply(dx,dx)+np.multiply(dy,dy))
Dmax = np.max(baselinelengths)
print Dmax
baselinelengths = None

eps = Dmax**2/Atotal
print eps

#Apply rotation matrix onto baseline vector and maps onto fourier plane.
UVcount,uvcoveragepercentage = func.rotationmatrix(dx, dy, dz, scaling, H, dH, timestepsneeded, delta, size)
UVcount = UVcount.astype(np.float32, copy=False)

image2 = plt.imshow(UVcount, extent=(-100*(lda/1000)/scaling,100*(lda/1000)/scaling,-100*(lda/1000)/scaling,100*(lda/1000)/scaling), interpolation='nearest',cmap='binary')
#plt.colorbar( orientation='vertical')
plt.xlabel('Kilo Wavelengths')
plt.ylabel('Kilo Wavelengths')
plt.title('%s UV Coverage'%Array)
plt.savefig('UV%s.png' %Array)
plt.clf()

psf = np.zeros((size,size))

#takes non-zero values to 1 as concerned with distribution in uv not count after rotations
for i in range (size):
    for j in range (size):
        if UVcount[i][j] != 0:
            psf[i][j]=1

image2 = plt.imshow(psf, extent=(-100*(lda/1000)/scaling,100*(lda/1000)/scaling,-100*(lda/1000)/scaling,100*(lda/1000)/scaling), interpolation='nearest',cmap='binary')
#plt.colorbar( orientation='vertical')
plt.xlabel('Kilo Wavelengths')
plt.ylabel('Kilo Wavelengths')
plt.title('%s UV Coverage'%Array)
plt.savefig('UVBinary%s.png' %Array)
plt.clf()

#UVcount = UVcount/UVcount

func.psf(dl, psf,size, 'PSF%s.png'%Array)
func.psfcrosssection(dl, psf,size,Array )
