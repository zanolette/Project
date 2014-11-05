__author__ = 'fredpiech'

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func

#Define size of view and resolution
dtheta = 15./(60.**2)   #these will be given by 21cm later on
theta = 4. #gives 960 points
size = int(theta/dtheta)
ci=int(size/2)
eps = 1.    #this is instrument efficiency
z = 20
B = 10000        #Bandwidth (Hz) - taking this estimate for now, from MG Santos

#remember generic print function - func.printgraph (image, xrange, yrange, xlabel, ylabel)

#creates fourier space matrix of circle radius 100
zhat = func.circlematrix(size, 100)

#Define Wavelength
lda=10

# Now we import array positions
(dx,dy,dz)=func.importarray('vla.b.cfg',lda)

#plt.scatter(dx, dy)
#plt.show()

#DEFINE EXPERIMENT PARAMETERS
H = 0.
tint = 60.      #interval in seconds
dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
totalintegrationtime = 2    #total time in hours
timestepsneeded= int(totalintegrationtime * 60 * 60 / tint) # unitlessmeasurement of number of steps needed
delta = 30./180. * np.pi    #declination angle
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
tsyst = 50 + 60*((1+z)/4.73)**2.55  #this is from "Probing . . . with the SKA" MG Santos

for i in range (size):
    for j in range (size):
        if image[i][j] != 0:
            sigma = tsyst/(eps*np.sqrt(UVcount[i][j]*tint*B))       #error eqn according to NRAO course + Pritchard
            real=np.random.normal(np.real(image[i][j]), sigma, 1)
            imaginary = np.random.normal(np.imag(image[i][j]), sigma, 1)
            image[i][j]=real[0] + imaginary[0]*1j



#THIS IS TO FIND THE PSF
func.psf(dtheta, image)

func.invert(image, dtheta)
#shows sample fourier image



#i'm testing this
