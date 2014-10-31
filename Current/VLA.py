__author__ = 'fredpiech'

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func

#Define size of view and resolution
dtheta = 15./(60.**2)   #these will be given by 21cm later on
theta = 8. #gives 960 points
size = int(theta/dtheta)
ci=int(size/2)



#creates fourier space matrix of circle radius 100
zhat = func.circlematrix(size, 100)

#Define Wavelength
lda=1

# Now we import array positions
(dx,dy,dz)=func.importarray('vla.b.cfg',lda)

#plt.scatter(dx, dy)
#plt.show()

#DEFINE EXPERIMENT PARAMETERS
H = 0.
dH = 1. / (60.* 24.) * np.pi
integrationtime = 30
delta = 30./180. * np.pi
scaling = 1./(size*dtheta)


#CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
#spatialfreq=np.fft.fftfreq(size, dtheta)
#print spatialfreq[4] - spatialfreq[3]
#print scaling


#Apply rotation matrix onto baseline vector and maps onto fourier plane.
(image, UVcount) = func.rotationmatrix(dx, dy, dz, scaling, H, dH, integrationtime, delta, size, zhat)

####################################MEASURE##################################################

# for loop that goes through the fourier space matrix and adds noise according to the number of times the signal gets sampled
# need to make sure this is correct according to pritchard - noise on complex thing
for i in range (size):
    for j in range (size):
        if image[i][j] != 0:
            sigma = 50000/np.sqrt(UVcount[i][j])
            real=np.random.normal(np.real(image[i][j]), sigma, 1)
            imaginary = np.random.normal(np.imag(image[i][j]), sigma, 1)
            image[i][j]=real[0] + imaginary[0]*1j





#shows sample fourier image
func.invert(image, dtheta, theta, size)
#func.invert(UVcount)

#i'm testing this