__author__ = 'fredpiech'

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo

#remember generic print function - func.printgraph (image, xrange, yrange, xlabel, ylabel)



#Define size of view and resolution
z = 20
theta = 4.

size = 200 # box size - maybe write a code to get this out of the title of the 21cmfast files
ci = 100

dtheta = float(theta/size)
print dtheta

eps = 1    #this is instrument efficiency

#DEPENDS ON Z! FIND THIS
B = 8000000        #Bandwidth (Hz) - taking this estimate for now, from MG Santos








zhat= func.twentyonecmmatrix('delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc',theta/2)







#Define Wavelength
lda=15

# Now we import array positions
(dx,dy,dz)=func.importarray('vla.a.cfg',lda) #'MWAcoordinate.txt''vla.a.cfg'

#plt.scatter(dx, dy)
#plt.show()

#DEFINE EXPERIMENT PARAMETERS
H = 0.

tint = 60.      #interval in seconds
dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
totalintegrationtime = 120    #total time in hours
timestepsneeded= int(totalintegrationtime * 60 * 24 / tint) # unitlessmeasurement of number of steps needed
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
