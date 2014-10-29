__author__ = 'fredpiech'

import numpy as np
import matplotlib.pyplot as plt
import Functions as func

#Define size of view and resolution
dtheta = 15./(60.**2)   #these will be given by 21cm later on
theta = 8. #gives 960 points
size = int(theta/dtheta)
ci=int(size/2)

#creates fourier space matrix of circle radius 100
zhat = func.circlematrix(size, 100)

#Define Wavelength
lda=10

# Now we import array positions
(dx,dy,dz)=func.importarray('vla.b.cfg',lda)

plt.scatter(dx, dy)
plt.show()

#DEF
H = 0.
dH = 1. / (60.* 24.) * np.pi
integrationtime = 6*60
delta = 30./180. * np.pi
scaling = 1./(0.5*size*dtheta)


uvplane=func.rotationmatrix(dx, dy, dz, scaling, H, dH, integrationtime, delta, ci)

plt.scatter(uvplane[:][0], uvplane[:][1])
plt.show()

####################################MEASURE##################################################

image = np.zeros((2*ci,2*ci),'complex')
PSFimage = np.zeros((2*ci,2*ci),'complex')

countvar = 0.
#this uses baseline positions to sample - dx/dy[i] gives baseline positions, which then sample the fourier space zhat
for i in range(len(uvplane[0])):
    if uvplane[0][i] < size and uvplane[1][i] < size and uvplane[0][i] > 0 and uvplane[1][i] > 0:   #this is to stop crashing
        image[uvplane[0][i], uvplane[1][i]] += zhat[uvplane[0][i], uvplane[1][i]]
        #PSFimage counts how often a point is traced in the uv plane
        PSFimage[uvplane[0][i], uvplane[1][i]] += 1
    else:
        #print ("error, image out of zhat bounds, baseline is", uvplane[0][i] , uvplane[1][i])
        countvar += 1.

print ("percentage of baselines ignored", (100*countvar/(len(uvplane[0]))))


#shows sample fourier image
func.invert(image)
func.invert(PSFimage)