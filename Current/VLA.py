__author__ = 'fredpiech'

######################### From Zanolette #########################
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
(dx,dy,dz)=func.importarray('vla.a.cfg',lda)

#this is putting this into uv
H = 0.
dH = 1. / (60.* 24.) * np.pi
integrationtime = 6*60
delta = 30./180. * np.pi
lenbase = len(dx)
uvplane = np.zeros((3,integrationtime*lenbase ))


scaling = 1./(0.5*size*dtheta)

#This measures the UV plane AND maps it onto our fourier plane
for t in range (integrationtime):
    for i in range(lenbase):
        uvplane[0][i + t*lenbase] = int(ci + scaling*(np.sin(H)*dx[i] + np.cos(H)*dy[i])) # have UV matrix with enough space to get 24 hours of integration
        #t*lenbase makes sure we dont overwrite previous hours of integration
        uvplane[1][i + t*lenbase] = int(ci + scaling*(-np.sin(delta)*np.cos(H)*dx[i] + np.sin(delta)*np.sin(H)*dy[i] + np.cos(delta)*dz[i]))
        uvplane[2][i + t*lenbase] = np.cos(delta)*np.cos(H)*dx[i] - np.cos(delta)*np.sin(H)*dy[i] + np.sin(delta)*dz[i]
    H += dH

plt.scatter(uvplane[:][0], uvplane[:][1])
plt.show()

####################################MEASURE##################################################
#create the measured array
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

