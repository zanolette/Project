__author__ = 'fredpiech'

######################### From Zanolette #########################
import numpy as np
import matplotlib.pyplot as plt

dtheta = 15./(60.**2)   #these will be given by 21cm later on
theta = 8. #gives 960 points
size = int(theta/dtheta)

# specify circle parameters: centre ij and radius
z = np.zeros((size,size))
ci,cj=int(size/2),(size/2) #centre point
cr=100  #arbitary for now

# Create index arrays to z
I,J=np.meshgrid(np.arange(z.shape[0]),np.arange(z.shape[1]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
z[np.where(np.sqrt((I-ci)**2+(J-cj)**2)<cr)]=100


#Fast Fourier Transform of z
zhat=np.fft.fft2(z)
zhat=np.fft.fftshift(zhat)
y = np.fft.ifft2(zhat)    #y is check of re-inversed
y = np.real(y)
zhatabs=abs(zhat)
#zhat = np.log(zhat)
#print zhat

#######################################################################

# Now we import array positions

positions = []
with open('vla.a.cfg') as configa:
    for line in configa:
        positions.append(line.strip().split('\t'))

dx=[]
dy=[]
dz=[]

lenpos=len(positions)


lda=10


#Set up dx dy dz vector - this is the baseline vector for
for i in range(lenpos):
    for j in range(lenpos):
        if i != j:
            dx.append(((float(positions[i][0])-float(positions[j][0]))/(lda)))  #doesn't take into account rotation of earth or angle of object
            dy.append(((float(positions[i][1])-float(positions[j][1]))/(lda)))
            dz.append(((float(positions[i][2])-float(positions[j][2]))/(lda)))  #do we need the lda or scaling!!!!!

#float(positions[i][0]) = x coord of ith element
#float(positions[i][1]) = y coord of ith element
plt.scatter(dx, dy)
plt.show()

#this is putting this into uv
H = 0.
dH = 1. / (60.* 24.) * np.pi
integrationtime = 60
delta = 30./180. * np.pi
lenbase = len(dx)
uvplane = np.zeros((3,integrationtime*lenbase ))


scaling = 1./(0.5*size*dtheta)

#This measures the UV plane AND maps it onto our fourier plane
for t in range (integrationtime):
    for i in range(lenbase):
        uvplane[0][i + t*lenbase] = int(ci + scaling*(np.sin(H)*dx[i] + np.cos(H)*dy[i])) # have UV matrix with enough space to get 24 hours of integration
        #t*lenbase makes sure we dont overwrite previous hours of integration
        uvplane[1][i + t*lenbase] = int(cj + scaling*(-np.sin(delta)*np.cos(H)*dx[i] + np.sin(delta)*np.sin(H)*dy[i] + np.cos(delta)*dz[i]))
        uvplane[2][i + t*lenbase] = np.cos(delta)*np.cos(H)*dx[i] - np.cos(delta)*np.sin(H)*dy[i] + np.sin(delta)*dz[i]
    H += dH

plt.scatter(uvplane[:][0], uvplane[:][1])
plt.show()

####################################MEASURE##################################################
#create the measured array
image = np.zeros((2*ci,2*cj),'complex')

countvar = 0.
#this uses baseline positions to sample - dx/dy[i] gives baseline positions, which then sample the fourier space zhat
for i in range(len(uvplane[0])):
    if uvplane[0][i] < size and uvplane[1][i] < size and uvplane[0][i] > 0 and uvplane[1][i] > 0:   #this is to stop crashing
        image[uvplane[0][i], uvplane[1][i]] = zhat[uvplane[0][i], uvplane[1][i]]
    else:
        #print ("error, image out of zhat bounds, baseline is", uvplane[0][i] , uvplane[1][i])
        countvar += 1.

print ("percentage of baselines ignored", (100*countvar/(len(uvplane[0]))))

#shows sample fourier image
image2 = np.log(abs(image)+1)
image2 = plt.imshow(image2, cmap='hot')


plt.colorbar(image2, orientation='horizontal')
plt.show()

#shows inverse image
imageinv = np.fft.ifft2(image)
imageinv2 = abs(imageinv)
imageinv2 = plt.imshow(imageinv2, cmap='hot')
plt.colorbar(imageinv2, orientation='horizontal')
plt.show()