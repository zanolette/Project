__author__ = 'fredpiech'

import numpy as np
import matplotlib.pyplot as plt

# specify circle parameters: centre ij and radius
ci,cj=500,500 #centre point
z = np.zeros((2*ci,2*cj))
cr=50

# Create index arrays to z
I,J=np.meshgrid(np.arange(z.shape[0]),np.arange(z.shape[1]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
z[np.where(np.sqrt((I-ci)**2+(J-cj)**2)<cr)]=1

ri,rj=800,800 #centre point
rr=10
z[np.where(np.sqrt((I-ri)**2+(J-rj)**2)<rr)]=2

ri,rj=200,400 #centre point
rr=30
z[np.where(np.sqrt((I-ri)**2+(J-rj)**2)<rr)]=1

ri,rj=300,300 #centre point
rr=50
z[np.where(np.sqrt((I-ri)**2+(J-rj)**2)<rr)]=1

ri,rj=800,150 #centre point
rr=10
z[np.where(np.sqrt((I-ri)**2+(J-rj)**2)<rr)]=2


#Fast Fourier Transform of z
zhat=np.fft.fft2(z)
zhat=np.fft.fftshift(zhat)
y = np.fft.ifft2(zhat)    #y is check of re-inversed
y = np.real(y)
zhatabs=abs(zhat)
#zhat = np.log(zhat)
#print zhat

'''
im = plt.imshow(zhat, cmap='hot')
plt.colorbar(im, orientation='horizontal')
plt.show()

new = plt.imshow(y, cmap='hot')
plt.colorbar(new, orientation='horizontal')
plt.show()
'''

lda = 0.18
#positions in km - ARBITARY LAMBDA = 0.18m

xpos = (np.random.rand(1,30)-0.5)*10
ypos = (np.random.rand(1,30)-0.5)*10
#xpos = [10, 10, -10, -10, 10, 0, -10, 0, 0, 20, 20, -20, -20, 20, 0, -20, 0, 5, 5, -5, -5, 5, 0, -5, 0, 10, 35, -25, -35, 7, 12, -24, 8]
#ypos = [10, -10, 10, -10, 0, 10, 0, -10, 0, 20, -20, 20, -20, 0, 20, 0, -20, 5, -5, 5, -5, 0, 5, 0, -5, 35, -10, 35, -25, 0, -23, 13, -32]
dx = []
dy = []

#this takes the positions and turns them into baselines (no rotation included yet), shifted so they correspond with centre of image
for i in range(np.size(xpos)):
    for j in range(np.size(ypos)):
        if i != j:
            dx.append(ci + int((xpos[0,i]-xpos[0,j])/lda))  #doesn't take into account rotation of earth or angle of object
            dy.append(cj + int((ypos[0,i]-ypos[0,j])/lda))

print dx
print dy

#create the measured array
image = np.zeros((2*ci,2*cj),'complex')

#this uses baseline positions to sample - dx/dy[i] gives baseline positions, which then sample the fourier space zhat
for i in range(len(dx)):
    image[dx[i], dy[i]] = zhat[dx[i], dy[i]]

#shows sample fourier image
image2 = abs(image)
image2 = plt.imshow(image2, cmap='hot')
plt.colorbar(image2, orientation='horizontal')
plt.show()

#shows inverse image
imageinv = np.fft.ifft2(image)
imageinv2 = abs(imageinv)
imageinv2 = plt.imshow(imageinv2, cmap='hot')
plt.colorbar(imageinv2, orientation='horizontal')
plt.show()