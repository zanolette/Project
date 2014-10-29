__author__ = 'fredpiech'
import numpy as np
import matplotlib.pyplot as plt

#creates fourier space matrix of circle radius "radius"
def circlematrix(size, radius):
    # specify circle parameters: centre ij and radius
    z = np.zeros((size,size))

    ci,cj=int(size/2),(size/2) #centre point
    cr=radius  #arbitary for now

    # Create index arrays to z
    I,J=np.meshgrid(np.arange(z.shape[0]),np.arange(z.shape[1]))

    # np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
    # Assign value of 1 to those points where distance to centre is less that cr:
    z[np.where(np.sqrt((I-ci)**2+(J-cj)**2)<cr)]=100

    #Fast Fourier Transform of z
    zhat=np.fft.fft2(z)
    zhat=np.fft.fftshift(zhat)

    #NOT SURE WHAT THESE WERE FOR
    #y = np.fft.ifft2(zhat)    #y is check of re-inversed
    #y = np.real(y)
    #zhatabs=abs(zhat)

    return zhat


#Imports and array from a text file and returns the baseline dx, dy and dz. depends on lambda
def importarray(filename,lda):
    positions = []
    with open(filename) as configa:
        for line in configa:
            positions.append(line.strip().split('\t'))

    dx=[]
    dy=[]
    dz=[]

    lenpos=len(positions)

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

    return (dx,dy,dz)

#This function takes matrix relating to the fourier space and inverts it. It plots both the fourier image and the real space image
def invert(image):
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