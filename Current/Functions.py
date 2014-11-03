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

    return (dx,dy,dz)

#Maps baselines onto UV plane along with rotational matrix
def rotationmatrix(dx, dy, dz, scaling, H, dH, integrationtime, delta, size, zhat):
    ci = int(size/2)
    lenbase = len(dx)

    image = np.zeros((size,size),'complex')
    UVcount = np.zeros((size,size),'complex')

    countvar = 0.

    #This measures the UV plane AND maps it onto our fourier plane
    for t in range (integrationtime): #iterates through time

        uvplane = np.zeros((3,lenbase ))

        for i in range(lenbase): #one application of the rotation matrix
            uvplane[0][i] = int(ci + scaling*(np.sin(H)*dx[i] + np.cos(H)*dy[i])) # have UV matrix with enough space to get 24 hours of integration
            #t*lenbase makes sure we dont overwrite previous hours of integration
            uvplane[1][i] = int(ci + scaling*(-np.sin(delta)*np.cos(H)*dx[i] + np.sin(delta)*np.sin(H)*dy[i] + np.cos(delta)*dz[i]))
            uvplane[2][i] = np.cos(delta)*np.cos(H)*dx[i] - np.cos(delta)*np.sin(H)*dy[i] + np.sin(delta)*dz[i]

        #Sample fourier space positions (zhat) for this given time
        for i in range(len(uvplane[0])):
            if uvplane[0][i] < size and uvplane[1][i] < size and uvplane[0][i] > 0 and uvplane[1][i] > 0 :   #this is to stop crashing
                if image[uvplane[0][i], uvplane[1][i]]==0: # this to decrease run time by not going over same point twice
                    image[uvplane[0][i], uvplane[1][i]] = zhat[uvplane[0][i], uvplane[1][i]]
                #PSFimage counts how often a point is traced in the uv plane
                UVcount[uvplane[0][i], uvplane[1][i]] += 1
            else:
                #print ("error, image out of zhat bounds, baseline is", uvplane[0][i] , uvplane[1][i])
                countvar += 1.

        H += dH

    print ("UV Plane Scan Complete, percentage of baselines ignored", (100*countvar/(integrationtime*lenbase)))
    return (image, UVcount)

#This function takes matrix relating to the fourier space and inverts it. It plots both the fourier image and the real space image
def invert(image, dtheta):

    RangeinComplexImage = len(image[1])/dtheta

    image2 = np.log(abs(image)+1)
    image2 = plt.imshow(image2, extent=(-RangeinComplexImage,RangeinComplexImage,-RangeinComplexImage,RangeinComplexImage), interpolation='nearest',cmap='hot')

    plt.colorbar( orientation='horizontal')
    plt.show()

    #shows inverse image
    imageinv = np.fft.ifft2(image)
    imageinv = abs(imageinv)

    #ax = imageinv2.add_subplot(111)

    #imageinv2 = plt.imshow(imageinv2, cmap='hot')

    #plt.colorbar(imageinv2, orientation='horizontal')
    #plt.show()
    RangeinRealImage = (len(image[1])*dtheta)/2

    fig = plt.figure(figsize=(6, 3.2))

    ax = fig.add_subplot(111)
    plt.imshow(imageinv,extent=(-RangeinRealImage,RangeinRealImage,-RangeinRealImage,RangeinRealImage),  interpolation='nearest', cmap='hot')
    ax.set_aspect('equal')

    #cax = fig.add_axes([-RangeinRealImage,RangeinRealImage,-RangeinRealImage,RangeinRealImage])

    plt.colorbar(orientation='vertical')
    plt.show()