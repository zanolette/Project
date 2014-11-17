
import numpy as np
import pylab as pl
from matplotlib import ticker
import matplotlib.pyplot as plt
import boxio as boximport

def printgraph (image, xrange, yrange, xlabel, ylabel, scalemin,scalemax):   #generic print function with ranges and labels
    if scalemin == 'None':
        if scalemax == 'None':  #seems convoluted but allows us to use scalemin or max or not
            image2 = plt.imshow(image, extent=(-xrange,xrange,-yrange,yrange), interpolation='nearest',cmap='jet')
        else:
            image2 = plt.imshow(image, extent=(-xrange,xrange,-yrange,yrange),vmax=scalemax, interpolation='nearest',cmap='jet')
    else:
        if scalemax == 'None':
            image2 = plt.imshow(image, extent=(-xrange,xrange,-yrange,yrange),vmin=scalemin,interpolation='nearest',cmap='jet')
        else:
            image2 = plt.imshow(image, extent=(-xrange,xrange,-yrange,yrange),vmin=scalemin,vmax=scalemax, interpolation='nearest',cmap='jet')

    #70 is chosen arbitarily from z=14 21cmbox which goes to 64, may need to change

    plt.colorbar( orientation='vertical')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def twentyonecmmatrix(filename,theta):  #imports 21cm box and takes single z slice then FFT's
    box=boximport.readbox(filename)

    ci=int(box.dim/2)

    twenty1=box.box_data[199]

    printgraph(twenty1, theta,theta,"theta","theta", 0, 70)

    zhat=np.fft.fft2(twenty1)     #then FFT
    zhat=np.fft.fftshift(zhat)
    return zhat,twenty1     #as twenty1 will be used in rms calculation


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


def radial_data(data,width): #here width means how many pixels wide each band is
    size = len(data[0])
    rmax = int(size/2.)  #also same as centre

    #starts at rmax, adds all pixel counts for sqrt(x^2 + y^2) < r and number of pixels counted

    countarray = np.zeros((2,int(rmax/width)))    #this is due to slits of width 'width' and radius is rmax

    #this goes over all points in image, and adds the cumulative counts and number of counts for each band of k, starting with < r=1
    for r in range (rmax):
        for i in range (-rmax,rmax,1):
            for j in range (-rmax,rmax,1):
                if np.sqrt(i**2+j**2) < r+1 and np.sqrt(i**2+j**2) > r:    #r+1 as starts at 0 and goes to rmax-1
                    countarray[1][int(r/width)-1] += 1  #not the -1 is a fudge to get it to stay in bounds, but only determining where data is put
                    countarray[0][int(r/width)-1] += data[rmax + i][rmax + j]

    return countarray[0]/countarray[1]


#Maps baselines onto UV plane along with rotational matrix

#FRED CHANGED THIS!!!  - Now this functions merely takes all of the locations and maps out the UV plane completely together with a count
#We will have to mask it onto the image later.
#Pritchards suggestion to make this faster...

def rotationmatrix(dx, dy, dz, scaling, H, dH, integrationtime, delta, size):

    ci = int(size/2)
    lenbase = len(dx)

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
                UVcount[uvplane[0][i], uvplane[1][i]] += 1
            else:
                #print ("error, image out of zhat bounds, baseline is", uvplane[0][i] , uvplane[1][i])
                countvar += 1.
        #plt.scatter(uvplane[0],uvplane[1])
        #plt.show()
        H += dH

    print ("UV Plane Scan Complete, percentage of baselines ignored", (100*countvar/(integrationtime*lenbase)))
    return (UVcount)


def psf(dtheta,image):
    size = len(image[0])

    psf = np.zeros((size,size))

    for i in range (size):
        for j in range (size):
            if image[i][j] != 0:
                psf[i][j]=1

    imageinv = np.fft.ifft2(psf)
    imageinv = np.fft.fftshift(imageinv)
    imageinv = abs(imageinv)

    RangeinRealImage = size*dtheta/2.

    #fig = plt.figure(figsize=(6, 3.2))

    printgraph(imageinv, RangeinRealImage,RangeinRealImage,"theta","theta", 'None', 'None')

    #plt.imshow(imageinv,extent=(-RangeinRealImage,RangeinRealImage,-RangeinRealImage,RangeinRealImage),  interpolation='nearest', cmap='gist_stern')
    #plt.colorbar( orientation='horizontal')
    #plt.show()


def psfcrosssection(dtheta, image):

    size = len(image[0])
    ci=int(size/2)

    psf = np.zeros((size,size))

    for i in range (size):
        for j in range (size):
            if image[i][j] != 0:
                psf[i][j]=1

    psf = np.fft.ifft2(psf)
    psf = np.fft.fftshift(psf)
    psf = abs(psf)

    RangeinRealImage = size*dtheta/2.

    xaxis=np.arange(-RangeinRealImage, RangeinRealImage, dtheta)
    #xaxis=((xaxis/size)-0.5)*RangeinRealImage*2
    #xaxis=xaxis*RangeinRealImage/(size/2)

    crosssection =psf[ci]
    psfimage=plt.plot(xaxis,crosssection)
    plt.ylabel('PSF Amplitude')
    plt.xlabel('theta')
    plt.show()


#This function takes matrix relating to the fourier space and inverts it. It plots both the fourier image and the real space image
def invert(image, dtheta):

    RangeinComplexImage = 0.5/dtheta

    image2 = np.log(abs(image)+1)   #!!!LOGGED!!! to make clear

    printgraph(image2, RangeinComplexImage,RangeinComplexImage,"1/theta","1/theta",'None','None')

    #image2 = plt.imshow(image2, extent=(-RangeinComplexImage,RangeinComplexImage,-RangeinComplexImage,RangeinComplexImage), interpolation='nearest',cmap='hot')
    #plt.colorbar( orientation='horizontal')
    #plt.xlabel("1./theta")
    #plt.ylabel("1./theta")
    #plt.show()


    #shows inverse image
    imageinv = np.fft.ifft2(image)
    imageinv = abs(imageinv)

    RangeinRealImage = (len(image[1])*dtheta)/2


    printgraph(imageinv, RangeinRealImage,RangeinRealImage,"theta","theta", 0,70)


    #plt.imshow(imageinv,extent=(-RangeinRealImage,RangeinRealImage,-RangeinRealImage,RangeinRealImage),  interpolation='nearest', cmap='hot')
    #plt.colorbar( orientation='horizontal')
    #plt.xlabel("theta")
    #plt.ylabel("theta")
    #plt.show()

def rmscalc (twenty1cm,image):
    max = len(image[0])
    squaredcount = 0    #this counts x**2 + y**2 . . . then will average and sqrt at end

    for i in range(max):
        for j in range(max):
            squaredcount += (np.real(image[i][j]) - twenty1cm[i][j])**2     #or should this be real(image)

    squaredcount = squaredcount/(max**2)
    return np.sqrt(squaredcount)
