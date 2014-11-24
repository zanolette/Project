
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

    twenty1=box.box_data[199]

    #printgraph(twenty1, theta,theta,"theta","theta", 0, 70)

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

def powerspectrumevolution(data,width,size,dtheta):    #this is to show 2D power spectrum evolving with z
    rmax = int(size/2.)
    rsize = 1+int(np.sqrt(2.*rmax**2)/width)      #this is largest r box (taking into account shell widths)


    psevolution = np.zeros((size,rsize))    #this will store each 2D powerspectrum with its z

    #xaxis needs to be k
    spatialfreq=np.fft.fftfreq(2*rsize, dtheta)
    spatialfreq=spatialfreq[:rsize]


    for i in range (size):
        temp = powerspectrum2D(data[i],width,size)
        psevolution[i] = temp

        #saves the 2Dpowerspectrum file
        plt.plot(spatialfreq,temp)
        plt.ylim((12,40))
        plt.xlim((0,46))    #will need to adjust these manually if major changes made to input
        #plt.yticks(range(10,21,1))
        plt.xlabel('k')
        plt.ylabel('P(k)')
        plt.savefig('test/%i' %i)
        plt.clf()

        print i

##2D version##
def powerspectrum2D(data,width,size): #here width means how many pixels wide each band is

    rmax = int(size/2.)  #also same as centre

    #starts at rmax, adds all pixel counts for sqrt(x^2 + y^2) < r and number of pixels counted

    countarray = np.zeros((2,1 + int(np.sqrt(2.*rmax**2)/width)))    #this is due to slits of width 'width' and radius is rmax, so longest r is sqrt(r**2 + r**2)

    #this goes over all points in image, and adds the cumulative counts and number of counts for each band of k, starting with < r=1
    for i in range(size):
        for j in range (size):
            r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2)/width)   #works out how far it is from the centre, r/width so can save larger wedges

            countarray[1][r] += 1   #adds 1 to count,
            countarray[0][r] += data[i][j]

    return countarray[0]/countarray[1]

##3D version##
def powerspectrum3D(image3Dinv,width,size,dtheta): #here width means how many pixels wide each band is

    rmax = int(size/2.)  #also same as centre

    '''
    plt.imshow(image3D[50],  interpolation='nearest',cmap='jet')
    plt.colorbar( orientation='vertical')
    plt.show()
    '''

    #these 3 take image into fourier space and make |p(k)|^2
    #image3D = np.fft.fftn(image3D)
    #image3D = np.fft.fftshift(image3D)

    image3Dinv = np.abs(image3Dinv)**2


    print 'fourier transform done'

    countarray = np.zeros((2,1+int(np.sqrt(3.*rmax**2)/width)))

    print len(countarray[0])

    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2 + (k-rmax)**2)/width)   #works out how far it is from the centre


                countarray[1][r] += 1   #adds 1 to count, r/width so can save larger wedges
                countarray[0][r] += image3Dinv[i][j][k]

    PowerSpectrum = countarray[0]/countarray[1]

    rsize = int(np.sqrt(3.*rmax**2)/width)
    #need an array that represents kmax (to the corner of the cube) = np.sqrt(3*(1/dtheta)**2)
    kmax = np.sqrt(3*(1/dtheta)**2)
    kaxis = np.arange(0,kmax+1,(kmax)/rsize) # rmax steps on the kaxis - ranging from 0 to kmax

    print len(kaxis)
    print len(PowerSpectrum)
    #kaxis = range(0,rmax,width)
    #kaxis = (kaxis/rmax)


    return kaxis, PowerSpectrum







'''
    #this goes over all points in image, and adds the cumulative counts and number of counts for each band of k, starting with < r=1
    for r in range (rmax):
        for i in range (-r,r,1):    #from -(r+1) to r+1 as then sqrt(x^2 + y^2 + z^2) must be < r^2
            for j in range (-r,r,1):
                for k in range (-r,r,1):
                    if np.sqrt(i**2+j**2+k**2) < r+1 and np.sqrt(i**2+j**2+k**2) >= r:    #r+1 as starts at 0 and goes to rmax-1
                        countarray[1][int(r/psdwidth)-1] += 1  #not the -1 is a fudge to get it to stay in bounds, but only determining where data is put
                        countarray[0][int(r/psdwidth)-1] += image3D[rmax + i][rmax + j][rmax + k]
        print r

    return countarray[0]/countarray[1]
    '''

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


def psf(dtheta,image,size):

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


def psfcrosssection(dtheta, image,size):

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
def invert2D(image, dtheta):

    RangeinComplexImage = 0.5/dtheta

    image2 = np.log(abs(image)+1)   #!!!LOGGED!!! to make clear

    #printgraph(image2, RangeinComplexImage,RangeinComplexImage,"1/theta","1/theta",'None','None')

    #image2 = plt.imshow(image2, extent=(-RangeinComplexImage,RangeinComplexImage,-RangeinComplexImage,RangeinComplexImage), interpolation='nearest',cmap='hot')
    #plt.colorbar( orientation='horizontal')
    #plt.xlabel("1./theta")
    #plt.ylabel("1./theta")
    #plt.show()


    #shows inverse image
    imageinv = np.fft.ifft2(image)
    imageinv = abs(imageinv)

    RangeinRealImage = (len(image[1])*dtheta)/2


    #printgraph(imageinv, RangeinRealImage,RangeinRealImage,"theta","theta", 0,70)


    #plt.imshow(imageinv,extent=(-RangeinRealImage,RangeinRealImage,-RangeinRealImage,RangeinRealImage),  interpolation='nearest', cmap='hot')
    #plt.colorbar( orientation='horizontal')
    #plt.xlabel("theta")
    #plt.ylabel("theta")
    #plt.show()

    return imageinv

def rmscalc (twenty1cm,image3D,max):

    squaredcount = 0    #this counts x**2 + y**2 . . . then will average and sqrt at end

    for i in range(max):
        for j in range(max):
            for k in range(max):
                squaredcount += (np.real(image3D[i][j][k]) - twenty1cm[i][j][k])**2     #or should this be real(image)

    squaredcount = squaredcount/(max**3)
    return np.sqrt(squaredcount)

#this function takes ONE image box and represents it as a 2d plot of horizontal average temperature against z
def visualizereionization(image, size, z, theta):

    yvszimage = np.zeros((size,size));


    for t in range(size):
        for i in range(size):
            for j in range(size):
                yvszimage[i][j]=np.average(image[i][j])



    fig = plt.figure()
    a1=fig.add_subplot(1,2,1)
    imgplot = plt.imshow(twenty1[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
    a1.set_title('21cmfast')
    plt.xlabel('X axis in $^\circ$s')
    plt.ylabel('Y axis in $^\circ$s')
    a2=fig.add_subplot(1,2,2)
    plt.setp( a2.get_yticklabels(), visible=False)
    imgplot = plt.imshow(image[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
    plt.xlabel('X axis in $^\circ$s')
    a2.set_title('SKA Image')



    #plt.imshow(yvsz, extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
    #plt.xlabel('X axis in $^\circle$s')
    #plt.ylabel('Y axis in $^\circle$s')
    #plt.colorbar( orientation='vertical')
    plt.savefig('Image/image%i.png'%t)
    plt.close(fig)


#this function prints out all the slices of 2 boxes to be compared - create gif on freds computer using "convert -delay 10 image*.png animated.gif"
def visualizereionizationslicebyslice(image,twenty1, size, z, theta):

    for t in range(size):
        fig = plt.figure()
        a1=fig.add_subplot(1,2,1)
        imgplot = plt.imshow(twenty1[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
        a1.set_title('21cmfast')
        plt.xlabel('X axis in $^\circ$s')
        plt.ylabel('Y axis in $^\circ$s')
        a2=fig.add_subplot(1,2,2)
        plt.setp( a2.get_yticklabels(), visible=False)
        imgplot = plt.imshow(image[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
        plt.xlabel('X axis in $^\circ$s')
        a2.set_title('SKA Image')

        plt.savefig('Image/image%i.png'%t)
        plt.close(fig)





