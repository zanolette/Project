
import numpy as np
import pylab as pl
from matplotlib import ticker
import matplotlib.pyplot as plt
import boxio as boximport
import Cosmo as Cosmo

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
def powerspectrum3D(image3Dinv,width,size,dtheta, dx, z): #here width means how many pixels wide each band is


    #UNSURE ABOUT K to R MAPPING
    #CosmoUnits=Cosmo.CosmoUnits() # REMOVE THIS FROM FUNCTIONS LATER MAYBE
    # to get scalefactor of inverse space in MPc^-1 instead of theta^1
    # Ignore - READ BELOW - THIS MIGHT BE AN ISSUE BUT WE ARE CONFIDENT... FOR NOW
    #kmax = np.sqrt(3*(1/(dtheta*kscalefactor)**2))#NO actually maybe the z component is different - we are investigating it in terms of distance already.
    #kmax = np.sqrt(2*(1/(dtheta*kscalefactor)**2)+(1/dx)**2)

    # RECAP - WE GET A BOX THAT HAS REAL SPACE SMALL UNITS OF DX (MPC), HENCE K SPACE LARGEST SIZE IS 1/DX (MPC^-1) - HENCE THE FOLLOWING FOR KMAX:
    kmax = np.sqrt(3*(1/(2*float(dx)))**2) #finds kmax (from centre to outer corner) for a k space that is 1/2dx large

    rspacemaxradius=np.sqrt(3*(size/2)**2)
    ktorratio=kmax/rspacemaxradius

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


    countarray = np.zeros((3,1+int(np.sqrt(3.*rmax**2)/width)))

    print len(countarray[0])

    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2 + (k-rmax)**2)/width)   #works out how far it is from the centre

                if image3Dinv[i][j][k] >= 0:    #this is as masked points are set to -1. No original points can be < 0
                    countarray[1][r] += 1   #adds 1 to count, r/width so can save larger wedges
                    countarray[0][r] += image3Dinv[i][j][k]
                    countarray[2][r] += (r*width*ktorratio)**3 * image3Dinv[i][j][k] #/(2*np.pi**2)

    PowerSpectrum = countarray[0]/(countarray[1]*size**3) # FUDGE (2*np.pi)**3/(2*np.pi**2)  # have to devide by V ???
    DelDel = countarray[2]/(countarray[1]*size**3)


    rsize = int(np.sqrt(3.*rmax**2)/width)

    #need an array that represents kmax (to the corner of the cube) = np.sqrt(3*(1/dtheta)**2)



    #print 'the smallest steps are equal to '
    #print dtheta*kscalefactor

    kaxis = np.arange(0,kmax+(kmax)/rsize,(kmax)/rsize) # rmax steps on the kaxis - ranging from 0 to kmax
    #kaxis=kaxis/kscalefactor

    return kaxis, PowerSpectrum, DelDel # delta(k-ko) gives a factor of V - assuming no (2pi)**3 factor - depends on the three dimensional fourier convention -


##################3D version WITH WEDGE REMOVED#############################
def powerspectrum3Dwedge(image3Dinv,width,size,dtheta, dx, z): #here width means how many pixels wide each band is
    #FOR THIS METHOD: we will first go over all of 3D box and set to zero any elements within wedge, then method will continue as usual


    #UNSURE ABOUT K to R MAPPING
    #CosmoUnits=Cosmo.CosmoUnits() # REMOVE THIS FROM FUNCTIONS LATER MAYBE
    # to get scalefactor of inverse space in MPc^-1 instead of theta^1
    # Ignore - READ BELOW - THIS MIGHT BE AN ISSUE BUT WE ARE CONFIDENT... FOR NOW
    #kmax = np.sqrt(3*(1/(dtheta*kscalefactor)**2))#NO actually maybe the z component is different - we are investigating it in terms of distance already.
    #kmax = np.sqrt(2*(1/(dtheta*kscalefactor)**2)+(1/dx)**2)



    # RECAP - WE GET A BOX THAT HAS REAL SPACE SMALL UNITS OF DX (MPC), HENCE K SPACE LARGEST SIZE IS 1/DX (MPC^-1) - HENCE THE FOLLOWING FOR KMAX:
    kmax = np.sqrt(3*(1/(2*float(dx)))**2) #finds kmax (from centre to outer corner) for a k space that is 1/2dx large

    print 'kmax is', kmax

    rspacemaxradius=np.sqrt(3*(size/2)**2)
    ktorratio=kmax/rspacemaxradius

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


    countarray = np.zeros((3,1+int(np.sqrt(3.*rmax**2)/width)))

    print len(countarray[0])

    #these will have an equation eventually
    #we want to write these in units of "index" rather than per Mpc
    #as its easier to work with our image box in terms of indices
    gradient = 0.5
    perpkcutoff = 0.055/ktorratio    #this is at what point the wedge becomes important - note i've taken these initial estimates from wedge I paper
    parrkcutoff = 0.5 /ktorratio  #this is the size of contaminated kparrallels when not needing wedge
    #constant = perpkcutoff - parrkcutoff/gradient   #(this is just working out C for y = mx + C)
    print 'perpkcutoff is', perpkcutoff

    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2 + (k-rmax)**2)/width)   #works out how far it is from the centre
                kperp = np.sqrt((i-rmax)**2 + (j-rmax)**2)

                #condition for < perpkcutoff is k < parrkcutoff
                # if i=x,j=y,k=z, then x**2 + y**2 = perp**2
                #so condition once past perpkcutoff is if k > gradient*sqrt(i**2 + j**2) + constant
                if kperp < perpkcutoff:
                    if np.abs(k-rmax) > parrkcutoff:
                        countarray[1][r] += 1   #adds 1 to count, r/width so can save larger wedges
                        countarray[0][r] += image3Dinv[i][j][k]
                        countarray[2][r] += (r*width*ktorratio)**3 * image3Dinv[i][j][k] #/(2*np.pi**2)

                elif np.abs(k-rmax) > np.sqrt(kperp-perpkcutoff)+parrkcutoff:    #this neglects wedge elements, which is parabola when not in loglog, abs as kparrallel = abs(kz)
                    countarray[1][r] += 1
                    countarray[0][r] += image3Dinv[i][j][k]
                    countarray[2][r] += (r*width*ktorratio)**3 * image3Dinv[i][j][k] #/(2*np.pi**2)


    PowerSpectrum = countarray[0]/(countarray[1]*size**3) # FUDGE (2*np.pi)**3/(2*np.pi**2)  # have to devide by V ???
    DelDel = countarray[2]/(countarray[1]*size**3)


    rsize = int(np.sqrt(3.*rmax**2)/width)

    #need an array that represents kmax (to the corner of the cube) = np.sqrt(3*(1/dtheta)**2)



    #print 'the smallest steps are equal to '
    #print dtheta*kscalefactor

    kaxis = np.arange(0,kmax+(kmax)/rsize,(kmax)/rsize) # rmax steps on the kaxis - ranging from 0 to kmax
    #kaxis=kaxis/kscalefactor

    return kaxis, PowerSpectrum, DelDel # delta(k-ko) gives a factor of V - assuming no (2pi)**3 factor - depends on the three dimensional fourier convention -





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
def visualizereionizationagainstz(image, size, z, theta):

    yvszimage = np.zeros((size,size));


    for t in range(size):
        for i in range(size):
            for j in range(size):
                yvszimage[i][j]=np.average(image[i][j])

    a1=fig.add_subplot(1,2,1)
    imgplot = plt.imshow(twenty1[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
    a1.set_title('21cmfast')
    plt.xlabel('Z axis in Mpc')
    plt.ylabel('Average over Y axis in Mpc')
    a2=fig.add_subplot(1,2,2)
    plt.setp( a2.get_yticklabels(), visible=False)
    imgplot = plt.imshow(image[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
    plt.xlabel('Z axis in Mpc')
    a2.set_title('SKA Image')

    plt.show()


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


#counts and compares the number of corresponding angles (x axis is 21cmfast and y axis is our image)
def phasecomparison(twenty1, image, size):

    phasearray=np.zeros((629,629))

    for i in range(size):
        for j in range(size):
            for k in range(size):
                twenty1phase=int(100*(np.angle(twenty1[i][j][k])+np.pi))

                imagephase = int(100*(np.angle(image[i][j][k])+np.pi))

                phasearray[twenty1phase][imagephase]+=1

    plt.imshow(phasearray, extent=(0,2*np.pi, 2*np.pi,0))
    plt.xlabel('21cmfast Phase (radians)')
    plt.ylabel('Image Phase (radians)')
    plt.show()

# This is a binning algo - linear binning - might have to move on to a log binning if we want clearer data...
def binningforbubblesizedist(distribution, binsizes):

    distribution = np.trim_zeros(distribution, 'b') #this trims the trailing zeros

    distsize=len(distribution)
    print distsize

    binneddist=np.zeros((2,int(distsize/binsizes)+1))

    for i in range(distsize):
        if i%binsizes == 0:
            binneddist[0][i/binsizes] = float(i + float(binsizes)/2)
        binneddist[1][int(i/binsizes)] = distribution[i]

    return binneddist[0], binneddist[1]

def logbinningforbubblesizedist(distribution, size,dl, powerfactor=1.5):

    distribution = np.trim_zeros(distribution, 'b') #this trims the trailing zeros
    largestbubble = len(distribution)

    numberofbins = int(np.ceil(np.log(largestbubble)/np.log(powerfactor))) #safetyfactor

    bins = np.zeros((2, numberofbins))

    for i in range(numberofbins):
        bins[0][i]=powerfactor**i #need to account for loads of bins in range 1-10

    for i in range(numberofbins-1):
        lowtemp = int(np.floor(bins[0][i]))
        hightemp = int(np.floor(bins[0][i+1]))
        for j in range(lowtemp,hightemp):
            bins[1][i]+=distribution[j]

        bins[0][i]=np.sqrt(float(lowtemp*hightemp))
        #bin[0][i]=np.sqrt(bin[0][i]*bin[0][i+1]) # relabelling the correct log position of each bin

    numberofcounts = 0  #this is so we know the normalisation factor
    for i in range(len(bins[1])):
        numberofcounts += bins[1][i]


    return dl**3*bins[0], bins[1]/numberofcounts    # factor of dl**3 as in volume (element^3) but need MPc^3




def bubblesizedistribution(imageoriginal, size,dl,thresholdfraction,imagename):

    image = imageoriginal #we dont want to change the original

    distribution=np.zeros(size*size*size) # this is the maximum size of a bubble so thats what ive made the distribution array equal to

    averagetemp=np.average(image)

    print 'printing average temp'
    print averagetemp
    temperature=0

    for i in range(size):
        for j in range(size):
            for k in range(size):

                #set up the ones and zeros according to some cut-off temperature
                if image[i][j][k] < thresholdfraction*averagetemp: # this is checking is a point is below a threshold temperature that we define as the temperature required to describe a space ionized
                    image[i][j][k] = 0
                else:
                    image[i][j][k] = 1

    print '%s neutral fraction is' %imagename
    print np.average(image)

    for i in range(size):
        for j in range(size):
            for k in range(size):

                if image[i][j][k]==0:
                    image, count = numberofnearestneighbours(image, i, j, k, size)
                    distribution[count]+=1


    return logbinningforbubblesizedist(distribution,size,dl, 1.5) #this automatically logbins




#this function goes through all the unoccupied nearest neighbours and returns the modified image and the count of that bubble size
def numberofnearestneighbours(image, i, j, k, size):

    bubblesize=1
    image[i][j][k]=1

    unvisitedi=list()
    unvisitedj=list()
    unvisitedk=list()

    unvisitedi.append(i)
    unvisitedj.append(j)
    unvisitedk.append(k)

    while len(unvisitedi) != 0: # while there are still sights which havent been checked for occupied nearest neighbours - keep checking

        i=unvisitedi.pop()
        j=unvisitedj.pop()
        k=unvisitedk.pop()

        #These functions check all the nearest neightbours - count them also sets equal to 1
        #also pushes onto unvisited

        if i!= size-1 and image[i+1][j][k] == 0:           # i+1
            image[i+1][j][k]=1
            bubblesize+=1
            unvisitedi.append(i+1)
            unvisitedj.append(j)
            unvisitedk.append(k)
        if  i!=0 and image[i-1][j][k] == 0:           # i-1
            image[i-1][j][k]=1
            bubblesize+=1
            unvisitedi.append(i-1)
            unvisitedj.append(j)
            unvisitedk.append(k)
        if j!= size-1 and image[i][j+1][k] == 0:           # j+1
            image[i][j+1][k]=1
            bubblesize+=1
            unvisitedi.append(i)
            unvisitedj.append(j+1)
            unvisitedk.append(k)
        if j!=0 and image[i][j-1][k] == 0:           # j-1
            image[i][j-1][k]=1
            bubblesize+=1
            unvisitedi.append(i)
            unvisitedj.append(j-1)
            unvisitedk.append(k)
        if k!= size-1 and image[i][j][k+1] == 0:           # k+1
            image[i][j][k+1]=1
            bubblesize+=1
            unvisitedi.append(i)
            unvisitedj.append(j)
            unvisitedk.append(k+1)
        if k!= 0 and image[i][j][k-1] == 0:           # k-1
            image[i][j][k-1]=1
            bubblesize+=1
            unvisitedi.append(i)
            unvisitedj.append(j)
            unvisitedk.append(k-1)

    return image, bubblesize

def secondbubbledistributioncalculator(image,size, thresholdfraction,dl,iterations = 10000):   #dl is dim[box]/size which give how many MPc each element is

    cutoff = thresholdfraction*np.average(image)     #this is cutoff threshold based on average temp

    print cutoff

    rmax = np.sqrt(3*size**2)
    meanfreepathdistribution = np.zeros(rmax)   #this is rmax long as that is largest possible bubble size

    counter = 0

    while counter < iterations:
        #random element of image
        i = np.floor(np.random.random()*size) #floor not int as don't want it rounded to size, which would be out of bounds
        j = np.floor(np.random.random()*size)
        k = np.floor(np.random.random()*size)

        if image[i][j][k] < cutoff:

            di = np.random.randn() #floor not int as don't want it rounded to size, which would be out of bounds
            dj = np.random.randn()
            dk = np.random.randn()

            if di == 0 and dj == 0 and dk == 0:     #safeguard
                di = 1

            length = np.sqrt(di**2 + dj**2 + dk**2) #so we can normalise

            di = di/length  #normailised
            dj = dj/length
            dk = dk/length

            freepath = 1

            if int(i + di) > 0 and int(i + di) < size-1 and int(j + dj) > 0 and int(j + dj) < size-1 and int(k + dk) > 0 and int(k + dk) < size-1:  #checking already not going to go out of bounds
                while image[int(i+freepath*di)][int(j+freepath*dj)][int(k+freepath*dk)] < cutoff:
                    if int(i + freepath*di) > 0 and int(i + freepath*di) < size-1 and int(j + freepath*dj) > 0 and int(j + freepath*dj) < size-1 and int(k + freepath*dk) > 0 and int(k + freepath*dk) < size-1:
                            freepath += 1
                    else:
                        break   #this stops us going out of image bounds

            meanfreepathdistribution[freepath] += 1   #this is adding 1 to the element which's index corresponds to the free path radius - ie estimate of radius

            counter += 1



    meanfreepathdistribution = np.trim_zeros(meanfreepathdistribution, 'b') #this trims the trailing zeros

    arraysize = len(meanfreepathdistribution)
    volume = np.zeros(arraysize)  #first column is xposition, second is count at that x, we'll rescale xpositions to 4pi * r^3

    for i in range (arraysize):
        volume[i] = dl**3*4/3*np.pi*(float(i)-0.5)**3     #this converts to volume, but r-0.5 as last point we don't know how far between ionised and notionised point, so guess halfway. FACTOR of dl**3: as want in MPc^3

    #this is a normalised weighted mean free path probability, so divided by r^3 to get relation to chance of sampling that size of bubble. divided by N to get probability
    #meanfreepathdistribution = meanfreepathdistribution/volume  #has to be done outside of return so that sum() function works
    return volume, meanfreepathdistribution/sum(meanfreepathdistribution)   #no need for 1/N factor as taken into account with sum
