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
    return zhat,twenty1     #want both as twenty1 will be used in rms calculation


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
                dz.append(((float(positions[i][2])-float(positions[j][2]))/(lda)))

    #float(positions[i][0]) = x coord of ith element
    #float(positions[i][1]) = y coord of ith element

    return (dx,dy,dz)

#this is to show 2D power spectrum evolving with z
def powerspectrumevolution(data,width,size,dtheta):
    rmax = int(size/2.)
    rsize = 1+int(np.sqrt(2.*rmax**2)/width)      #this is largest r box (taking into account shell widths)

    psevolution = np.zeros((size,rsize))    #this will store each 2D powerspectrum with its z

    #xaxis needs to be k, this gives the units/ratio of units
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

##2D PS Calculation##
#starts at rmax, adds all pixel counts for sqrt(x^2 + y^2) < r and number of pixels counted
def powerspectrum2D(data,width,size): #here width means how many pixels wide each band is

    rmax = int(size/2.)  #also same as centre

    #slits of width 'width' and radius is rmax, so longest r is sqrt(r**2 + r**2)
    countarray = np.zeros((2,1 + int(np.sqrt(2.*rmax**2)/width)))

    #goes over all image points, and adds the cumulative counts & number of counts for each band of k
    for i in range(size):
        for j in range (size):
            r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2)/width)   #distance from centre, r/width so can save larger wedges

            countarray[1][r] += 1   #adds 1 to count, used to average
            countarray[0][r] += data[i][j]

    return countarray[0]/countarray[1]

##3D PS Calculation##
#IF YOU WANT WEDGE MASKED, CALL func.EORWINDOW BEFORE SENDING IMAGE INTO THIS FUNCTION
def powerspectrum3D(image3Dinv,width,size,dtheta, dx, z): #here width means how many pixels wide each band is

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

    image3Dinv = np.abs(image3Dinv)**2  #giving us |P(k)|**2

    countarray = np.zeros((3,1+int(np.sqrt(3.*rmax**2)/width))) #storing k wedges of averaged |P(k)|**2

    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2 + (k-rmax)**2)/width)   #works out how far it is from the centre

                #saving both |P(k)|**2 and k**3 |P(k)|**2 (/2pi**2)
                if image3Dinv[i][j][k] > 0:    #this is as masked points are set to -1. No original points can be < 0
                    countarray[1][r] += 1   #adds 1 to count, r/width so can save larger wedges
                    countarray[0][r] += image3Dinv[i][j][k]
                    countarray[2][r] += (r*width*ktorratio)**3 * image3Dinv[i][j][k]  #/(2*np.pi**2)

    PowerSpectrum = countarray[0]/(countarray[1]*size**3)   # have to divide by V    # FUDGE (2*np.pi)**3/(2*np.pi**2)
    DelDel = countarray[2]/(countarray[1]*size**3)

    #need an array that represents kmax (to the corner of the cube) = np.sqrt(3*(1/dtheta)**2)
    rsize = int(np.sqrt(3.*rmax**2)/width)

    #print 'the smallest steps are equal to '
    #print dtheta*kscalefactor

    kaxis = np.arange(0,kmax+(kmax)/rsize,(kmax)/rsize) # rmax steps on the kaxis - ranging from 0 to kmax

    return kaxis, PowerSpectrum, DelDel # delta(k-ko) gives a factor of V - assuming no (2pi)**3 factor - depends on the three dimensional fourier convention -


#Method: takes all of the locations and maps out the UV plane completely together with a count
#We will have to mask it onto the image later.
def rotationmatrix(dx, dy, dz, scaling, H, dH, integrationtime, delta, size):

    ci = int(size/2)
    lenbase = len(dx)

    UVcount = np.zeros((size,size),'complex')

    countvar = 0.   #this will give us the % of baselines ignored (outside image)

    #This measures the UV plane AND maps it onto our fourier plane
    for t in range (integrationtime): #iterates through time

        uvplane = np.zeros((3,lenbase))

        for i in range(lenbase): #one application of the rotation matrix
            uvplane[0][i] = int(ci + scaling*(np.sin(H)*dx[i] + np.cos(H)*dy[i])) # have UV matrix with enough space to get 24 hours of integration
            uvplane[1][i] = int(ci + scaling*(-np.sin(delta)*np.cos(H)*dx[i] + np.sin(delta)*np.sin(H)*dy[i] + np.cos(delta)*dz[i]))
            uvplane[2][i] = np.cos(delta)*np.cos(H)*dx[i] - np.cos(delta)*np.sin(H)*dy[i] + np.sin(delta)*dz[i]

        #Sample fourier space positions (zhat) for this given time
        for i in range(len(uvplane[0])):
            if uvplane[0][i] < size and uvplane[1][i] < size and uvplane[0][i] > 0 and uvplane[1][i] > 0 :   #this is to stop crashing
                UVcount[uvplane[0][i], uvplane[1][i]] += 1
            else:
                countvar += 1.

        #this is the scatter plot without rotation
        #plt.scatter(uvplane[0],uvplane[1])
        #plt.show()

        H += dH

    print ("UV Plane Scan Complete, percentage of baselines ignored", (100*countvar/(integrationtime*lenbase)))
    return (UVcount)

#Method: takes uvplane and takes all non-zero values to 1 to give simple psf
def psf(dtheta,image,size):

    psf = np.zeros((size,size))

    #takes non-zero values to 1 as concerned with distribution in uv not count after rotations
    for i in range (size):
        for j in range (size):
            if image[i][j] != 0:
                psf[i][j]=1

    imageinv = np.fft.ifft2(psf)    #psf is in real space
    imageinv = np.fft.fftshift(imageinv)
    imageinv = abs(imageinv)

    RangeinRealImage = size*dtheta/2.   #rescales the xaxis

    printgraph(imageinv, RangeinRealImage,RangeinRealImage,"theta","theta", 'None', 'None')

#Method: takes 3D fourier space image (central) slice and takes cross section of psf
def psfcrosssection(dtheta, image,size):

    ci=int(size/2)

    psf = np.zeros((size,size))

    #takes non-zero values to 1 as concerned with distribution in uv not count after rotations
    for i in range (size):
        for j in range (size):
            if image[i][j] != 0:
                psf[i][j]=1

    psf = np.fft.ifft2(psf) #psf is in real space
    psf = np.fft.fftshift(psf)
    psf = abs(psf)

    #rescales the xaxis
    RangeinRealImage = size*dtheta/2.
    xaxis=np.arange(-RangeinRealImage, RangeinRealImage, dtheta)

    crosssection =psf[ci]   #takes central part of psf (as centre is of most interst/circular symmetry)
    psfimage=plt.plot(xaxis,crosssection)
    plt.ylabel('PSF Amplitude')
    plt.xlabel('theta')
    plt.show()

#calcuates the rms difference between image and 21cmFAST file
def rmscalc (twenty1cm,image3D,max):

    squaredcount = 0    #counts x**2 + y**2

    for i in range(max):
        for j in range(max):
            for k in range(max):
                squaredcount += (np.real(image3D[i][j][k]) - twenty1cm[i][j][k])**2

    squaredcount = squaredcount/(max**3)    #divided by volume to get average
    return np.sqrt(squaredcount)


###################UNFINISHED FRED###############################
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


#Method: prints out all the slices of 2 boxes to be compared - creates gif on freds computer - "convert -delay 10 image*.png animated.gif"
def visualizereionizationslicebyslice(image,twenty1, size, z, theta):

    for t in range(size):
        fig = plt.figure()
        a1=fig.add_subplot(1,2,1)
        imgplot = plt.imshow(twenty1[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet')#,vmin=0,vmax=70)
        a1.set_title('21cmfast')
        plt.xlabel('X axis in $^\circ$s')
        plt.ylabel('Y axis in $^\circ$s')
        a2=fig.add_subplot(1,2,2)
        plt.setp( a2.get_yticklabels(), visible=False)
        imgplot = plt.imshow(image[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet')#,vmin=0,vmax=70)
        plt.xlabel('X axis in $^\circ$s')
        a2.set_title('SKA Image')

        plt.savefig('Image/z%iimage%i.png'%(z,t))
        plt.close(fig)


#Method: counts and compares the number of corresponding angles (x axis is 21cmfast and y axis is our image)
def phasecomparison(twenty1, image, size):

    phasearray=np.zeros((629,629))

    for i in range(size):
        for j in range(size):
            for k in range(size):
                twenty1phase=int(100*(np.angle(twenty1[i][j][k])+np.pi)) #turns phase onto grid

                imagephase = int(100*(np.angle(image[i][j][k])+np.pi))

                phasearray[twenty1phase][imagephase]+=1

    plt.imshow(phasearray, extent=(0,2*np.pi, 2*np.pi,0))
    plt.xlabel('21cmfast Phase (radians)')
    plt.ylabel('Image Phase (radians)')
    plt.show()

# This is a binning algorithm (linear), might have to move on to a log binning if we want clearer data...
def binningforbubblesizedist(distribution, binsizes):

    distribution = np.trim_zeros(distribution, 'b') #this trims the trailing zeros
    distsize=len(distribution)

    binneddist=np.zeros((2,int(distsize/binsizes)+1))

    for i in range(distsize):
        if i%binsizes == 0:
            binneddist[0][i/binsizes] = float(i + float(binsizes)/2)
        binneddist[1][int(i/binsizes)] = distribution[i]

    return binneddist[0], binneddist[1]

def logbinningforbubblesizedist(distribution, size,dl, powerfactor=1.5):

    distribution = np.trim_zeros(distribution, 'b') #this trims the trailing zeros
    largestbubble = len(distribution)

    numberofbins = int(np.ceil(np.log(largestbubble)/np.log(powerfactor)))

    bins = np.zeros((2, numberofbins))

    #designating bin limits
    for i in range(numberofbins):
        bins[0][i]=powerfactor**i #need to account for loads of bins in range 1-10

    #puts all distribution points into correct bins
    for i in range(numberofbins-1):
        lowtemp = int(np.floor(bins[0][i]))
        hightemp = int(np.floor(bins[0][i+1]))
        for j in range(lowtemp,hightemp):
            bins[1][i]+=distribution[j]

        bins[0][i]=np.sqrt(float(lowtemp*hightemp)) #this gives x value for bin (log position)

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

    for i in range(size):
        for j in range(size):
            for k in range(size):

                #set up the ones and zeros according to some cut-off temperature
                #threshold temperature defined as the required temperature for ionized space
                if image[i][j][k] < thresholdfraction*averagetemp:
                    image[i][j][k] = 0  #0 values will be counted, 1's ignored
                else:
                    image[i][j][k] = 1

    print '%s neutral fraction is' % imagename
    print np.average(image)

    for i in range(size):
        for j in range(size):
            for k in range(size):

                #for all points above threshold it calculates the bubble size into a bubble size distribution
                if image[i][j][k]==0:
                    image, count = numberofnearestneighbours(image, i, j, k, size)
                    distribution[count]+=1

    return logbinningforbubblesizedist(distribution,size,dl, 1.5) #this automatically logbins




#Method: goes through all the unoccupied nearest neighbours, returns the modified image and the  bubble size count
def numberofnearestneighbours(image, i, j, k, size):

    bubblesize=1
    image[i][j][k]=1    #by setting to 1 this means it will be ignored in future, and not recounted

    unvisitedi=list()
    unvisitedj=list()
    unvisitedk=list()

    unvisitedi.append(i)
    unvisitedj.append(j)
    unvisitedk.append(k)

    #keep checking while there are still sites which haven't been checked for occupied nearest neighbours
    while len(unvisitedi) != 0:

        #takes the last co-ordinate from the list
        i=unvisitedi.pop()
        j=unvisitedj.pop()
        k=unvisitedk.pop()

        #These functions check all the nearest neighbours - count them also sets equal to 1
        #pushes unvisited sites onto list to be checked for nearest neighbours
        #checks that nearest neighbours aren't outside of matrix before checking
        if i!= size-1 and image[i+1][j][k] == 0:           # i+1
            image[i+1][j][k]=1  #so it isn't double checked
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

#Method: this calculates the distribution from mean free paths of random points in an ionized region
def secondbubbledistributioncalculator(image,size, thresholdfraction,dl,iterations = 10000):
    #dl is dim[box]/size which give how many MPc each element is

    cutoff = thresholdfraction*np.average(image)     #this is cutoff threshold based on average temp

    print 'temperature cutoff is', cutoff

    rmax = np.sqrt(3*size**2)
    meanfreepathdistribution = np.zeros(rmax)   #this is rmax long as that is largest possible bubble size

    counter = 0

    while counter < iterations:
        #generates a random element of image
        #floor not int as don't want it rounded to size, which would be out of bounds
        i = np.floor(np.random.random()*size)
        j = np.floor(np.random.random()*size)
        k = np.floor(np.random.random()*size)

        if image[i][j][k] < cutoff:
            #randn gives uniform distribution, this is our direction of iteration from element
            di = np.random.randn()
            dj = np.random.randn()
            dk = np.random.randn()

            if di == 0 and dj == 0 and dk == 0:     #safeguard, very unlikely
                di = 1

            length = np.sqrt(di**2 + dj**2 + dk**2) #so we can normalise
            di = di/length  #normailised
            dj = dj/length
            dk = dk/length

            freepath = 1    #counts how far we've moved

            #checking not going to go out of bounds, then iterating while the values are below the cutoff
            if int(i + di) > 0 and int(i + di) < size-1 and int(j + dj) > 0 and int(j + dj) < size-1 and int(k + dk) > 0 and int(k + dk) < size-1:
                while image[int(i+freepath*di)][int(j+freepath*dj)][int(k+freepath*dk)] < cutoff:
                    if int(i + freepath*di) > 0 and int(i + freepath*di) < size-1 and int(j + freepath*dj) > 0 and int(j + freepath*dj) < size-1 and int(k + freepath*dk) > 0 and int(k + freepath*dk) < size-1:
                            freepath += 1
                    else:
                        break   #this stops us going out of image bounds

            #this is adding 1 to the element who's index corresponds to the free path radius
            meanfreepathdistribution[freepath] += 1
            counter += 1



    meanfreepathdistribution = np.trim_zeros(meanfreepathdistribution, 'b') #this trims the trailing zeros

    arraysize = len(meanfreepathdistribution)
    volume = np.zeros(arraysize)  #first column is xposition, second is that positions count

    for i in range (arraysize):
        #converts to volume, but r-0.5 as we don't know how far between ionised and non-ionised point,
        #therefore estimated at halfway. FACTOR of dl**3 as want in MPc^3
        volume[i] = dl**3*4/3*np.pi*(float(i)-0.5)**3

    #this is a normalised weighted mean free path probability, so divided by r^3 to get relation to chance of sampling that size of bubble. divided by N to get probability
    #meanfreepathdistribution = meanfreepathdistribution/volume  #has to be done outside of return so that sum() function works
    return volume, meanfreepathdistribution/sum(meanfreepathdistribution)   #no need for 1/N factor as taken into account with sum

#Method: Takes image and returns masked window image, including replacing centre after masking to avoid flipped image
def EORWINDOW(Windowedimage, size, dl,z,B): #units of r are in terms of the index while k is per Mpc

    CosmoUnits=Cosmo.CosmoUnits()
    #FROM The Epoch of Reionization Window I: Mathematical Formalism by A. Liu et al.
    E = np.sqrt(CosmoUnits.omega_m*((1+z)**3)+CosmoUnits.omega_l)
    Dz = CosmoUnits.Dcomovingrad(z)
    H0 = CosmoUnits.H0 # in km / s /Mpc
    c = CosmoUnits.C
    theta0 = 0.000070 #this is 4 degrees in radians (ther is a 'conservative' of 1 radian according to the paper which seems way too large)
    #NEED TO TALK TO PRITCHARD ABOUT THIS!!! WHAT IS THE CHARACTERISTIC BEAM THICKNESS!!!?


    kmax = np.sqrt(3*(1/(2*float(dl)))**2) #finds kmax (from centre to outer corner) for a k space that is 1/2dx large
    rspacemaxradius=np.sqrt(3*(size/2)**2)
    ktorratio=kmax/rspacemaxradius #ratio between our indexes and kspace
    rmax = int(size/2.)


    parrkcutoff = 0.3#/ktorratio #(H0/(((1+z)**2)*B))/ktorratio #the size of contaminated kparrallels with no wedge
    print parrkcutoff


    centre = np.zeros((125),dtype=complex)
    counter = 0

    #this is saving the centre cube, so it is safe from masking
    for i in range (rmax -2,rmax +3,1):
        for j in range (rmax -2,rmax +3,1):
            for k in range (rmax -2,rmax +3,1):
                centre[counter] = Windowedimage[i][j][k]
                print Windowedimage[i][j][k]
                counter += 1
    print centre


    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-rmax)**2 + (j-rmax)**2 + (k-rmax)**2))   #works out how far it is from the centre
                kperp = np.sqrt((i-rmax)**2 + (j-rmax)**2)

                if np.abs(k-rmax) < parrkcutoff:
                    Windowedimage[i][j][k]=0.+0.j
                elif np.abs(k-rmax) < kperp*H0*Dz*E*theta0/(c*(1+z)*ktorratio):    #abs as kparrallel = abs(kz)
                    Windowedimage[i][j][k]=0.+0.j

    return Windowedimage

##########################################printing function################################################

# This function compares the powerspectra of the image, the twenty1cmsignal and the error
# realps refers to the powerspectrum as provided by 21cmfast and can be uncommented to compare our results to this
def printpowerspectrum(image3Dinverse, sigma3Dinverse, twenty1inverse,psdwidth,size,dtheta, dx, z):

    #realps = np.loadtxt('ps_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_100_200Mpc_v2.txt', delimiter='\t')

    imagek, imagepowerspectrum , imagedeldel= powerspectrum3Dwedge(image3Dinverse,psdwidth,size,dtheta,dx, z) # this is the size of steps in real space dx=float(box_info['dim'])/float(box_info['BoxSize'])
    print 'done imagepowerspectrum'
    sigmak, sigmapowerspectrum , sigmadeldel= powerspectrum3Dwedge(sigma3Dinverse,psdwidth,size,dtheta, dx, z)
    print 'done sigmapowerspectrum'
    twenty1k, twenty1powerspectrum, twenty1deldel= powerspectrum3D(twenty1inverse,psdwidth,size,dtheta,dx, z)
    print 'done twenty1powerspectrum'

    #plots the compared powerspectra
    plt.loglog(imagek,imagedeldel)
    plt.loglog(sigmak,sigmadeldel)
    plt.loglog(twenty1k,twenty1deldel)
    #plt.loglog(realps[:,0],realps[:,1])
    plt.ylim(0.00001,100)
    plt.xlim(0.02,3)

    plt.xlabel('k (MPc$^{-1}$)')
    plt.ylabel('k$^3$ P(k)/2$\pi^2$')
    plt.savefig('DELDEL POWERSPEC for z = %i' %z)
    plt.clf()

    plt.loglog(imagek,imagepowerspectrum)
    plt.loglog(sigmak,sigmapowerspectrum)
    plt.loglog(twenty1k,twenty1powerspectrum)
    #plt.loglog(realps[:,0],realps[:,1]/(realps[:,0]**3)) #Important: here, as with other plot, we have no 2Pi**2 factor
    plt.ylim(0.02,100000)
    plt.xlim(0.02,3)

    plt.xlabel('k (MPc$^{-1}$)')
    plt.ylabel('P(k)')
    plt.savefig('POWERSPEC for z = %i' %z)
    plt.clf()

def printbubblesizedist(image3D, twenty1, size, dl, cutoff):

    distimagex, distimagey= bubblesizedistribution(image3D, size,dl,cutoff,'image')
    dist21x, dist21y= bubblesizedistribution(twenty1, size,dl,cutoff,'twenty1')

    figure = plt.loglog(distimagex,distimagey)
    plt.loglog(dist21x, dist21y)

    plt.xlabel('Ionised Volume (MPc$^{3}$)')
    plt.ylabel('Probability')
    plt.xlim(2,100000)
    plt.ylim(0.000007,1)
    plt.show()


def printmeanfreepathdist(image3D, twenty1, size, dl, cutoff, iterations):

    imagemeanpathx,imagemeanpathdist = secondbubbledistributioncalculator(image3D,size,cutoff,dl,iterations)
    twenty1meanpathx,twenty1meanpathdist = secondbubbledistributioncalculator(twenty1,size,cutoff,dl,iterations)


    figure= plt.loglog(imagemeanpathx,imagemeanpathdist)
    plt.loglog(twenty1meanpathx,twenty1meanpathdist)

    plt.xlabel('Ionised Volume (MPc$^{3}$)')
    plt.ylabel('Probability')
    plt.xlim(2,100000)
    plt.ylim(0.000007,1)
    plt.show()

