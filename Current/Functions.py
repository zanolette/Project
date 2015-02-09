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

#######
#do we still need this?
#######
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

#this takes a 3D fourier space, plots k parrallel vs k perpendicular (remember kperp = sqrt(kx**2 + ky**2),kparr = kz)
def kperpvskparrgraph(data,width,size,dl):

    length = 1 + int(np.sqrt(.5*size**2/width))
    print length
    dataarray = np.zeros((length,size))  #creates our kperp vs kparr values, value at element is the average PS at that point

    for i in range (size):
        dataarray[:][i] = powerspectrum2D(data[:][:][i], width, size)   #this returns the 1D array of averaged k perpendicular values, for that kparr

    #want axis not in element steps, but in k steps (both axis), but both axis different lengths
    #rparrsize = size    #or could be size/2. if we want abs(kparr) . . .
    kparrmax = 1./dl
    #kparraxis = np.arange(0,kparrmax+kparrmax/rparrsize,kparrmax/rparrsize) # rmax steps on the kparrallelaxis

    #rperpsize = int(np.sqrt(.5*size**2/width))  #or + 1??
    kperpmax = np.sqrt(2*((1./(2*float(dl)))**2))


    #kperpaxis = np.arange(0,kperpmax+kperpmax/rperpsize,kperpmax/rperpsize) # rmax steps on the kperpendicularaxis

    fig = plt.figure()

    ax1 = fig.add_subplot(111) # x (z) and y axis

    ax1.imshow(dataarray,extent=(0,kperpmax,0,1/(2*dl)), cmap='jet',interpolation='nearest',)    #is this right?
    ax1.set_xlabel('k Parrallel (MPc$^{-1}$)')
    ax1.set_ylabel('k Perpendicular (MPc$^{-1}$)')
    #resize axis how?

    #########This commented section is the print NF thing from below copied to help maybe
    '''
    nf_axis_ticklocations = np.array([1, 5, 9, 13, 17]) # in terms of z
    nfindexes=np.searchsorted(redshift, nf_axis_ticklocations) # finds corresponding indices

    ax2 = ax1.twiny() # further x axis corresponding to the same y axis
    ax2.set_xticks(nf_axis_ticklocations) # ticks at desires z locations.
    ax2.set_xticklabels(neutralfractions[nfindexes]) # prints nf for each z location
    ax2.set_xlabel("Un-ionized Fraction")
    '''




    #plt.savefig('KperpvsKparr')
    plt.show()


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
            #do we want deldel??
    return countarray[0]/countarray[1]

##3D PS Calculation##
#IF YOU WANT WEDGE MASKED, CALL func.EORWINDOW BEFORE SENDING IMAGE INTO THIS FUNCTION
def powerspectrum3D(fourier3Dbox,width,size,dtheta, dl, z,cornercorrection): #here width means how many pixels wide each band is

    #cornercorretion == 1 means that we do not count the outermost corners of k space to our powerspectrum
    # this is because there a very few points and statistics are poor.

    # RECAP - WE GET A BOX THAT HAS REAL SPACE SMALL UNITS OF DX (MPC), HENCE K SPACE LARGEST SIZE IS 1/DX (MPC^-1) - HENCE THE FOLLOWING FOR KMAX:
    kmax = np.sqrt(3*((1./(2*float(dl)))**2)) #finds kmax (from centre to outer corner) for a k space that is 1/2dl large

    ktorratio=1./(dl*size) # ratio between k space and indexes (r space in our minds) in per mega parsecs
    #print 'ktorratio', ktorratio

    ci = int(size/2.)  #centre point

    '''
    plt.imshow(image3D[50],  interpolation='nearest',cmap='jet')
    plt.colorbar( orientation='vertical')
    plt.show()
    '''

    absimage3Dinv = np.abs(fourier3Dbox)**2  #giving us |d(k)|**2

    #need an array that represents kmax (to the corner of the cube) = np.sqrt(3*(1/dtheta)**2)
    rsize = int(np.sqrt(3.*ci**2)/width)    #rsize is number of integer steps to get to kmax (ie if width isn't 1)
    countarray = np.zeros((3,1+rsize)) #storing k wedges of averaged |d(k)|**2

    for i in range(size):
        for j in range (size):
            for k in range (size):

                r = int(np.sqrt((i-ci)**2 + (j-ci)**2 + (k-ci)**2)/width)   #works out how far it is from the centre

                #saving both |P(k)|**2 and k**3 |P(k)|**2 (/2pi**2)
                if absimage3Dinv[i][j][k] != 0:    #Do no measure deleted or unsampled points.
                    countarray[1][r] += 1   #adds 1 to count, r/width so can save larger wedges
                    countarray[0][r] += absimage3Dinv[i][j][k]
                    countarray[2][r] += (r*width*ktorratio)**3 * absimage3Dinv[i][j][k]


    PowerSpectrum = countarray[0]/(countarray[1]*size**3)   # have to divide by V    #  (2*np.pi)**3/(2*np.pi**2)
    DelDel = countarray[2]/(countarray[1]*size**3)

    #print 'the smallest steps are equal to '
    #print dtheta*kscalefactor

    kaxis = np.arange(0,kmax+(kmax)/rsize,(kmax)/rsize) # rmax steps on the kaxis - ranging from 0 to kmax

    # removes the highest k points from powerspectrum if cornercorrection == 1, this is due to their statistical inaccuracy

    if cornercorrection == 1:
        kaxis = np.delete(kaxis, (len(kaxis)-1))
        PowerSpectrum = np.delete(PowerSpectrum, (len(PowerSpectrum)-1))
        DelDel = np.delete(DelDel, (len(DelDel)-1))

    return logbinningforPowerspectrum(kaxis, PowerSpectrum, DelDel, dl, 1.5) # delta(k-ko) gives a factor of V - assuming no (2pi)**3 factor - depends on the three dimensional fourier convention -


def logbinningforPowerspectrum(K, PofK, DelDel, dl, powerfactor=1.2):

    binlimit = 0.01 #this is first bin limit in k space
    kbinssize = np.array(binlimit) #this will hold the bin edges

    numberofbins = 0

    while binlimit < np.amax(K):
        if binlimit < 0.99:
            binlimit = binlimit**(1/(1.1)) # change bin limit logrithmically
            kbinssize = np.append(kbinssize, binlimit) #saves new bin limit
            numberofbins+=1
        elif binlimit <= 1: # jumps accross 1 for log asymptote
            binlimit = 1.01
            kbinssize = np.append(kbinssize, binlimit) #saves new bin limit
            numberofbins+=1
        else:
            binlimit = binlimit**(powerfactor) # change bin limit logrithmically
            kbinssize = np.append(kbinssize, binlimit) #saves new bin limit
            numberofbins+=1

    kbinssize= np.insert(kbinssize, 0, 0.) # insert 0 point.

    PofKBins = np.zeros(numberofbins+1) #stores new summed P(k) # needs one more as we have inserted 0.
    DelDelBins = np.zeros(numberofbins+1)  ##stores new summed DelDel

    i=0 #this is the index for K
    j=0
    countsperbin=0
    while j < len(kbinssize)-1 and i < len(K):   #goes over each bin
        if K[i]>=kbinssize[j] and K[i]<kbinssize[j+1]: #if K[i] is in bin, put it in
            #its never here!!!
            PofKBins[j]+=PofK[i]
            DelDelBins[j]+=DelDel[i]
            i += 1  #advance i so that we can look at next k
            countsperbin+=1
        else:
            kbinssize[j] = np.sqrt(kbinssize[j]*kbinssize[j+1])    #when we move onto the next bin, defines new intermediate k value
            if countsperbin != 0:
                PofKBins[j] = PofKBins[j]/countsperbin # find average value in a bin once finished adding to that bin.
                DelDelBins[j] = DelDelBins[j]/countsperbin
                countsperbin=0 #reset the number of counts for an added bin.
            j += 1 # move to next bin

    k=0 #k is the integer counter
    while k < (len(PofKBins)):
        if PofKBins[k] == 0.:
            kbinssize=np.delete(kbinssize, k)  #if PofK[k] is zero then want to delete all 3 elements, leaving only non-0 values
            PofKBins=np.delete(PofKBins, k)
            DelDelBins=np.delete(DelDelBins, k)
        else:
            k += 1  #if not 0 then looks at the next PofK[k]

    kbinssize=np.delete(kbinssize, (np.size(kbinssize)-1)) # delete the top box which doesnt refer to anything
    return kbinssize, PofKBins, DelDelBins  # factor of dl**3 as in volume (element^3) but need MPc^3



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
    tempcounter = 0 #to get percentage of UV coverage
    for i in range(size):
        for j in range(size):
            if UVcount[i][j] != 0:
                tempcounter += 1
    tempcounter = tempcounter/(size**3) #gives percentage

    print ("UV Plane Scan Complete, percentage of baselines ignored", (100*countvar/(integrationtime*lenbase)))
    return (UVcount), tempcounter

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

#calcuates the rms  between each pixel away from average value image and 21cmFAST file - dimensionless
def rmscalc (twenty1cm,image3D,max):

    squaredcount = 0    #counts x**2 + y**2
    twenty1average = np.average(twenty1cm)  #these are our average values
    twenty1var = np.var(twenty1cm)
    imagevar = np.var(image3D)
    imageaverage = np.average(image3D)

    for i in range(max):
        for j in range(max):
            for k in range(max):
                squaredcount += ((image3D[i][j][k]-imageaverage)/imagevar - (twenty1cm[i][j][k]-twenty1average)/twenty1var)**2

    squaredcount = squaredcount/(max**3)    #divided by volume to get average
    return np.sqrt(squaredcount)

def PearsonR(twenty1cm,image3D,max):

    #Pearson's R is a correlation indicator also known as a sample correlation coefficient
    #the equation is sum((x-<x>)(y-<y>))/(sqrt(sum((x-<x>)^2))*sqrt(sum((y-<y>)^2)))

    #sum((x-<x>)(y-<y>)) =
    Covariance = 0
    #sum((x-<x>)^2) =
    twenty1selfcovariance = 0
    #sum((y-<y>)^2) =
    imageselfcovariance = 0

    twenty1average = np.average(twenty1cm)  #these are our average values
    imageaverage = np.average(image3D)

    for i in range(max):
        for j in range(max):
            for k in range(max):
                Covariance += (image3D[i][j][k]-imageaverage)*(twenty1cm[i][j][k]-twenty1average)
                twenty1selfcovariance += (twenty1cm[i][j][k]-twenty1average)**2
                imageselfcovariance += (image3D[i][j][k]-imageaverage)**2

    PearsonsR = Covariance/(np.sqrt(twenty1selfcovariance)*np.sqrt(imageselfcovariance))

    return PearsonsR
'''
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
'''

#Method: prints out all the slices of 2 boxes to be compared - creates gif on freds computer - "convert -delay 10 image*.png animated.gif"
def visualisereionizationslicebyslice(image,twenty1, size, z, theta):

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

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.3, 0.03, 0.4])
        cbar=plt.colorbar(imgplot, orientation='vertical',cax=cbar_ax)# shrink=0.6)

        cbar.set_label('Temperature (mK)', size = 13)

        plt.savefig('Image/z%simage%03i.png'%(z,t))
        plt.close(fig)


#Method: counts and compares the number of corresponding angles (x axis is 21cmfast and y axis is our image)
def phasecomparison(twenty1, image, size):

    phasearray=np.zeros((629,629))

    for i in range(size):
        for j in range(size):
            for k in range(size):
                twenty1phase=int(100*(np.angle(twenty1[i][j][k]))) #turns phase onto grid

                imagephase = int(100*(np.angle(image[i][j][k])))

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

    #print 'temperature cutoff is', cutoff

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
            #if counter%1000==0:
             #   print counter



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
def EORWINDOW(Windowedimageinv, size, dl,z,B): #units of r are in terms of the index while k is per Mpc

    CosmoUnits=Cosmo.CosmoUnits()
    #FROM The Epoch of Reionization Window I: Mathematical Formalism by A. Liu et al.
    E = np.sqrt(CosmoUnits.omega_m*((1+z)**3)+CosmoUnits.omega_l)
    Dz = CosmoUnits.Dcomovingrad(z)
    H0 = CosmoUnits.H0 # in km / s /Mpc
    c = CosmoUnits.C
    theta0 = 0.000070 #this is 4 degrees in radians (ther is a 'conservative' of 1 radian according to the paper which seems way too large)
    nu21 = 1420000000 # Hz
    Bband = 8000000 # Hz
    #NEED TO TALK TO PRITCHARD ABOUT THIS!!! WHAT IS THE CHARACTERISTIC BEAM THICKNESS!!!?


    kmax = np.sqrt(3*(1/(2*float(dl)))**2) #finds kmax (from centre to outer corner) for a k space that is 1/2dl large
    rspacemaxradius=np.sqrt(3*(size/2)**2)
    ktorratio=kmax/rspacemaxradius #ratio between our indexes and kspace
    print 'ktorratio', ktorratio
    ci = int(size/2.)


    #CUTOFF SHOULD BE IN INDEX UNITS AS THIS IS HOW THE CONDITIONAL IS USED LATER
    parrkcutoff = (H0*2*np.pi*E*nu21/(((1+z)**2)*Bband*c))/ktorratio #the size of contaminated kparrallels with no wedge

    #print 'parkis', parrkcutoff

    Niquist = Windowedimageinv[ci][ci][ci]
    '''
    centre = np.zeros((27),dtype=complex)
    counter = 0

    #this is saving the centre cube, so it is safe from masking

    for i in range (ci -1,ci +2,1):
        for j in range (ci -1,ci +2,1):
            for k in range (ci -1,ci +2,1):
                centre[counter] = Windowedimageinv[i][j][k]
                #print Windowedimageinv[i][j][k]

                counter += 1
    '''

    for i in range(size):
        for j in range (size):
            for k in range (size):

                kperp = np.sqrt((i-ci)**2 + (j-ci)**2)

                if np.abs(k-ci) < parrkcutoff:
                    Windowedimageinv[i][j][k]=0.+0.j#np.nan
                elif np.abs(k-ci) < kperp*H0*Dz*E*theta0/(c*(1+z)*ktorratio):    #abs as kparrallel = abs(kz)
                    Windowedimageinv[i][j][k]=0.+0.j#np.nan
    '''
    #replaces the central DC cube
    counter=0
    for i in range (ci -1,ci +2,1):
        for j in range (ci -1,ci +2,1):
            for k in range (ci -1,ci +2,1):
                Windowedimageinv[i][j][k] = centre[counter]
                counter += 1
    '''
    Windowedimageinv[ci][ci][ci]=Niquist
    return Windowedimageinv

#IMPORTANT: calcuates the rms between the two PS,  divided by 21cm PS to get unitless ratio of rms to real value
def PSrmscalc(onex,oney,twox,twoy):
    #currently (2.2.15) one is 21cm, two is windowed image
    #print oney
    #print twoy

    ''' #This is commented as it concatenated the two x axis' but now don't need to as terns out both axis are the same
    xaxis = np.zeros((1)) #just to start with, remember that first point is empty

    #this concaterates onex and oney, just not in a nice way
    for i in range(len(onex)):
        for j in range(len(twox)):
            if twox[j] == onex[i]:
                xaxis = np.append(xaxis,twox[j])
                #order won't matter as both are ordered lists without repeated numbers

    xaxis = np.delete(xaxis,0) # delete intitial 0


    print 'fraction of values compared = ', float(len(xaxis))/len(onex)
    '''

    nanbool = np.isnan(twoy)   #makes an array where all nan index's are True, not is False
    nanbool[0] = True  #manually set the first element, which is wrong to true so it#ll be ignored by rmscalc

    counter = 0 #this will save (y-y')**2 values to be averaged
    Falsecounter = 0    #this gives the number of averaged points, for rms calculation
    #we now want to go to each shared kvalue, find the relevant index, then compare the values at that index
    for i in range (len(twoy)):  #this is -1 as first point in xaxis is empty
        #indexone = np.nonzero(onex == xaxis[i])[0][0] #onex.index(xaxis[i]) #this finds the index where the shared k value is in both arrays
        #indextwo = np.nonzero(twox == xaxis[i])[0][0] #twox.index(xaxis[i])
        if nanbool[i] == False:
            counter += (1.- twoy[i]/oney[i])**2 #this is rms difference divided by 21cm value to give dimensionless
            Falsecounter += 1

    counter =float(counter)/Falsecounter  #this averages the values
    return np.sqrt(counter)  #this is rms value

##########################################printing function################################################

# This function compares the powerspectra of the image, the twenty1cmsignal and the error
# realps refers to the powerspectrum as provided by 21cmfast and can be uncommented to compare our results to this

def printpowerspectrum(oneinverse, twoinverse, threeinverse, fourinverse, psdwidth,size,dtheta, dx, z,rmsornot):


    realps = np.loadtxt('PowerSpectrumFiles/PS%s'%z, delimiter='\t')

    onek, onepowerspectrum , onedeldel= powerspectrum3D(oneinverse,psdwidth,size,dtheta,dx, z,0) # this is the size of steps in real space dx=float(box_info['dim'])/float(box_info['BoxSize'])
    print 'done imagepowerspectrum'
    twok, twopowerspectrum , twodeldel= powerspectrum3D(twoinverse,psdwidth,size,dtheta, dx, z,1)
    print 'done sigmapowerspectrum'
    threek, threepowerspectrum, threedeldel= powerspectrum3D(threeinverse,psdwidth,size,dtheta,dx, z,0)
    print 'done twenty1powerspectrum'
    #fourk, fourpowerspectrum, fourdeldel= powerspectrum3D(fourinverse,psdwidth,size,dtheta,dx, z,0)
    #print 'done twenty1powerspectrum'

    if rmsornot == 1:   #will pass in either 1 or 0 as to wether we want to calculate this
        rmsvalue = PSrmscalc(twok,twopowerspectrum,threek,threepowerspectrum)   #sends in 21cm then windowed image


    #plots the compared powerspectra
    plt.loglog(onek,onedeldel)
    plt.loglog(twok,twodeldel)
    plt.loglog(threek,threedeldel)
    #plt.loglog(fourk,fourdeldel)
    plt.loglog(realps[:,0],realps[:,1])

    plt.ylim(0.00009,300)
    plt.xlim(0.01,3)

    plt.xlabel('k (MPc$^{-1}$)')
    plt.ylabel('k$^3$ P(k)/2$\pi^2$ (mK$^2$)')
    plt.savefig('ComparingEorforPS/DELDEL_POWERSPEC_z%s.png' %z)
    plt.clf()


    plt.loglog(onek,onepowerspectrum)
    plt.loglog(twok,twopowerspectrum)
    plt.loglog(threek,threepowerspectrum)
    #plt.loglog(fourk,fourpowerspectrum)
    plt.loglog(realps[:,0],realps[:,1]/(realps[:,0]**3)) #Important: here, as with other plot, we have no 2Pi**2 factor

    plt.ylim(0.003,450000)
    plt.xlim(0.01,3)

    plt.xlabel('k (MPc$^{-1}$)')
    plt.ylabel('P(k) (mK$^2$ Mpc$^3$)')
    plt.savefig('ComparingEorforPS/POWERSPEC_z%s.png' %z)
    plt.clf()

    #saves the windowed image deldelpowerspectrum seperately, so we can choose which ones to plot together against z
    plt.loglog(threek,threedeldel)
    plt.ylim(0.003,450000)
    plt.xlim(0.01,3)
    plt.xlabel('k (MPc$^{-1}$)')
    plt.ylabel('k$^3$ P(k)/2$\pi^2$ (mK Mpc$^{-3}$)')
    plt.savefig('Powerspectrums/DELDEL_POWERSPEC_z%s.png' %z)
    plt.clf()

    if rmsornot == 1:
        return rmsvalue #this can then be saved in an array to compare rms error changing with z
    else:
        return 0

def printvaluesvsz(YVALUE,redshift,neutralfractions,labelname):

    fig = plt.figure()

    ax1 = fig.add_subplot(111) # x (z) and y axis

    ax1.plot(redshift,YVALUE)
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel('RMS Error in', labelname)

    nf_axis_ticklocations = np.array([1, 5, 9, 13, 17]) # in terms of z
    nfindexes=np.searchsorted(redshift, nf_axis_ticklocations) # finds corresponding indices

    ax2 = ax1.twiny() # further x axis corresponding to the same y axis
    ax2.set_xticks(nf_axis_ticklocations) # ticks at desires z locations.
    ax2.set_xticklabels(neutralfractions[nfindexes]) # prints nf for each z location
    ax2.set_xlabel("Un-ionized Fraction")

    plt.savefig('Statisticalvaluesvsz/%s_vsz' %labelname)   #this saves the graph using the string labelname
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
    #imagemeanpathdist = imagemeanpathdist/imagemeanpathx  #has to be done outside of return so that sum() function worksa
    twenty1meanpathx,twenty1meanpathdist = secondbubbledistributioncalculator(twenty1,size,cutoff,dl,iterations)
    #twenty1meanpathdist = twenty1meanpathdist/twenty1meanpathx  #has to be done outside of return so that sum() function worksa


    figure= plt.loglog(imagemeanpathx,imagemeanpathdist)
    plt.loglog(twenty1meanpathx,twenty1meanpathdist)

    '''
    plt.xlabel('Ionised Volume (MPc$^{3}$)')
    plt.ylabel('Probability')
    plt.xlim(2,100000)
    plt.ylim(0.000007,1)
    plt.show()
    '''
    #mean21,median21,uqmean21 = meanfreepathstatistics(twenty1meanpathx,twenty1meanpathdist)  #this sends the bubble distribution yaxis to be averaged
    #meanimage, medianimage,uqmeanimage = meanfreepathstatistics(imagemeanpathx,imagemeanpathdist)
    #returned below in the form seen above

    return meanfreepathstatistics(twenty1meanpathx,twenty1meanpathdist),meanfreepathstatistics(imagemeanpathx,imagemeanpathdist)     #this gives us back the mean for this z

def meanfreepathstatistics(xdata,ydata): #data is just the bubble size x values

    length = len(xdata)

    medianfound=False
    mean = 0
    weightedmean = 0 #walleymean is a weighted average again. i.e. a walleymean
    upperquartilemean = 0.
    upperquartilecounter = 0.
    distributioncounter=0.

    for i in range (length):
        mean += xdata[i]*ydata[i] #this is usual mean - y data is the distribution so weighting
        weightedmean += xdata[i]*xdata[i]*ydata[i]
        distributioncounter += ydata[i]
        if distributioncounter > 0.75: # start finding mean once we are in the upper quartile
            upperquartilemean += xdata[i]*ydata[i]
            upperquartilecounter += ydata[i]
        if distributioncounter >= 0.5 and medianfound == False:
        # define median as first time distribution counter is above 0.5 // not quite correct but should be fine.
            median = xdata[i]
            medianfound = True
            print 'at distribution', distributioncounter, 'our medianish is ', median

    upperquartilemean = upperquartilemean/upperquartilecounter

    return mean, median, upperquartilemean, weightedmean/mean