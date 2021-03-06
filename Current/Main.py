

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo
import glob

# defines our universe
CosmoUnits=Cosmo.CosmoUnits()

#remember generic print function - func.printgraph (image, xrange, yrange, xlabel, ylabel,scalemin,scalemax)

#getting 21cm box information
#fname = 'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'#'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'  #
path = "LARGEBOXES/*"
for fname in glob.glob(path):

    box_info = boximport.parse_filename(fname)
    box = boximport.readbox(fname)

    #Define size of view and resolution
    z = box_info['z']
    theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
    print theta

    #so box_info['BoxSize'] is the Mpc of box length, whereas box_info['dim'] is index size of box length
    size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files

    ci = int(size/2)

    dtheta = float(theta/size)
    print dtheta

    arrayname = 'MWA128.txt'

    #Define Wavelength - find this out from z!!
    lda=0.21106*(1+z)   #in meters

    if (arrayname=='LOFAR66CoreHB.txt'):
        Aeff = min(lda**2/3.,1.5625)
        N = 16*24*66 #number of tiles in a station * number of stations in an array
    else:   #hence assuming its MWA
        Aeff = min(lda**2/3.,14.5/16)
        N = 16*128

    Atotal = Aeff*N  #instrument efficiency: here Aeff is effective area per dipole, N is no. dipoles per station


    #DEPENDS ON Z! FIND THIS
    dl = box_info['BoxSize']/size   #so this is how many Mpc per index of box
    psdwidth = 1    #can change this!


    #DEFINE EXPERIMENT PARAMETERS
    H = 0.
    tint = 300.      #interval in seconds
    dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
    totalintegrationtime = 24    #total time in hours

    daysofscanning = 1# keep as 1 unless you want more than one day.
    timestepsneeded= 1#int(totalintegrationtime * 60 * 60 / tint) # unitlessmeasurement of number of steps needed

    delta = 90./180. * np.pi    #declination angle
    scaling = 1./(size*dtheta) #the length of one interval in inverse space

    #CODE THAT CAN BE USED TO CHECK OUR SPATIAL FREQUENCY IS IN THE RIGHT UNITS
    #spatialfreq=np.fft.fftfreq(size, dtheta)
    #print spatialfreq[1] - spatialfreq[0]
    #print scaling


    # Now we import array positions
    (dx,dy,dz)=func.importarray(arrayname,lda) #this assumes lda doesn't change too much over slices!!! 'MWAcoordinate.txt''vla.a.cfg' 'MWAhalved.txt''MWA128.txt'

    #this calculated the longest baseline
    baselinelengths=np.sqrt(np.multiply(dx,dx)+np.multiply(dy,dy))
    Dmax = np.max(baselinelengths)
    print Dmax
    baselinelengths = None

    eps = Dmax**2/Atotal
    print eps

    #Apply rotation matrix onto baseline vector and maps onto fourier plane.
    UVcount,uvcoveragepercentage = func.rotationmatrix(dx, dy, dz, scaling, H, dH, timestepsneeded, delta, size)
    print 'UV coverage percentage at z = %s is' %z, uvcoveragepercentage

    # need to make sure this is correct according to pritchard - noise on complex thing
    tsyst = 50000 + 60000*((1+z)/4.73)**2.55  #(mK) this is from "Probing . . . with the SKA" MG Santos
    B = 1420.41e6*dl/((1+z)*CosmoUnits.Dcomovingrad(z))        # is this Bband or Bchannel - think its already Bchannel
    #Bandwidth (Hz) - this is given by frequency inteval. DeltaF = f1-f2, seperated by deltaz = z* deltaL/Dcomovingrad


    #imports 21cm box and takes single z slice
    #zhat,twenty1 = func.twentyonecmmatrix(fname,theta/2)    #z will be used later to compute rms of image and measured sky
    box=boximport.readbox(fname)
    twenty1 = box.box_data  #so we have it seperate - this is 3D!


    twenty1inverse = np.fft.fftn(twenty1)   #gives 3D FFT of 21cm box!
    twenty1inverse = np.fft.fftshift(twenty1inverse)

    image3Dinverse = np.zeros((size,size,size),'complex')   #so we can save into a 3D box in fourier space
    sigma3Dinverse = np.zeros((size,size,size))

    ##########################This is where we introduce a slice################################

    #using UV count - this now merges the UVcoverage and the Image

    centre = np.zeros((27),dtype=complex)
    counter = 0

    #this is saving the centre cube, so it is safe from masking - SAVES NIQUIST FREQUENCIES

    for i in range (ci,ci +1,1):
            for j in range (ci,ci +1,1):
                for k in range (ci,ci +1,1):
                    centre[counter] = twenty1inverse[i][j][k]
                    #print Windowedimageinv[i][j][k]
                    print counter

                    counter += 1

    for slice in range(size):   #iterates over all slices

        print 'z', z, ' Slice: ', slice

        zhat = twenty1inverse[slice]    #takes a 2D slice of 21cm in fourier space

        ####################################MEASURE##################################################

        # for loop that goes through the fourier space matrix and adds noise according to the number of times the signal gets sampled

        for i in range (size):
            for j in range (size):
                if UVcount[i][j] != 0:
                    sigma3Dinverse[slice][i][j] = (tsyst*eps)/(np.sqrt(UVcount[i][j]*tint*daysofscanning*B))       #saved seperately to calculate power spectrum seperately, error eqn according to NRAO course + Pritchard
                    real=np.random.normal(np.real(zhat[i][j]), sigma3Dinverse[slice][i][j]/np.sqrt(2), 1)     #sqrt(2) here as real and imag components share it
                    imaginary = np.random.normal(np.imag(zhat[i][j]), sigma3Dinverse[slice][i][j]/np.sqrt(2), 1)
                    #####################!!!!!!!!this saves as [k][i][j]!!!!!!#######################
                    image3Dinverse[slice][i][j]=real[0] + imaginary[0]*1j

    #replaces the central DC cube
    counter=0
    for i in range (ci,ci +1,1):
        for j in range (ci,ci +1,1):
            for k in range (ci,ci +1,1):
                print image3Dinverse[i][j][k]
                image3Dinverse[i][j][k] = centre[counter]
                counter += 1

    #Could have psf stuff here

    #How we get our image from the fourier transform
    #image3D = np.fft.ifftn(image3Dinverse)
    #image3D = np.abs(image3D)

    #now we have the actual image in k space but we want to add a k space window onto it.
    Windowedimageinverse=np.copy(image3Dinverse) # create new variable early, function was changing original variable
    Windowedimageinverse=func.EORWINDOW(Windowedimageinverse, size, dl,z,B)

    #Windowedimage = np.fft.ifftn(Windowedimageinverse)
    #Windowedimage = np.abs(Windowedimage)   #abs or real?

    #func.visualisereionizationslicebyslice(Windowedimage,twenty1, size, z, theta,True, 'Windowed')

    #This function compared the phases of the real and imaginary
    #func.phasecomparison(twenty1inverse, Windowedimageinverse, size)

    #func.printpowerspectrum(image3Dinverse, twenty1inverse, Windowedimageinverse, sigma3Dinverse, psdwidth,size,dtheta,dl, z,1)    ##will pass in either 1 or 0 as to whether we want to calculate this

    #The cutoff refers to the fraction of the average temperature at which the code defines a point to be ionised
    #cutoff = 0.65
    #iterations = 10000

    #These functions print the different size distribution analysis methods which we have worked on
    #THIS FUNCTION SHOWS THE PLOT RATHER THAN SAVING - NEED TO CHANGE BEFORE AUTOMATION
    #func.printmeanfreepathdist(image3D, twenty1, size, dl, cutoff, iterations)
    #func.printbubblesizedist(image3D, twenty1, size, dl, cutoff)

    #func.printmeanfreepathdist(Windowedimage, twenty1, size, dl, cutoff, iterations)
    #func.printbubblesizedist(Windowedimage, twenty1, size, dl, cutoff)

    #This function compares two different 21cmboxes (maybe change it so you can insert
    #func.visualisereionizationslicebyslice(image3D,twenty1, size, z, theta)

    # This function compares the powerspectra of the image, the twenty1cmsignal and the error
    # func.printpowerspectrum(image3Dinverse, sigma3Dinverse, twenty1inverse, psdwidth,size,dtheta,dl, z)

    #This function shows the kperpendicular vs kparallel graph
    #func.kperpvskparrgraph(Windowedimageinverse,psdwidth,size,dl,z,'Windowed')
    #func.kperpvskparrgraph(twenty1inverse,psdwidth,size,dl,z,'21cmimage')
    #func.kperpvskparrgraph(sigma3Dinverse,psdwidth,size,dl,z,'Sigma')

    #IF YOU WANT TO SAVE BOXES FOR LATER ANALYSIS - USE THESE
    #np.save('Experiment/image3D_z%s'%(z),image3D)
    #del image3D
    np.save('Experiment_MWA1Shot/image3Dinv_z%s'%(z),image3Dinverse)
    del image3Dinverse
    np.save('Experiment_MWA1Shot/sigma3Dinv_z%s'%(z),sigma3Dinverse)
    del sigma3Dinverse
    np.save('Experiment_MWA1Shot/windowedinv_z%s'%(z),Windowedimageinverse)
    del Windowedimageinverse
    #np.save('Experiment/finalrealimage_z%s'%(z),Windowedimage)
    #del Windowedimage

    #np.save('image3Darraydim%s,%sMpc,z%s,time%s'%(size,box_info['BoxSize'],z,timestepsneeded),image3D)
    #np.save('sigma3Darraydim%s,%sMpc,z%s,time&s'%(size,box_info['BoxSize'],z,timestepsneeded),sigma3D)
