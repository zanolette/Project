import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo
import glob

###############################Taken from Main#################################

psdwidth = 1

#Define size of view and resolution
# defines our universe
CosmoUnits=Cosmo.CosmoUnits()

#These arrays for plotting variable/statistic changes with z/neutral fraction
PSrmsarray = np.zeros((2, 28))   #this is for z=7.5 to z=18 runs in half steps, plus some .25 steps before 10.5
averagetemp = np.zeros((3, 28))  #this is the average temp: 21cm,image,Windowed
rmserrorintemp = np.zeros((2, 28))   #this gives rms difference between two images at every point
PearsonRarray = np.zeros((3, 28)) #this saves z,image vs 21cm Pearson,windowed vs 21cm Parson
MeanValues = np.zeros((12,28))  #will save all of our different mean values, and for 21cm,image and windowed

redshift = np.zeros(28)
neutralfractions = np.zeros(28)

counter=0   #this is just to help allocate saved values into the arrays

##########Taking each z box seperately#############

path = "LARGEBOXES2/*"
for fname in glob.glob(path):

    box_info = boximport.parse_filename(fname)
    box = boximport.readbox(fname)

    twenty1 = box.box_data  #so we have it seperate - this is 3D!

    dl=float(box_info['dim'])/float(box_info['BoxSize'])
    z = box_info['z']
    theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
    size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files
    dtheta = float(theta/size)

    redshift[counter] = z   #saved for later to put labels on axes
    neutralfractions[counter] = box_info['nf']  #saves neutral fraction of this z

    print 'Now starting z =',z

    #############Load in Files for this z########################################################

    image3Dinverse = np.load('Experiment/image3Dinv_z%s.npy' %z)
    sigma3Dinverse = np.load('Experiment/sigma3Dinv_z%s.npy' %z)
    Windowedimageinverse = np.load('Experiment/windowedinv_z%s.npy' %z)

    func.kperpvskparrgraph(image3Dinverse,psdwidth,size,dl,z,'image',  1e7, 2e11 )
    func.kperpvskparrgraph(Windowedimageinverse,psdwidth,size,dl,z,'Windowed',1e7, 2e11)
    func.kperpvskparrgraph(sigma3Dinverse,psdwidth,size,dl,z,'Sigma',10, 1e4)

    ###############Calculating fft's################################

    twenty1inverse = np.fft.fftn(twenty1)   #gives 3D FFT of 21cm box!
    twenty1inverse = np.fft.fftshift(twenty1inverse)

    #How we get our image from the fourier transform
    image3D = np.fft.ifftn(image3Dinverse)
    image3D = np.abs(image3D)

    Windowedimage = np.fft.ifftn(Windowedimageinverse)
    Windowedimage = np.abs(Windowedimage)   #abs or real?


    ##################Calculating and saving statistical values############################

    #Calcules PearsonR Value for non-windowed then windowed with z
    PearsonRarray[0][counter] = z
    PearsonRarray[1][counter] = func.PearsonR(twenty1,image3D,size)
    PearsonRarray[2][counter] = func.PearsonR(twenty1,Windowedimage,size)


    #Temperature statistics
    averagetemp[0][counter] = np.average(twenty1)  #this saves the average 21cm temperature
    averagetemp[1][counter] = np.average(image3D)  #this saves the average image temperature
    averagetemp[2][counter] = np.average(Windowedimage)   #this saves the average windowed image temperature
    rmserrorintemp[0][counter] = func.rmscalc(twenty1,image3D,size)
    rmserrorintemp[1][counter] = func.rmscalc(twenty1,Windowedimage,size)
    #print 'z', z,'nf',box_info['nf'],  'temp',average21cmtemp[counter],averagewindowedimagetemp[counter], 'rms', rmserrorintemp[counter]


    #func.visualisereionizationslicebyslice(Windowedimage, twenty1, size, z, theta, True,'Windowed')
    #func.visualisereionizationslicebyslice(image3D, twenty1, size, z, theta, True,'Non-Windowed') #don't think we'll need this


    #This function compared the phases of the real and imaginary
    #func.phasecomparison(twenty1inverse, Windowedimageinverse, size, z)


    #!!printpowerspectrum is for comparing the windowed,non-windowed and 21cm powerspectrums on one graph, but also saves deldelPS's seperately for comparison!!
    #The first of these is the non-windowed PS rms, the second is the windowed PS rms
    PSrmsarray[0][counter],PSrmsarray[1][counter]=func.printpowerspectrum(image3Dinverse, twenty1inverse, Windowedimageinverse, sigma3Dinverse, psdwidth,size,dtheta,dl, z,1)    ##,saves rms error between windowed and 21cm.
    #print 'z', z, 'PSrms', PSrmsarray[counter]


    #The cutoff refers to the fraction of the average temperature at which the code defines a point to be ionised
    cutoff = 0.65
    iterations = 100000


    #These functions print the different size distribution analysis methods which we have worked on
    #Saves the text files for the distributions and their statistics plus prints distributions together
    #these are: mean21,median21,uqmean21,weightedmean21,meanimage, medianimage,uqmeanimage,weightedmeanimage,meanwindowed,medianwindowed,uqmeanwindowed,weightedmeanwindowed
    MeanValues[:][counter] = func.printmeanfreepathdist(image3D,Windowedimage,twenty1, size, dl, cutoff, iterations,z)
    func.printbubblesizedist(image3D,Windowedimage, twenty1, size, dl, cutoff,z)

    ''' - Temperature distribution code and how to plot it
    distributionplot=plt.figure()
    plt.plot(func.temperaturedistribution(Windowedimage, size))
    plt.plot(func.temperaturedistribution(twenty1, size))
    plt.savefig('TempDistribution/TempHisto%s.png'%z)
    plt.clf()
    '''
    counter += 1    #this increases the counter so the next values are stored in the next element slots of the arrays


#print all of the arrays against z/neutral fraction
#printing z and neutral fraction on two x axes
#ASSUMES - 3 numpy arrays - z, neutrlfractions and YVALUE
func.printvaluesvsz(PSrmsarray,redshift,neutralfractions,'RMS in PowerSpectrum (mK$^2$ Mpc$^3$)')   #error in powerspectrum vs z
func.printvaluesvsz(PearsonRarray,redshift,neutralfractions,'Pearson R Value in Temperature (mK)')  #error in temp vs z
func.printvaluesvsz(rmserrorintemp,redshift,neutralfractions,'RMS in Temperature (mK)')  #error in temp vs z


###############This shows both 21cm and windowed image average temperatures varying with time
func.averagetempvsz(averagetemp,neutralfractions,redshift)


############Saving all arrays so they can be outputted#########################
np.savetxt('TextFiles/PSrmsOutputFile',PSrmsarray, delimiter='\t')
np.savetxt('TextFiles/AverageTempOutputFile',averagetemp, delimiter='\t')
np.savetxt('TextFiles/RMSErrorsinTempOutputFile',rmserrorintemp, delimiter='\t')
np.savetxt('TextFiles/PearsonROutputFile',PearsonRarray, delimiter='\t')



