import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import ticker
import Functions as func
import boxio as boximport
import Cosmo
import glob

###############################Taken from Main#################################

psdwidth = 2

#Define size of view and resolution
# defines our universe
CosmoUnits=Cosmo.CosmoUnits()

#These arrays for plotting variable/statistic changes with z/neutral fraction
PSrmsarray = np.zeros(23)   #this is for z=7 to z=18 runs in half steps
average21cmtemp = np.zeros(23)  #this is the average temp
averagewindowedimagetemp = np.zeros(23) #compared to the above
rmserrorintemp = np.zeros(23)   #this gives rms difference between two images at every point

redshift = np.zeros(23)
neutralfractions = np.zeros(23)

counter=0   #this is just to help allocate saved values into the arrays

##########Taking each z box seperately#############

path = "LARGEBOXES/*"
for fname in glob.glob(path):

    box_info = boximport.parse_filename(fname)
    box = boximport.readbox(fname)

    twenty1 = box.box_data  #so we have it seperate - this is 3D!

    dl=float(box_info['dim'])/float(box_info['BoxSize'])
    z = box_info['z']
    theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
    size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files
    dtheta = float(theta/size)

    print z

    redshift[counter] = z   #saved for later to put labels on axes
    neutralfractions[counter] = box_info['nf']  #saves neutral fraction of this z

    #############Load in Files for this z########################################################

    image3Dinverse = np.load('Experiment/image3Dinv_z%s.npy' %z)
    sigma3Dinverse = np.load('Experiment/sigma3Dinv_z%s.npy' %z)
    Windowedimageinverse = np.load('Experiment/windowedinv_z%s.npy' %z)

    ###############Calculating fft's################################

    twenty1inverse = np.fft.fftn(twenty1)   #gives 3D FFT of 21cm box!
    twenty1inverse = np.fft.fftshift(twenty1inverse)

    #How we get our image from the fourier transform
    #image3D = np.fft.ifftn(image3Dinverse)
    #image3D = np.abs(image3D)

    #Windowedimage = np.fft.ifftn(Windowedimageinverse)
    #Windowedimage = np.abs(Windowedimage)   #abs or real?


    ##################Calculating and saving statistical values############################

    #average21cmtemp[counter] = np.average(twenty1)  #this saves the average 21cm temperature
    #averagewindowedimagetemp[counter]=np.average(Windowedimage)   #this saves the average image temperature

    #rmserrorintemp[counter] = func.rmscalc(twenty1,Windowedimage,size)

    #func.visualisereionizationslicebyslice(Windowedimage,twenty1, size, z, theta)

    #This function compared the phases of the real and imaginary
    #func.phasecomparison(twenty1inverse, Windowedimageinverse, size)

    #!!printpowerspectrum is for comparing the windowed,non-windowed and 21cm powerspectrums on one graph, but also saves deldelPS's seperately for comparison!!
    PSrmsarray[counter]=func.printpowerspectrum(image3Dinverse, twenty1inverse, Windowedimageinverse, sigma3Dinverse, psdwidth,size,dtheta,dl, z,1)    ##,saves rms error between windowed and 21cm.

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

    counter += 1    #this increases the counter so the next values are stored in the next element slots of the arrays

#print all of the arrays against z/neutral fraction
#printing z and neutral fraction on two x axes
#ASSUMES - 3 numpy arrays - z, neutrlfractions and YVALUE

func.printvaluesvsz(PSrmsarray,redshift,neutralfractions,'rmserrorforPS')   #error in powerspectrum vs z
'''
func.printvaluesvsz(rmserrorintemp,redshift,neutralfractions,'rmserrorinTemp')  #error in temp vs z


###############This shows both 21cm and windowed image average temperatures varying with time
fig = plt.figure()

ax1 = fig.add_subplot(111) # x (z) and y axis

ax1.plot(redshift,averagewindowedimagetemp,redshift,average21cmtemp)
ax1.set_xlabel("Redshift")
ax1.set_ylabel("Y QUANTITY")

nf_axis_ticklocations = np.array([7, 9, 11, 13, 15, 17]) # in terms of z
nfindexes=np.searchsorted(redshift, nf_axis_ticklocations) # finds corresponding indices

ax2 = ax1.twiny() # further x axis corresponding to the same y axis
ax2.set_xticks(nf_axis_ticklocations) # ticks at desires z locations.
ax2.set_xticklabels(neutralfractions[nfindexes]) # prints nf for each z location
ax2.set_xlabel("Un-ionized Fraction")

plt.savefig('Statisticalvaluesvsz/averagetemperaturecomparisson')   #this saves the graph using the string labelname
plt.clf()

'''