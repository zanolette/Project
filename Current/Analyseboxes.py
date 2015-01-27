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

fname = 'delta_T_v2_no_halos_nf0.926446_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.68_100_200Mpc'#'delta_T_v2_no_halos_nf0.932181_z14.00_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb30.80_200_400Mpc'#
box_info = boximport.parse_filename(fname)

#Define size of view and resolution
# defines our universe
CosmoUnits=Cosmo.CosmoUnits()
z = box_info['z']
theta = CosmoUnits.thetaboxsize(z,box_info['BoxSize'])
size = box_info['dim'] # box size - maybe write a code to get this out of the title of the 21cmfast files
dtheta = float(theta/size)

path = "LARGEBOXES/*"
for fname in glob.glob(path):
    box=boximport.readbox(fname)
    twenty1 = box.box_data  #so we have it seperate - this is 3D!

    dx=float(box_info['dim'])/float(box_info['BoxSize'])

    ###############################################################################
    #Load in Files for this z:


    image3Dinverse = np.load('experiment/image3Dinv_z%s.npy' %z)
    sigma3Dinverse = np.load('experiment/sigma3Dinv_z%s.npy' %z)
    Windowedimageinverse = np.load('experiment/windowedinv_z%s.npy' %z)

    ###############################################


    #How we get our image from the fourier transform
    image3D = np.fft.ifftn(image3Dinverse)
    image3D = np.abs(image3D)

    Windowedimage = np.fft.ifftn(Windowedimageinverse)
    Windowedimage = np.abs(Windowedimage)   #abs or real?

    #func.visualisereionizationslicebyslice(Windowedimage,twenty1, size, z, theta)

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

