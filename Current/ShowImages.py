#This code is made to print as many images next to each other as possible. - Think it would be good for comparing what the EOR window does

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

    redshift[counter] = z   #saved for later to put labels on axes
    neutralfractions[counter] = box_info['nf']  #saves neutral fraction of this z

    image3Dinverse = np.load('Experiment2/image3Dinv_z%s.npy' %z)
    Windowedimageinverse = np.load('Experiment2/windowedinv_z%s.npy' %z)


    #How we get our image from the fourier transform
    image3D = np.fft.ifftn(image3Dinverse)
    image3D = np.abs(image3D)

    Windowedimage = np.fft.ifftn(Windowedimageinverse)
    Windowedimage = np.abs(Windowedimage)   #abs or real?

    for t in range(size):
        fig = plt.figure()
        a1=fig.add_subplot(1,3,1)
        imgplot = plt.imshow(twenty1[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
        a1.set_title('21cmfast', size = 8)
        plt.xlabel('X axis in $^\circ$s', size = 8)
        plt.ylabel('Y axis in $^\circ$s', size = 8)
        a2=fig.add_subplot(1,3,2)
        plt.setp( a2.get_yticklabels(), visible=False)
        imgplot = plt.imshow(image3D[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
        plt.xlabel('X axis in $^\circ$s', size = 8)
        a2.set_title('SKA Image, Full k Space',size = 8)
        a3=fig.add_subplot(1,3,3)
        plt.setp( a3.get_yticklabels(), visible=False)
        imgplot = plt.imshow(Windowedimage[t],extent=(-theta/2,theta/2,-theta/2,theta/2),interpolation='nearest',cmap='jet',vmin=0,vmax=70)
        plt.xlabel('X axis in $^\circ$s', size = 8)
        a3.set_title('SKA Image with EOR Window Avoidance', size = 8)

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.3, 0.03, 0.4])
        cbar=plt.colorbar(imgplot, orientation='vertical',cax=cbar_ax)# shrink=0.6)

        cbar.set_label('Temperature (mK)', size = 8)

        plt.savefig('Compare3Images/z%simage%03i.png'%(z,t))
        plt.close(fig)