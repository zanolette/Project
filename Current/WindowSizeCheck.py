

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

    #DEPENDS ON Z! FIND THIS
    dl = box_info['BoxSize']/size   #so this is how many Mpc per index of box
    psdwidth = 1    #can change this!

    #DEFINE EXPERIMENT PARAMETERS
    H = 0.
    tint = 300.      #interval in seconds
    dH = tint*(2.*np.pi) / (60.*60.* 24.)     #this is 2pi/time - converts from time interval (seconds) to angle interval
    totalintegrationtime = 24    #total time in hours
    daysofscanning = 1 # keep as 1 unless you want more than one day.
    timestepsneeded= int(totalintegrationtime * 60 * 60 / tint) # unitlessmeasurement of number of steps needed
    delta = 90./180. * np.pi    #declination angle
    scaling = 1./(size*dtheta) #the length of one interval in inverse space

    # need to make sure this is correct according to pritchard - noise on complex thing
    tsyst = 50000 + 60000*((1+z)/4.73)**2.55  #(mK) this is from "Probing . . . with the SKA" MG Santos
    B = 1420.41e6*dl/((1+z)*CosmoUnits.Dcomovingrad(z))        #Bandwidth (Hz) - this is given by frequency inteval. DeltaF = f1-f2, seperated by deltaz = z* deltaL/Dcomovingrad


    Windowedimageinverse=np.ones((size,size,size), complex)
    Windowedimageinverse=func.EORWINDOW(Windowedimageinverse,size,dl, z,B) # create new variable early, function was changing original variable
    #print Windowedimageinverse
    #print 'made windowed image'

    func.kperpvskparrgraph(Windowedimageinverse,psdwidth,size,dl)

    print 'z', z, 'inversespacecovered', np.average(Windowedimageinverse)