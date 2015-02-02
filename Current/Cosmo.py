
# A class that lays the ground works for getting the right cosmological units.

import numpy as np

class CosmoUnits:   #all except thetaboxsize are ratio's or in MPc, but for angles need to remember rad - degrees
    #Fitted to same values as 21cmfast Cosmology.h
    def __init__(self,
                 H0=70.0,           # Hubble Constant
                 omega_m=0.28,     # Mass Density
                 omega_l=0.72):    # Dark energy Density

        self.H0 = float(H0)
        self.omega_m = float(omega_m)
        self.omega_l = float(omega_l)
        self.omega_k = float(1 - omega_m - omega_l)

        self.C=2.99792458e5 # speed of light km per second because hubble was weird
        self.CMpcperGigyear=3.06601394e2


    #Using the Wikipedia calculation for all for these

    def inv_Efunc(self, z):
        #This is the Hubble Parameter
        return 1/np.sqrt(self.omega_m*((1+z)**3) + self.omega_k*((1+z)**2) + self.omega_l)

    def Dhubble(self):       # Hubble Distance

        return self.C/self.H0

    def Dcomovingrad(self, z):    # Comoving Distance - radial

        from scipy.integrate import quad
        return self.Dhubble() * quad(self.inv_Efunc, 0, z)[0]

    def Dcomovingtrans(self, z): #Transverse Comoving Distance between two objects (i.e. cancels expansion of the universe)

        if self.omega_k>0:
            return self.Dhubble() * np.sinh( np.sqrt(self.omega_k) * self.Dcomovingrad(z) / self.Dhubble() )/np.sqrt(self.omega_k)
        if self.omega_k==0:
            return self.Dcomovingrad(z)
        if self.omega_k<0:
            return self.Dhubble() * np.sin( np.sqrt(np.mod(self.omega_k)) * self.Dcomovingrad(z) / self.Dhubble() )/np.sqrt(np.mod(self.omega_k))

    def Dangular(self,z): # Distance given an observed angle theta - get this explained better.

        return self.Dcomovingtrans(z)/(1+z)

    def Dlum(self,z): # luminosity distance

        return self.Dcomovingtrans(z)*(1+z)

    def mod_inv_Efunc(self, z): #needed for light travel calc below
            #This is the Hubble Parameter
            return 1./(1.+z) * 1/np.sqrt(self.omega_m*((1+z)**3) + self.omega_k*((1+z)**2) + self.omega_l)

    def Dlighttravel(self, z):

        from scipy.integrate import quad
        return self.Dhubble() * quad( self.mod_inv_Efunc, 0, z)[0]

    def Tlighttravel(self, z): # returns result in gigayears
        return self.Dlighttravel(z)/self.CMpcperGigyear

    def thetaboxsize(self,z,MPc):
        return 360.*MPc/(self.Dcomovingtrans(z)*2.*np.pi)

#x=CosmoUnits(70,0.3,0.7)

#z=3.

#print 'given a redshift of',z

#print 'light travel time:', x.Tlighttravel(z),'Gigayears'
#print 'comoving radial distance:',x.Dcomovingrad(z) , 'Mpc'
#print 'angular distance:', x.Dangular(z), 'Mpc'
#print 'luminosity distance:', x.Dlum(z), 'Mpc'
