Project
=======
import numpy as np
import matplotlib.pyplot as plt

z = np.zeros((400,400))

# specify circle parameters: centre ij and radius
ci,cj=200,200
cr=150

# Create index arrays to z
I,J=np.meshgrid(np.arange(z.shape[0]),np.arange(z.shape[1]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
z[np.where(np.sqrt((I-ci)**2+(J-cj)**2)<cr)]=1

#Fast Fourier Transform of z
zhat=np.fft.fft2(z)
zhat=abs(zhat)

print zhat

# show result in a simple plot
#fig=plt.figure()
#ax=fig.add_subplot(111)
#ax.pcolormesh(z)
#ax.set_aspect('equal')
#plt.show()

#fig2=plt.figure()
#ax=fig2.add_subplot(111)


im = plt.imshow(zhat, cmap='hot')
plt.colorbar(im, orientation='horizontal')
plt.show()
