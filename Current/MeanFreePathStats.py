import numpy as np
import Functions as func
import matplotlib.pyplot as plt

size = 100


radius = 5

# specify circle parameters: centre ij and radius
image = np.ones((size,size,size))
image=image*2.

ci,cj=int(size/4),(size/2) #centre point

cr=radius  #arbitary for now

# Create index arrays to z
I,J,K=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]), np.arange(image.shape[2]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
image[np.where(np.sqrt((I-ci)**2+2*(J-cj)**2+(K-ci)**2)<cr)]=0.

#doing a second sphere
cr = 15
ci,cj=int(3*size/4),(size/2) #centre point

# Create index arrays to z
I,J,K=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]), np.arange(image.shape[2]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
image[np.where(np.sqrt((I-ci)**2+(J-cj)**2+(K-ci)**2)<cr)]=0.

#doing a second sphere
cr = 3
ci,cj,ck=int(3*size/4),int(3*size/4),int(5*size/7) #centre point

# Create index arrays to z
I,J,K=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]), np.arange(image.shape[2]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
image[np.where(np.sqrt((I-ci)**2+(J-cj)**2+(K-ck)**2)<cr)]=0.

#doing a second sphere
cr = 3
ci,cj,ck=int(1*size/4),int(1*size/4),int(2*size/7) #centre point

# Create index arrays to z
I,J,K=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]), np.arange(image.shape[2]))

# np.sqrt((I-ci)**2+(J-cj)**2) calculates distance of all points to centre
# Assign value of 1 to those points where distance to centre is less that cr:
image[np.where(np.sqrt((I-ci)**2+(J-cj)**2+(K-ck)**2)<cr)]=0.





print 'go'

x,y=func.secondbubbledistributioncalculator(image,size, 0.5,1, 100000)

plt.loglog(x,y)
plt.ylim(0.0000000001,1)
plt.show()