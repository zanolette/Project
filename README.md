Project
=======
import numpy as np
import matplotlib as mpl

e = numpy.zeros((100,100))
c = 50
ra = 3
j = -r
i = -floor(sqrt(ra**2 - j**2))


while j <= r:
    while i <= ceil(sqrt(ra**2 - j**2)):
        e[i+x0][j+y0] = 1
        i += 1
    j += 1


        
  '''
a = array([[1,2],[3,4]])
b = fft.fft2(a)

c = array ([[0,0,0],[0,1,0],[0,0,0]])
d = fft.fft2(c)

for j in range(2*r):    
    for i in range(1 + 2*ceil(sqrt(r**2 - (j-r)**2))):
        e[i + x0 - ceil(sqrt(r**2 - (j-r)**2))][j+y0-r] = 1
  
[[e[i][j] for j in range(len(e)) ]
                        for i in range(len(e))]

if  abs(s - f) < floor(sqrt(f^2 - j^2))
  '''


s = array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])

