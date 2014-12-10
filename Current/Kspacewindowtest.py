import numpy as np
import matplotlib.pyplot as plt

testarray = np.zeros((2,100))
testarray[0] = np.arange(0,0.1,0.001)

counter = 0.001

for i in range(100):
    if testarray[0][i]<0.055:
        testarray[1][i]=0.5
    else:
        testarray[1][i]=np.sqrt(testarray[0][i]-0.055)+0.5


plt.plot(testarray[0],testarray[1])
plt.show()


plt.loglog(testarray[0],testarray[1])
plt.show()


