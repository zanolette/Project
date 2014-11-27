

import numpy as np
import matplotlib.pyplot as plt
import Functions as func

size=100
cube=np.random.randint(2, size=(size, size,size))

x, dist= func.bubblesizedistribution(cube, size)

plt.loglog(x, dist)
plt.show()

