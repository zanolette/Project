

import numpy as np
import matplotlib.pyplot as plt
import Functions as func

size=20
cube=np.random.randint(2, size=(size, size,size))

dist= func.bubblesizedistribution(cube, size)

plt.loglog(dist)
plt.show()

