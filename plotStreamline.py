import numpy as np
import matplotlib.pyplot as plt

data = np.fromfile('lidCavityResult.dat', dtype='double')
x, y = data[::4], data[1::4]
#y, x = np.mgrid[0.005:0.995:100j, 0.005:0.995:100j]
u, v = data[2::4], data[3::4]

u = np.reshape(u, (100, 100))
v = np.reshape(v, (100, 100))
x = np.reshape(x, (100, 100))
y = np.reshape(y, (100, 100))

#speed = np.sqrt(u * u + v * v)

plt.figure()
plt.streamplot(x, y, u, v, density=(3, 3), color='k')
plt.show()
