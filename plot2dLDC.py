#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import math as math


data = np.fromfile('lidCavityResult.dat', dtype='double')
x, y = data[::3], data[1::3]
umag = data[2::3]

size = math.sqrt(len(umag))
umag = np.reshape(umag, (size, size))
x = np.reshape(x, (size, size))
y = np.reshape(y, (size, size))

plt.figure()
levels = np.linspace(0, 1, 11)
im = plt.contourf(x, y, umag, levels=levels)
#plt.contour(x, y, umag)
t = plt.title("Contour of Velocity Magnitude for Lid Driven Cavity at Re = 1000 MRT")
t.set_y(1.015)
plt.colorbar(im)
#plt.show()
plt.savefig('ldcFilledContourMRT.eps')
plt.close()

#Plot vertical velocity
d1 = np.loadtxt('verticalGhia.dau')
y1 = (d1[:,0]-1)/128
u1 = d1[:,1]
d2 = np.loadtxt('ldcVerticalCenterline.dat')
y2 = (d2[:,0])
u2 = d2[:,1]
plt.figure()
plt.plot(y1, u1, 'x', label='Ghia et. al.');
plt.plot(y2, u2, label='LBM');
plt.legend()
plt.savefig('ldcVerticalComparisonMRT.eps')
plt.close()

#Plot horizontal velocity
d1 = np.loadtxt('horizGhia.dau')
y1 = (d1[:,0]-1)/128
u1 = d1[:,1]
d2 = np.loadtxt('ldcHorizontalCenterline.dat')
y2 = (d2[:,0])
u2 = d2[:,1]
plt.figure()
plt.plot(y1, u1, 'x', label='Ghia et. al.');
plt.plot(y2, u2, label='LBM');
plt.savefig('ldcHorizComparisonMRT.eps')
plt.close()
