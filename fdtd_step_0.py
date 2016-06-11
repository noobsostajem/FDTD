import numpy as np
import math as mt
import matplotlib.pyplot as plt

 size=100;maxTime=347
 -#for size=100 use time 347
 -#size=1000;maxTime=1048
ez=np.zeros((size))
hy=np.zeros((size))
ez_snap=np.zeros((maxTime))
hy_snap=np.zeros((maxTime))
imp0=377.0
snap_moment=size/2-1

source_width = 5
delay = 10*source_width

i=0

for qtime in np.arange(0,maxTime+1,1):
    print(qtime)
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2));

    for mm in range(0,size-1):
        hy[mm]=hy[mm]+(ez[mm+1]-ez[mm])/imp0
    for mm in range(1,size):
        ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0

    ez_snap[i-1]=ez[size/2-1]
    i=i+1

fig = plt.figure()

plt.subplot(2, 1, 1)
plt.plot(range(0,size,1),ez,'-')
plt.plot(range(0,size,1),hy*imp0,'-')
plt.ylabel('Ez-field, V/m')
plt.xlabel('Spatial')

#Ez (Spatial)
#plt.subplot(2, 1, 2)
#plt.plot(range(0,size,1),ez,'-')
#plt.ylabel('Ez-field, V/m')
#plt.xlabel('Spatial')
#plt.show()

#Ez (Time)
plt.subplot(2, 1, 2)
plt.plot(range(0,maxTime,1),ez_snap,'-')
plt.ylabel('Ez-field, V/m')
plt.xlabel('Time')
plt.show()
