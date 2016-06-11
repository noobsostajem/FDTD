import numpy as np
import math as mt
import matplotlib.pyplot as plt
plt.ion()

size=1000;maxTime=500
ez=np.zeros((size))
hy=np.zeros((size))
ez_snap=np.zeros((maxTime))
hy_snap=np.zeros((maxTime))
imp0=377.0
snap_moment=size/2-1
ABC_CELLS=5
eps=4

source_width = 5.0
delay = 10*source_width

i=0

fig = plt.figure()
for qtime in np.arange(0,maxTime+1,1):

    print(qtime)
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2));
    #ABS BOUNDARIES FOR H_FIELD
    hy[size-1:size-1-ABC_CELLS]=hy[size-2-ABC_CELLS]
    for mm in range(0,size-1):

        hy[mm]=hy[mm]+(ez[mm+1]-ez[mm])/imp0
    #ABS BOUNDARIES FOR E_FIELD
    ez[1:1+ABC_CELLS]=ez[2+ABC_CELLS]
    for mm in range(1,size):

        ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0/eps

    i=i+1

    eZ=plt.plot(range(0,size,1),ez,'-',color='b')
    hY=plt.plot(range(0,size,1),hy*imp0,'-',color='g')
    plt.ylabel('Ez-field, V/m')
    plt.xlabel('Spatial')

    plt.show()
    plt.pause(0.00001)
    l1 = eZ.pop(0)
    l1.remove()
    l2 = hY.pop(0)
    l2.remove()