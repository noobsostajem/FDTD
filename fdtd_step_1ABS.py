import numpy as np
import math as mt
import matplotlib.pyplot as plt
plt.ion()

size=1800;maxTime=1000
ez=np.zeros((size))
hy=np.zeros((size))
ez_snap=np.zeros((maxTime))
hy_snap=np.zeros((maxTime))
imp0=377.0
snap_moment=size/2-1
ABC_CELLS=5
eps=1
c = 1/np.sqrt(eps)
a = (c-1)/(c+1)
b = 2/(c + 1)


source_width = 5.0
delay = 10*source_width

i=0

# Left boundary
wl_nm1,wl_n,wl_np1 = 0,0,0 # Field at x=0 at time steps n-1, n, n+1
wlp1_nm1,wlp1_n,wlp1_np1 = 0,0,0 # Field at x=1 at time steps n-1, n, n+1
# Right boundary
wr_nm1,wr_n,wr_np1 = 0,0,0 # Field at x=size at time steps n-1, n, n+1
wrm1_nm1,wrm1_n,wrm1_np1 = 0,0,0 # Field at x=size-1 at time steps n-1, n, n+1

fig = plt.figure()
for qtime in np.arange(0,maxTime+1,1):

    print(qtime)
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2));
    #ABS BOUNDARIES FOR H_FIELD
    #hy[size-1:size-1-ABC_CELLS]=hy[size-2-ABC_CELLS]
    for mm in range(0,size):

        hy[mm]=hy[mm]+(ez[mm]-ez[mm])/imp0

        wrm1_np1 = hy[-2]
        wr_np1 = -wrm1_nm1 + a*(wrm1_np1+wr_nm1) + b*(wr_n+wrm1_n)
        hy[-1] = wr_np1

        wr_nm1, wrm1_nm1 = wr_n, wrm1_n
        wr_n, wrm1_n = wr_np1, wrm1_np1
    #ABS BOUNDARIES FOR E_FIELD
        #ez[1:1+ABC_CELLS]=ez[2+ABC_CELLS]
        ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0/eps
        #ez[source_x] += source(time, delay, source_width)
    #Evaluate Mur ABC value (eq. 6.35 Taflove)
        wlp1_np1 = ez[1]
        wl_np1 = -wlp1_nm1 + a*(wlp1_np1+wl_nm1) + b*(wl_n+wlp1_n)
        ez[0] = wl_np1
    #Cycle field values at boundary
        wl_nm1, wlp1_nm1 = wl_n, wlp1_n
        wl_n, wlp1_n = wl_np1, wlp1_np1        

    i=i+1

    eZ=plt.plot(range(0,size,1),ez,'-',color='b')
    hY=plt.plot(range(0,size,1),hy*imp0,'-',color='g')
    plt.ylabel('Ez-field, V/m')
    plt.xlabel('Spatial')

    plt.show()
    plt.pause(0.0000000000001)
    l1 = eZ.pop(0)
    l1.remove()
    l2 = hY.pop(0)
    l2.remove()
