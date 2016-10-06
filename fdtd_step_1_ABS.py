import numpy as np
import math as mt
import matplotlib.pyplot as plt

#GROUND CONDITIONS###########
size=500;
imp0=377.0
epsilon = 1.0
eps= np.ones(size)
eps[:] = epsilon
dt=1
dx=1
c = 1/np.sqrt(epsilon)
A = (c-dx/dt)/(c+dx/dt)
B = 2/(c + dx/dt)
#Source Parameters###########
source_width=10.0*np.sqrt(epsilon)
delay = 10*source_width
maxTime=int(1*(size+delay)*np.sqrt(epsilon)*10/13.9)
print(maxTime)

#############################
ez=np.zeros((size))
hy=np.zeros((size))
ez_snap=np.zeros((maxTime))
hy_snap=np.zeros((maxTime))
snap_moment=size/2-1

#############################

plt.ion()
fig = plt.figure()
for qtime in np.arange(0,maxTime+1,1):

    #print(qtime)
    #Source inizialization############################
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2))
    ##################################################

    hy[-1]=hy[-2]
    for x in range(0,size-1):
        hy[x]+=(ez[x+1]-ez[x])/imp0

    #########################################

    #ABC BOUNDARIES FOR H_FIELD###

    ez[0]=ez[1]
    for x in range(1,size):
        ez[x]+=(hy[x]-hy[x-1])*imp0/eps[x]

#PLOT GRAPH#####################
    eZ=plt.plot(range(0,size,1),ez,'-',color='b')
    hY=plt.plot(range(0,size,1),hy*imp0,'-',color='g')
    plt.ylabel('Ez-field, V/m')
    plt.xlabel('x, nm')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.set_title('Time = {0} ats'.format(qtime))
    plt.ylim([1.1*min([min(ez),min(hy*imp0)]),1.1*max([max(ez),max(hy*imp0)])])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.show()
    plt.pause(0.000000001)
    l1 = eZ.pop(0)
    l1.remove()
    l2 = hY.pop(0)
    l2.remove()
################################PRINT#####
eZ=plt.plot(range(0,size,1),ez,'-',color='b')
hY=plt.plot(range(0,size,1),hy*imp0,'-',color='g')
plt.ylabel('Ez-field, V/m')
plt.xlabel('x, nm')
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
ax.set_title('Time = {0} ats'.format(qtime))
plt.ylim([1.1*min([min(ez),min(hy*imp0)]),1.1*max([max(ez),max(hy*imp0)])])
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.show()
