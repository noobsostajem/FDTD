import numpy as np
import math as mt
import matplotlib.pyplot as plt

#GROUND CONDITIONS###########
size=500
imp0=377.0
epsilon = 5
eps= np.ones(size)
eps[:] = epsilon
dx=1
#############################

ez=np.zeros((size)) #Array's for ez and hy
hy=np.zeros((size))

#Source Parameters###########
source_width=10.0*np.sqrt(epsilon)
delay = 10*source_width
maxTime=int(1*(size+delay)*np.sqrt(epsilon)*10/13.9)
print(maxTime)

#############################

#CPML Constants##############
N_cell=120.0
cell_width=1.0
R0=1e-5
g=3.4 # Order of polynomial grading
mm = np.arange(0,int(N_cell))
##############Conductivities#
sxmax=-(g+1)*np.log(R0)/2/imp0/(N_cell*cell_width)
sx=np.zeros(size)
sxm=np.zeros(size)
Phx=np.zeros(size)
Pex=np.zeros(size)
sx[mm+1] = sxmax*((N_cell-mm-0.5)/N_cell)**g
sxm[mm] = sxmax*((N_cell-mm)/N_cell)**g  # Shifted to the right
sx[size-mm-1] = sxmax*((N_cell-mm-0.5)/N_cell)**g
sxm[size-mm-1] = sxmax*((N_cell-mm)/N_cell)**g

aex = np.exp(-sx*imp0)-1
bex = np.exp(-sx*imp0)
ahx = np.exp(-sxm*imp0)-1
bhx = np.exp(-sxm*imp0)
x = np.arange(0,size-1,1)
#############################

plt.ion()
fig = plt.figure()
for qtime in np.arange(0,maxTime+1,1):
    Phx[x] = bhx[x]*Phx[x] + ahx[x]*(ez[x+1] - ez[x])
    hy[x] +=((ez[x+1] - ez[x]) + Phx[x])/imp0
    Pex[x+1] = bex[x+1]*Pex[x+1] + aex[x+1]*(hy[x+1]-hy[x])
    ez[x+1] +=((hy[x+1]-hy[x]) +Pex[x+1])*imp0/eps[x]
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2))

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
    plt.pause(0.00001)
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

 #plt.pause(0.5)

