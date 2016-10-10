import numpy as np
import math as mt
import matplotlib.pyplot as plt

#GROUND CONDITIONS###########
epsilon_1 = 9
epsilon_2=9
imp0=377.0
size=500;
eps= np.ones(size)
eps[0:size/2] = epsilon_1
eps[size/2:size] = epsilon_2
source_width=10.0*np.sqrt(max(eps))
delay = 10*source_width

maxTime=((size/2)*np.sqrt(epsilon_1)-int(0.05*size)+delay)+max(2*source_width/np.sqrt(max(eps))*np.sqrt(epsilon_1),2*source_width/np.sqrt(max(eps))*np.sqrt(epsilon_2))

#############################
ez=np.zeros((size)) #Array's for ez and hy
hy=np.zeros((size))
#Source Parameters###########

#maxTime=int(size*3.0+delay)
print(maxTime)
ez=np.zeros((size)) #Array's for ez and hy
hy=np.zeros((size))
#############################
E_monitor_1=int(size*0.3)
E_monitor_2=int(size*0.7)
H_monitor_1=int(size*0.3)
H_monitor_2=int(size*0.7)
#CPML Constants##############
N_cell=10.0
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
    ez[int(0.05*size)] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2))

    E_monitor_1 = max(0, (max(np.abs(ez[0:int(size/2-1)]))))
    E_monitor_2 = max(0, (max(np.abs(ez[int(size/2-1):size]))))
    H_monitor_1 = max(0, (max(np.abs(hy[0:int(size/2-1)]))))
    H_monitor_2 = max(0, (max(np.abs(hy[int(size/2-1):size]))))
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

#PLOT GRAPH#####################

res_fresnel=1-np.abs((np.sqrt(epsilon_1)-np.sqrt(epsilon_2))/(np.sqrt(epsilon_1)+np.sqrt(epsilon_2)))**2
res_fdtd=(E_monitor_2*H_monitor_2)/(E_monitor_1*H_monitor_1+E_monitor_2*H_monitor_2)
#res_fdtd=(E_monitor_2*H_monitor_2)/(E_monitor_1*H_monitor_1)
error=np.abs((res_fdtd-res_fresnel)/res_fresnel)
################################PRINT#####
print("n1 = %f, n2 = %f || FDTD = %f, Fresnel = %f, Error = %f%%" % (np.sqrt(epsilon_1), np.sqrt(epsilon_2),res_fdtd,res_fresnel,error*100))