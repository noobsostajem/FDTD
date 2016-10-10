import numpy as np
import math as mt
import matplotlib.pyplot as plt


pld=0 ###### 1 - to show plot

    #GROUND CONDITIONS###########
size=500
imp0=377.0
epsilon_free = 2
epsilon_2=5
wavelength = int(size/5.0) #in host media
c_free=1/mt.sqrt(epsilon_free)
factor = 2 # for slab lambda/2
n2 = 1.0/int(wavelength*mt.sqrt(epsilon_free)/mt.sqrt(epsilon_2)/factor)*wavelength*mt.sqrt(epsilon_free)/factor
epsilon_2=n2**2
print(epsilon_2)
plate_thickness=wavelength*mt.sqrt(epsilon_free)/mt.sqrt(epsilon_2)/factor
eps= np.ones(size)
eps[:]= epsilon_free
eps[int(size/2.0)-int(plate_thickness/2):int(size/2.0)+int(plate_thickness)/2] = epsilon_2
eps1=np.ones(size)
eps1[:]=epsilon_free


    #############################

ez=np.zeros((size)) #Array's for ez and hy
hy=np.zeros((size))

#Source Parameters###########
source_width = wavelength*3.0*np.sqrt(max(eps))
delay = 10*source_width
source_x = int(1.0*size/10.0)  #Source position
maxTime=int(((size)*np.sqrt(epsilon_free)-int(0.05*size)+delay))
#maxTime=int(size*3.0+delay)
print("FDTD has launched: total time steps = %d"%(maxTime))
ez=np.zeros((size)) #Array's for ez and hy
hy=np.zeros((size))
ez_free=np.zeros((size)) #Array's for ez and hy
hy_free=np.zeros((size))
#############################

#CPML Constants##############
N_cell=20.0
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
Phx_free=np.zeros(size)
Pex_free=np.zeros(size)
sx[mm+1] = sxmax*((N_cell-mm-0.5)/N_cell)**g
sxm[mm] = sxmax*((N_cell-mm)/N_cell)**g
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
plt.axvline(size/2-int(plate_thickness/2), color='r',linewidth=3)
plt.axvline(size/2+int(plate_thickness/2), color='r',linewidth=3)
def source(current_time, delay, source_width):
    amp =  mt.exp(-(current_time-delay)**2/(2.0 * source_width**2))
    if current_time > delay:
        amp = 1.0
    return amp/np.sqrt(epsilon_free)*np.sin(2*np.pi*current_time*c_free/wavelength)
for qtime in np.arange(0,maxTime+1,1):
    Phx[x] = bhx[x]*Phx[x] + ahx[x]*(ez[x+1] - ez[x])
    Phx_free[x] = bhx[x]*Phx_free[x] + ahx[x]*(ez_free[x+1] - ez_free[x])
    hy[x] +=((ez[x+1] - ez[x]) + Phx[x])/imp0
    hy_free[x] +=((ez_free[x+1] - ez_free[x]) + Phx_free[x])/imp0
    Pex[x+1] = bex[x+1]*Pex[x+1] + aex[x+1]*(hy[x+1]-hy[x])
    Pex_free[x+1] = bex[x+1]*Pex_free[x+1] + aex[x+1]*(hy_free[x+1]-hy_free[x])
    ez[x+1] +=((hy[x+1]-hy[x]) +Pex[x+1])*imp0/eps[x]
    ez_free[x+1] +=((hy_free[x+1]-hy_free[x]) +Pex_free[x+1])*imp0/eps1[x]
    ez[int(0.05*size)] += source(qtime, delay, source_width)
    ez_free[int(0.05*size)] += source(qtime, delay, source_width)
#########Monitors###########
    ymin=1.1*min([min(ez),min(hy*imp0)])
    ymax=1.1*max([max(ez),max(hy*imp0)])
    ymin_free=1.1*min([min(ez_free),min(hy_free*imp0)])
    ymax_free=1.1*max([max(ez_free),max(hy_free*imp0)])
    ymin=min(ymin,ymin_free)
    ymax=max(ymax,ymax_free)
###########################
    if pld==1:
        eZ=plt.plot(range(0,size,1),ez,'-',color='b', label="Ez")
        hY=plt.plot(range(0,size,1),ez_free,'-',color='g', label="Ez_free")
        plt.legend(bbox_to_anchor=(0.75, 1.0), loc=2, ncol=1)
        plt.ylabel('Ez-field, V/m')
        plt.xlabel('x, nm')


        ax = fig.add_subplot(111)

        ax.set_title('Time = {0} ats'.format(qtime))
        plt.ylim([ymin,ymax])
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        plt.show()
        plt.pause(0.000000001)
        l1 = eZ.pop(0)
        l1.remove()
        l2 = hY.pop(0)
        l2.remove()



################################PRINT#####
plt.subplot(3, 1, 1)
plt.axvline(size/2-int(plate_thickness/2), color='r',linewidth=2)
plt.axvline(size/2+int(plate_thickness/2), color='r',linewidth=2)
eZ=plt.plot(range(0,size,1),ez,'-',color='b', label="Ez")
hY=plt.plot(range(0,size,1),ez_free,'-',color='g', label="Ez_free")
plt.legend(bbox_to_anchor=(0.8, 1.0), loc=2, ncol=1)
plt.ylabel('Ez-field, V/m')
plt.xlabel('x, nm')

ax = fig.add_subplot(311)

ax.set_title('Time = {0} ats'.format(qtime))
plt.ylim([ymin,ymax])
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.show()
plt.pause(0.000000001)
l1 = eZ.pop(0)
l2 = hY.pop(0)

plt.subplot(3,1,2)
plt.axvline(size/2-int(plate_thickness/2), color='r',linewidth=2)
plt.axvline(size/2+int(plate_thickness/2), color='r',linewidth=2)
eZ=plt.plot(range(0,size,1),ez-ez_free,'-',color='b', label="Ez-Ez_free")
plt.legend(bbox_to_anchor=(0.8, 1.0), loc=2, ncol=1)
plt.ylabel('Ez-field, V/m')
plt.xlabel('x, nm')

ax = fig.add_subplot(312)
ax.set_title('Time = {0} ats'.format(qtime))
plt.ylim([ymin,ymax])
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.show()
plt.pause(0.000000001)
l1 = eZ.pop(0)

plt.subplot(3,1,3)
plt.subplots_adjust(hspace=0.5)
eZ=plt.plot(range(0,size,1),ez-ez_free,'-',color='b', label="Ez-Ez_free")
plt.legend(bbox_to_anchor=(0.8, 1.0), loc=2, ncol=1)
plt.ylabel('Ez-field, V/m')
plt.xlabel('x, nm')
ymin1=0.05*min(ez[0:size/2-int(plate_thickness)-1])
ymax1=0.05*max(ez[0:size/2-int(plate_thickness)-1])

ax = fig.add_subplot(313)
ax.set_title('Time = {0} ats'.format(qtime))
plt.ylim([ymin1,ymax1])
plt.xlim([0,size/2-int(plate_thickness/2)])
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.show()
plt.pause(0.000000001)
l1 = eZ.pop(0)
plt.show()