
#--------------- Two symmetrically spaced 50 Ohm Microstrip transmission lines with open termination (ZL = inf) ------
#--------------- Cylindrical dielectric region in between the 2 striplines and feed on one strip ---------------------
#--------------- To plot standing waves (SW) on the 2 strip lines and observe the coupling ---------------------------


import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os

#----- Substrate and line dimensions for 50 Ohm----------
sub_l = 35e-3
sub_h = 0.794e-3
sub_w = 22e-3

strip_l = 2*sub_l/3
strip_w = 2.46e-3

#----- Substrate material and surrounding air ----------
eps0 = 8.854e-12
meu0 = 4*np.pi*1e-7
epsr = 2.2                   # RT/Duroid material
meur = 1
epss = eps0*epsr             # substrate
meus = meu0*meur
epsa = eps0                  # air
meua = meu0
epse = (epsa + epss)/2       # Average eps

#---- EM wave -----------
c = 3e8
vmin = c/np.sqrt(epsr*meur)
vmax = c
fmax = 15e9
lamda = vmin/fmax

#------ FDTD cell length and time step---------
#------- Substrate thickness is in Y-direction
dx = lamda/35
dy = dx/2
dz = dx 

dt = 1/(vmax*np.sqrt(1/np.power(dx,2) + 1/np.power(dy,2) + 1/np.power(dz,2)))

xoff = 4
yoff = 20
zoff = 10

Nx = int(round(sub_w/dx)) + xoff - 1
Ny = int(round(sub_h/dy)) + yoff                     # upward
Nz = int(round(sub_l/dz)) + zoff

#----------- Start and end coordinates of the subsrate 
sub_x1 = int(round(xoff/2))
sub_x2 = int(round(Nx-sub_x1))
sub_y1 = 1                                # thickness starts from index 1
sub_y2 = int(round(sub_h/dy)) + sub_y1
sub_z1 = int(round(zoff/2))
sub_z2 = int(round(Nz-sub_z1))

#----------- Start and end coordinates of strip 1
Nw = int(round(strip_w/dx))
strip_x11 = int(round(0.25*(Nx))-(Nw/2))

strip_x12 = strip_x11 + Nw
strip_y1 = sub_y1 + 1                            # strip is on the top surface of the substrate
strip_z11 = sub_z1 + 5                           # strip is starting 5 cells into the substrate
strip_z12 = strip_z11 + int(round(strip_l/dz))

#----------- Start and end coordinates of strip 2
strip_x21 = int(round(0.75*(Nx))-(Nw/2))

strip_x22 = strip_x21 + Nw
strip_y2 = sub_y1 + 1                            # strip is on the top surface of the substrate
strip_z21 = sub_z1 + 5                           # strip is starting 5 cells into the substrate
strip_z22 = strip_z21 + int(round(strip_l/dz))

dia = strip_x21-strip_x12                      #Circle 
rad = dia/2
x_centre = int(round(0.5*Nx))
z_centre = int(round(45 + rad))

#-------- Feed point of the strip 1
feed_x1 = strip_x11 + int(round(0.5*strip_w/dx))       # center of strip
feed_y1 = sub_y1
feed_z11 = strip_z11                              # feed starts at the strip edge
feed_z12 = strip_z12

#-------- Feed point of the strip 2
feed_x2 = strip_x21 + int(round(0.5*strip_w/dx))       # center of strip
feed_y2 = sub_y1
feed_z21 = strip_z21                              # feed starts at the strip edge
feed_z22 = strip_z22

#------- Initialise E,H,multipliers and ABC arrays ------------
ex = np.zeros((Nx,Ny,Nz))
ey = np.zeros((Nx,Ny,Nz))
ez = np.zeros((Nx,Ny,Nz))

hx = np.zeros((Nx,Ny,Nz))
hy = np.zeros((Nx,Ny,Nz))
hz = np.zeros((Nx,Ny,Nz))

const_ex = np.zeros((Nx,Ny,Nz))
const_ey = np.zeros((Nx,Ny,Nz))
const_ez = np.zeros((Nx,Ny,Nz))

const_hx = np.zeros((Nx,Ny,Nz))
const_hy = np.zeros((Nx,Ny,Nz))
const_hz = np.zeros((Nx,Ny,Nz))

#------- Default constants in update eqn
const_ex[0:Nx,0:Ny,0:Nz] = dt/(epsa*dx)
const_ey[0:Nx,0:Ny,0:Nz] = dt/(epsa*dy)
const_ez[0:Nx,0:Ny,0:Nz] = dt/(epsa*dz)

#--------------------------------------
const_hx[0:Nx,0:Ny,0:Nz] = dt/(meua*dx)
const_hy[0:Nx,0:Ny,0:Nz] = dt/(meua*dy)
const_hz[0:Nx,0:Ny,0:Nz] = dt/(meua*dz)

#----------- Substrate---------
const_ex[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epss*dx) 
const_ey[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epss*dy)
const_ez[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epss*dz)

#----------- Top substrate-air interface has eps_avg
const_ex[sub_x1-1:sub_x2, strip_y1-1, sub_z1-1:sub_z2] = dt/(epse*dx)
const_ez[sub_x1-1:sub_x2, strip_y1-1, sub_z1-1:sub_z2] = dt/(epse*dz)

#----------- Side substrate-air interface has eps_avg 
const_ex[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z1-1] = dt/(epse*dx)
const_ey[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z1-1] = dt/(epse*dy)
const_ex[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z2-1] = dt/(epse*dx)
const_ey[sub_x1-1:sub_x2, sub_y1-1:sub_y2, sub_z2-1] = dt/(epse*dy)

const_ey[sub_x1-1, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epse*dy)
const_ez[sub_x1-1, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epse*dz)
const_ey[sub_x2-1, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epse*dy)
const_ez[sub_x2-1, sub_y1-1:sub_y2, sub_z1-1:sub_z2] = dt/(epse*dz)

#----------- Metal boundaries have Etan=0 ---------
#------- Constants in E-update eqn
gx = np.ones((Nx,Ny,Nz))
gy = np.ones((Nx,Ny,Nz))
gz = np.ones((Nx,Ny,Nz))
fx = np.ones((Nx,Ny,Nz))
fy = np.ones((Nx,Ny,Nz))

#---- GND plane
gx[sub_x1-1:sub_x2, sub_y1-1, sub_z1-1:sub_z2] = 0
gz[sub_x1-1:sub_x2, sub_y1-1, sub_z1-1:sub_z2] = 0

#---- Stripline
gx[strip_x11-1:strip_x12, strip_y1-1, strip_z11-1:strip_z12] = 0
gz[strip_x11-1:strip_x12, strip_y1-1, strip_z11-1:strip_z12] = 0

gx[strip_x21-1:strip_x22, strip_y2-1, strip_z21-1:strip_z22] = 0
gz[strip_x21-1:strip_x22, strip_y2-1, strip_z21-1:strip_z22] = 0

for x in range(strip_x12, strip_x21):
    val1 = z_centre - round((np.sqrt((rad*rad)-(np.power((x-x_centre),2)))))
    val2 = z_centre + round((np.sqrt((rad*rad)-(np.power((x-x_centre),2)))))
    for z in range(int(val1),int(val2)+1):
        gx[x,strip_y1-1,z] = dt/(100*eps0*dx)
        gz[x,strip_y1-1,z] = dt/(100*eps0*dz)
        gx[x,strip_y2-1,z] = dt/(100*eps0*dx)
        gz[x,strip_y2-1,z] = dt/(100*eps0*dz)

#------ ABC constants
ex_abc = np.zeros((Nx,Ny,Nz))
ey_abc = np.zeros((Nx,Ny,Nz))
ez_abc = np.zeros((Nx,Ny,Nz))

constz_abc = (vmax*dt-dz)/(vmax*dt+dz)
consty_abc = (vmax*dt-dy)/(vmax*dt+dy)
constx_abc = (vmax*dt-dx)/(vmax*dt+dx)

#------ Plot gx at the surface of the substrate to verify the model
arr = np.copy(gx[:,strip_y1-1,:])
xdata = np.linspace(0,Nz-1,Nz, dtype=int)
ydata = np.linspace(0,Nx-1,Nx, dtype=int)
zdata = np.copy(arr[0:Nx, 0:Nz])
XDAT, YDAT = np.meshgrid(xdata, ydata)

fig = plt.figure()
ax = plt.axes(projection="3d")

ax.plot_wireframe(XDAT, YDAT, zdata)
plt.xlim(0,Nz)
plt.ylim(0,Nx)
plt.show()

#------ Gaussian pulse ------------
Ts = 10*dt                          # pulse width
t0 = 3*Ts                           # delay
N_steps = 1                         # default iteration steps
count = 0

voltage_in = np.zeros((1000000,1))
current_in = np.copy(voltage_in)
voltage_out = np.copy(voltage_in)
current_out = np.copy(voltage_in)
#current_in_corr = current_in

#************** Iteration loop ****************
while(int(round(N_steps))>0) : 
    N_steps = int(input('\nTime steps(0 to quit): '))
     
    for n in range (int(round(N_steps))):        
        time = count*dt
        count = count + 1

        print('\n\tTime step : ', count)  
        os.system('cls')
       
        # ----- Sinusoidal signal of 10 GHz---------
        pulse = np.sin(2*np.pi*10e9*time)
        ey[strip_x11-1:strip_x12, sub_y1-1:sub_y2, feed_z11-1] = pulse/dy   # Feed on strip 1
         
        #------------------------ compute H ------------------------- 
        for j in range(Ny-1):
            for k in range(Nz-1):      
                hx[:,j,k] = hx[:,j,k] + fx[:,j,k]*(const_hz[:,j,k]*(ey[:,j,k+1]-ey[:,j,k])-const_hy[:,j,k]*(ez[:,j+1,k]-ez[:,j,k]))    
        
        for i in range(Nx-1):
            for k in range(Nz-1):
                hy[i,:,k] = hy[i,:,k] + fy[i,:,k]*(const_hx[i,:,k]*(ez[i+1,:,k]-ez[i,:,k])-const_hz[i,:,k]*(ex[i,:,k+1]-ex[i,:,k]))
        
        for i in range(Nx-1):
            for j in range(Ny-1):
                hz[i,j,:] = hz[i,j,:] + const_hy[i,j,:]*(ex[i,j+1,:]-ex[i,j,:]) - const_hx[i,j,:]*(ey[i+1,j,:]-ey[i,j,:])
              
        #------------------------- compute E ------------------------      
        
        for j in range(1,Ny-1):
            for k in range(1,Nz-1):
                ex[:,j,k] = ex[:,j,k] + gx[:,j,k]*(const_ey[:,j,k]*(hz[:,j,k]-hz[:,j-1,k])-const_ez[:,j,k]*(hy[:,j,k]-hy[:,j,k-1]))  
               
        for i in range(1,Nx-1):
            for k in range(1,Nz-1):
                ey[i,:,k] = ey[i,:,k] + gy[i,:,k]*(const_ez[i,:,k]*(hx[i,:,k]-hx[i,:,k-1])-const_ex[i,:,k]*(hz[i,:,k]-hz[i-1,:,k])) 
              
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                ez[i,j,:] = ez[i,j,:] + gz[i,j,:]*(const_ex[i,j,:]*(hy[i,j,:]-hy[i-1,j,:])-const_ey[i,j,:]*(hx[i,j,:]-hx[i,j-1,:]))              
        
        #----------------- ABC---------------------  
        #---------- XY plane 
        ex[:,:,0] = ex_abc[:,:,1] + constz_abc*(ex[:,:,1]-ex[:,:,0]) 
        ey[:,:,0] = ey_abc[:,:,1] + constz_abc*(ey[:,:,1]-ey[:,:,0])        
                        
        ex[:,:,Nz-1] = ex_abc[:,:,Nz-2] + constz_abc*(ex[:,:,Nz-2]-ex[:,:,Nz-1])
        ey[:,:,Nz-1] = ey_abc[:,:,Nz-2] + constz_abc*(ey[:,:,Nz-2]-ey[:,:,Nz-1])     
        
        #---------- XZ plane       
        # bottom plane (index 1) is already defined GND        
        ex[:,Ny-1,:] = ex_abc[:,Ny-2,:] + consty_abc*(ex[:,Ny-2,:]-ex[:,Ny-1,:])
        ez[:,Ny-1,:] = ez_abc[:,Ny-2,:] + consty_abc*(ez[:,Ny-2,:]-ez[:,Ny-1,:]) 
         
        #---------- YZ plane
        ey[0,:,:] = ey_abc[1,:,:] + constx_abc*(ey[1,:,:]-ey[0,:,:])
        ez[0,:,:] = ez_abc[1,:,:] + constx_abc*(ez[1,:,:]-ez[0,:,:])       
                      
        ey[Nx-1,:,:] = ey_abc[Nx-2,:,:] + constx_abc*(ey[Nx-2,:,:]-ey[Nx-1,:,:])
        ez[Nx-1,:,:] = ez_abc[Nx-2,:,:] + constx_abc*(ez[Nx-2,:,:]-ez[Nx-1,:,:])  
       
        ex_abc = np.copy(ex)
        ey_abc = np.copy(ey)  
        ez_abc = np.copy(ez) 
        
    #---------- Plot Ey on the strips -----------------
    fig1 = plt.figure()
    e1 = np.zeros(Nz)
    e1[:] = abs(ey[feed_x1+2, strip_y1-1, :])
    k1 = np.linspace(strip_z11-1, strip_z12-1, strip_z12-strip_z11+1, dtype = int)   
    plt.plot(k1,e1[k1])
    plt.xlim(strip_z11-1, strip_z12)
    plt.ylim(0)
    plt.title('Ey on strip 1')   
    plt.show()

    fig2 = plt.figure()
    e2 = np.zeros(Nz)
    e2[:] = abs(ey[feed_x2+2, strip_y2-1, :])
    k2 = np.linspace(strip_z21-1, strip_z22-1, strip_z22-strip_z21+1, dtype = int)   
    plt.plot(k2,e2[k2])
    plt.xlim(strip_z21-1, strip_z22)
    plt.ylim(0)
    plt.title('Ey on strip 2')   
    plt.show()

