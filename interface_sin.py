
#------------------ Sine wave propagation through interface ------------------

import numpy as np
import matplotlib.pyplot as plt

length = 1

#----- Air----------
length_air = 0.5
eps0 = 8.854e-12
meu0 = 4*np.pi*1e-7
epsr = 1
meur = 1
eps = eps0*epsr
meu = meu0*meur

#------- Dielectric -----
epsrd = 4
epsd = eps0*epsrd

#---- Signal -----------
c = 3e8
freq = 3e9
vmax = c/np.sqrt(epsr*meur)
vmin = c/np.sqrt(epsrd*meur)
lamda_min = vmin/freq

#------ Cell length and time step---------
dz = lamda_min/10
dt = dz/vmax
N_cells = int(length/dz)

#------- Initialise arrays ------------
ex = np.zeros(N_cells)
hy = np.zeros(N_cells)
const_e = np.zeros(N_cells)
ex_abc = np.zeros(N_cells)

#------ Multiplying constants--------
const_e[0:N_cells] = dt/(eps*dz)
const_e[99:N_cells] = dt/(epsd*dz)
const_h = dt/(meu*dz)

const_abc = (vmax*dt-dz)/(vmax*dt+dz)
const_abc1 = (vmin*dt-dz)/(vmin*dt+dz)

#------ Sinusoidal pusle parameters------------
Ts = 10*dt          # pulse width
t0 = 3*Ts           # delay
N_steps = 350       # maximum iteration steps

#************** Iteration loop ****************

for n in range (N_steps):        
    time = n*dt
    
    #------- Sinusoidal pulse launched -----------
    pulse = np.sin(2*np.pi*freq*time)/dz            
    ex[4] = ex[4] + pulse
    
    #------------------------ compute H -------------------------
    k = np.linspace(0, N_cells-2, N_cells-1, dtype = int)
    hy[k] = hy[k] - const_h*(ex[k+1]-ex[k])
        
    #------------------------- compute E ------------------------
    k = np.linspace(1, N_cells-2, N_cells-2, dtype = int)
    ex[k] = ex[k] - const_e[k]*(hy[k]-hy[k-1]) 
    
    #------------------------- ABC ------------------------------
    ex[0] = ex_abc[1] + const_abc*(ex[1]-ex[0])
    ex[N_cells-1] = ex_abc[N_cells-2] + const_abc1*(ex[N_cells-2]-ex[N_cells-1])
    ex_abc = np.copy(ex)
   
    #------------------------ plot ------------------------------
    plt.plot(np.linspace(1, N_cells, N_cells, dtype = int),hy)
    plt.xlim(0,200)
    plt.ylim(-0.5,0.5)
    plt.grid()
    plt.show()
    