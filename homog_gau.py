
#------------------ Gaussian pulse propagation in a homogeneous medium

import numpy as np
import matplotlib.pyplot as plt

#----- Medium ----------
length = 2         
eps0 = 8.854e-12
meu0 = 4*np.pi*1e-7
epsr = 1
meur = 1
eps = eps0*epsr
meu = meu0*meur

#---- Signal -----------
c = 3e8
v = c/np.sqrt(epsr*meur)
freq = 3e9
lamda = v/freq

#------ Cell length and time step---------
dz = lamda/10
dt = dz/v
N_cells = int(length/dz)
Mid_cell = int(N_cells/2)

#------ Multiplying constants --------
const_e = dt/(eps*dz)
const_h = dt/(meu*dz)

#------- Initialise E and H arrays ------------
ex = np.zeros(N_cells)
hy = np.zeros(N_cells)

#------ Gaussian pulse ------------
Ts = 10*dt          # pulse width
t0 = 3*Ts           # delay
N_steps = 200       # maximum iteration steps

#************** Iteration loop ****************

for n in range (N_steps):        
    time = n*dt        
         
    #------- Gaussian pulse launched in the middle cell ---------
    pulse = (np.exp(-np.power(((time-t0)/Ts),2)))/dz             
    ex[Mid_cell-1] = pulse             
    
    #------------------------ compute H -------------------------
    k = np.linspace(0, N_cells-2, N_cells-1, dtype = int)
    hy[k] = hy[k]-const_h*(ex[k+1]-ex[k])
        
    #------------------------ compute E -------------------------
    k = np.linspace(1, N_cells-2, N_cells-2, dtype = int)
    ex[k] = ex[k]-const_e*(hy[k]-hy[k-1]) 
        
    #------------------------ plot ------------------------------
    plt.plot(np.linspace(1, N_cells, N_cells, dtype = int),ex)
    plt.xlim(0,200)
    plt.ylim((-150,150))
    plt.grid()
    plt.show()
