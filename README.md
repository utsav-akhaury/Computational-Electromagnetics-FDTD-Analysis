# Computational Electromagnetics (FDTD Analysis)

Codes for implementing the Finite Difference Time Domain (FDTD) method to solve the Maxwell's equations for modelling different electromagnetic structures.

Requires the following libraries/toolkits for the Python codes- 

1) **NumPy** &nbsp;- &nbsp;for mesh discretization

2) **matplotlib** &nbsp;- &nbsp;for general plots

3) **mpl_toolkits.mplot3d** &nbsp;- &nbsp;for surface plots <br /><br />

## Description of the Models Implemented <br />

1) [homog_gau.py](https://github.com/utsav-akhaury/Computational-Electromagnetics-FDTD-Analysis/blob/master/homog_gau.py)
Gaussian pulse propagation in a homogeneous medium

2) [homog_gau_abc.py](https://github.com/utsav-akhaury/Computational-Electromagnetics-FDTD-Analysis/blob/master/homog_gau_abc.py)
Gaussian pulse propagation in a homogeneous medium terminated with Absorbing Boundary Condition (ABC)

3) [interface_gau.py](https://github.com/utsav-akhaury/Computational-Electromagnetics-FDTD-Analysis/blob/master/interface_gau.py)   
Gaussian Pulse propagation through an interface with ABC at extreme boundaries

4) [interface_sin.py](https://github.com/utsav-akhaury/Computational-Electromagnetics-FDTD-Analysis/blob/master/interface_sin.py)  
Sine wave propagation through an interface with ABC at extreme boundaries

5) [Microstrip_SW.py](https://github.com/utsav-akhaury/Computational-Electromagnetics-FDTD-Analysis/blob/master/Microstrip_SW.py)  
50 Ω Microstrip transmission line with an open termination (ZL = inf) and ABC

6) ***Microstrip_SW 2.py***  
Two symmetrically spaced 50 Ω Microstrip transmission lines with open termination (ZL = inf) and ABC

7) ***Microstrip_SW 2 dielectric.py***  
Two symmetrically spaced 50 Ω Microstrip transmission lines with open termination (ZL = inf), a cylindrical dielectric region in between the striplines and feed on one strip

8) ***Two_port_sin.m***  
Two-port transmission line with sinusoidal feed at transmitting end. The corresponing signal at the receiving end is analysed in time & frequency domain

9) ***Two_port_gaussian.m***  
Two-port transmission line with gaussian feed at transmitting end. The corresponing signal at the receiving end is analysed in time & frequency domain

10) ***Two_port_modulated_gaussian.m***  
Two-port 50 Ω transmission line with gaussian signal modulated by a sinusoidal pulse at transmitting end

11) ***Two_port_Luebbers_Tapered_sin.m***  
Two-port 50 Ω transmission line with Luebber's source and Source Tapering (staircased FDTD mesh transition from ground plane to stripline) for a sinusoidal signal

12) ***Two_port_Luebbers_Tapered_gaussian.m***  
Two-port 50 Ω transmission line with Luebber's source and Source Tapering for a gaussian signal

13) ***Two_port_Luebbers_Tapered_mod_gaussian.m***  
Two-port 50 Ω transmission line with Luebber's source and Source Tapering for a gaussian signal modulated by a sinusoidal pulse
