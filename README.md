# Computational Electromagnetics (FDTD Analysis)

Python codes for implementing the Finite Difference Time Domain (FDTD) method to solve the Maxwell's equations for modelling different electromagnetic structures.

Requires the following libraries/toolkits - 

1) **NumPy** &nbsp;- &nbsp;for mesh discretization

2) **matplotlib** &nbsp;- &nbsp;for general plots

3) **mpl_toolkits.mplot3d** &nbsp;- &nbsp;for surface plots <br /><br />

## Description of the Models Implemented <br />

1) ***homog_gau.py***  
Gaussian pulse propagation in a homogeneous medium

2) ***homog_gau_abc.py***   
Gaussian pulse propagation in a homogeneous medium terminated with Absorbing Boundary Condition (ABC)

3) ***interface_gau.py***   
Gaussian Pulse propagation through an interface with ABC at extreme boundaries

4) ***interface_sin.py***  
Sine wave propagation through an interface with ABC at extreme boundaries

5) ***Microstrip_SW.py***  
50 Ohm Microstrip transmission line with an open termination (ZL = inf) and ABC

6) ***Microstrip_SW 2.py***  
Two symmetrically spaced 50 Ohm Microstrip transmission lines with open termination (ZL = inf) and ABC

7) ***Microstrip_SW 2 dielectric.py***  
Two symmetrically spaced 50 Ohm Microstrip transmission lines with open termination (ZL = inf), a cylindrical dielectric region in between the striplines and feed on one strip
