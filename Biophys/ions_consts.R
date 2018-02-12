#############################################################
## constants for simulating the buildup of the Nernst potential in a simple neuron


## the size of the cell
r <- 1e-3 # 10 um radius for a cell, expressed in cm
vi <- 4*r^3*pi/3 # intracell volume cm3
ve <- vi            # extracell volume - first approximation
a <- 4*r^2*pi # surface area in cm2


## electric properties
Cm <- 1 # membrane capacitance per unit area in  uF/cm2
cm <- Cm * a # membrane capacitance, uF
cm.default <- cm


gK.unit <- 0.36 # membrane conductance per unit area - in mS/cm2
gK <- gK.unit * a # membrane conductance, mS
gK.default <- gK

gNa.unit <- 0.0144 # membrane conductance per unit area - in mS/cm2
gNa <- gNa.unit * a # membrane conductance, mS
gNa.default <- gNa


## some useful constants
e=1.602*1e-19 #charge of an electron C
A=6.022*1e23  #Avogadro number db

Faraday <- 96500 # Faraday, C/mol
R <- 8.314 # universal gas constant J / K*mol


C.Ki.init <- 142.5 # mM, initial intracell K concentration 
C.Ke.init <- 3 # mM, initial EXTRAcell K concentration 
C.Ki.default <- C.Ki.init
C.Ke.default <- C.Ke.init