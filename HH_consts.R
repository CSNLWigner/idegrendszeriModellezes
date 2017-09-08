#############################################################
## constants for simulating a Hodgkin-Huxley neuron

## the size of the cell
r <- 1e-3 # 10 um radius for a cell, expressed in cm
vi <- 4*r^3*pi/3 # intracell volume cm3
ve <- vi            # extracell volume - first approximation

a <- 4*r^2*pi # surface area in cm2


## electric properties
Cm <- 1 # membrane capacitance per unit area in  uF/cm2
cm <- Cm * a # membrane capacitance, uF

gK.unit=36 # membrane conductance per unit area - in mS/cm2
gK=gK.unit * a # membrane conductance, mS
gNa.unit=120 # membrane conductance per unit area - in mS/cm2
gNa=gNa.unit * a # membrane conductance, mS
gL.unit=0.3 # membrane conductance per unit area - in mS/cm2
gL=gL.unit * a # membrane conductance, mS

E.Na <- 50 # mV, Na reversal potential
E.K <- -77 # mV, K reversal potential
E.L <- -54.4 # mV, leak reversal potential

I0.default <- 0 # input current in uA
I0 <- I0.default # input current in uA

