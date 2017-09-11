#############################################################
## functions for simulating the buildup of the equilibrium potential in a simple neuron

isa <- require(deSolve)
if(isa=='FALSE') {
	install.packages('deSolve')
	require(deSolve)
}

#############################################
## definition of constants 

source('ions_consts.R')

###############################################################
## calculates the reversal potential of an ion form the intracellular 
## and extracellular concentrations - Nernst equation

reversal <- function(Ci, Ce, z=1, T=300){
	# Ci, Ce: intracell and extracell concentration
	# z: charge of ion
	# T: temperature in Kelvin
	# R = 8.314 # universal gas constant J / K*mol
	# Faraday = 96500 # Faraday, C/mol
	Vrev <-  log(Ce/Ci) * R*T/(z*Faraday) * 1000 # in mV
	# (J / K*mol) * K  / (C / mol) = J / C = V
	Vrev
}


###############################################################
## simulates the flux of K ions in the presence of both electric and chemical gradients
## implemented as a differential equation with 3 states
## v: intracell voltage, C.Ki: intracell K concentration), V.Ke: extracell K concentration

sim.Nernst <- function(t, state, params){
	with(as.list(c(state, params)),{
		E.K <- reversal(C.Ki,C.Ke)
		# rate of change
		dv <- (gK*(E.K-v))/cm
		# mV/ms = (mS*mV) / uF
		# V/s = uA / uF
		# V/s = C/s / C/V
		# V/s = V/s

		dC.Ki <- (gK*(E.K-v))/(vi * 1000*Faraday)
		dC.Ke <- (-gK*(E.K-v))/(ve * 1000*Faraday)
		# (mM/dm3) / ms = mS * mV / (cm3*1000 * C/M)
		# (mM/dm3) / ms = mS * mV / (mm3 * C/M) 
		# Me-3 /(dm3 * se-3) = Ae-6 * M / (dm3e-6 * C) 
		# M /(dm3 * s) = C/s * M / (dm3 * C) 
		# M /(dm3 * s) = M / (dm3 * s) 
	    	# return the rate of change
 	   list(c(dv, dC.Ki, dC.Ke))		
	})
}

###############################################################
## simulates the flux of Na and K ions in the presence of both electric and chemical gradients
## and an active membrane transpot mechanism (Na/K exchanger)
## implemented as a differential equation with 5 states
## v: intracell voltage, C.Ki: intracell K concentration), V.Ke: extracell K concentration, C.Nai: intracell Na concentration), V.Nae: extracell Na concentration

sim.equilibrium <- function(t, state, params){
	with(as.list(c(state, params)),{
		E.K <- reversal(C.Ki,C.Ke)
		E.Na <- reversal(C.Nai,C.Nae)

		# rate of change
		dv <- (I.pump.Na + I.pump.K + gK*(E.K-v) + gNa*(E.Na-v))/cm
		dC.Nai <- (I.pump.Na + gNa*(E.Na-v)) / (vi*1000*Faraday)
		dC.Nae <- (-1) * (I.pump.Na + gNa*(E.Na-v)) / (ve*1000*Faraday)
		dC.Ki <- (I.pump.K + gK*(E.K-v)) / (vi*1000*Faraday)
		dC.Ke <- (-1) * (I.pump.K + gK*(E.K-v)) / (ve*1000*Faraday)
 	   list(c(dv, dC.Nai, dC.Nae, dC.Ki, dC.Ke))
 	 })
}

###############################################################
## realine text to number
