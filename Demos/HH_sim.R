## simulating the Hodgkin-Huxley equations

isa <- require(deSolve)
if(isa=='FALSE') {
	install.packages('deSolve')
	require(deSolve)
}

#############################################
## definition of constants 

source('HH_consts.R')

###############################################################
## the famous Hodgkin-Huxley equations - excitability of the squid giant axon

sim.HH <- function(t, state, params){
	with(as.list(c(state, params)),{

		## voltage-dependent opening and closing of the gates, m, h and n
		am <-  .1*(v+40)/(1-exp(-(v+40)/10))
		bm <-  4*exp(-(v+65)/18)
		ah <-  .07*exp(-(v+65)/20)
		bh <-  1/(1+exp(-(v+35)/10))
		an <-  .01*(v+55)/(1-exp(-(v+55)/10))
		bn <-  .125*exp(-(v+65)/80)

		# rate of change
		dv <- (I.ext(t) - gNa*h*(v-E.Na)*m^3-gK*(v-E.K)*n^4-gL*(v-E.L))/cm

		# first order kinetics of the gating variables
		dm <- am*(1-m)-bm*m
		dh <- ah*(1-h)-bh*h
		dn <- an*(1-n)-bn*n
		#units: mV, mS, uF, uA

 	   list(c(dv, dm, dh, dn), stim=I.ext(t))
 	 })
}


###################################################
# this function is used to generate external stimulus (current pulses) for the HH equation

set.input <- function(t1, t2, I, input){
	tt <- input[,1]
	if (t1>t2){
		t3 <- t1
		t1 <- t2
		t2 <- t3
	}
	if (t1 < max(tt)) x1 <- min(which(tt > t1)) else return(input)
	if (t2 > min(tt)) x2 <- max(which(tt < t2)) else return(input) 
	input[x1:x2,2] <- I
	input
}

