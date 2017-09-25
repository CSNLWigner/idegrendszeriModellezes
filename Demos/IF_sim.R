# a simple integrate and fire cell, simulated using the euler method 
sim.IF <- function(I, dt=1, tau=20, Tmax=1000, Rm=0.060, v.rest=-65, v.init=v.rest, v.reset=v.rest, v.threshold=v.rest+10, v.spike=v.rest+100){
	# I: input current, either a constant or a vector, pA
	# dt: timestep, in ms
	# tau: time constant of the integration, ms
	# Rm: membrane resistance, GOhm
	# Tmax: duration of the simulation (ms)
	# v.rest: resting potential, mV
	# v.init: initial value for the voltage
	# v.reset: reset after spike
	# v.spike: peak of the spike
	# v.threshold: spike threshold
	
	# returns: the voltage response of the cell
	
	if (is.vector(I)) Tmax <- length(I) else {
		I <- rep(I, Tmax/dt)
	}
	N <- length(I) # number of data points
	V <- rep(v.init, length(I))
	v <- v.init
	
	for (n in 1:N){
		dv.dt <- (v.rest - v)/tau + Rm*I[n]/tau
		v <- v + dv.dt * dt
		V[n] <- v
		if (v > v.threshold){
			V[n] <- v.spike
			v <- v.reset
		}
		n <- n + 1
	}
	V
}


# stim <- rnorm(1000, 150, 200)
# v <- sim.IF(stim)
# plot(v, t="l")