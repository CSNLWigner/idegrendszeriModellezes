# an Izhikevich neuron simulated using the euler method 
# for more details, see:
# Eugene M. Izhikevich: Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting
# https://www.izhikevich.org/publications/dsn.pdf
# Chapter 8.1.4.

sim.Izhikevich <- function(I, dt=1, tau=20, Tmax=1000, Rm=0.2, v.rest=-60, v.init=v.rest, v.reset=-50, v.threshold=-40, v.peak=35, a=0.03, b=-2, d=100){
	## we simulate an integrate and fire neuron using eulers method
	# I: input current, either a constant or a vector, pA
	# dt: timestep, in ms
	# tau: time constant of the integration, ms
	# Tmax: duration of the simulation (ms)
	# Rm: membrane resistance, GOhm
	# v.rest: resting potential, mV
	# v.init: initial value for the voltage
	# v.reset: reset after spike
	# v.threshold: spike threshold
	# v.peak: peak of the spike
	# a,b,c: parameters of the slow variable

	if (length(I) > 1) Tmax <- length(I) else {
		I <- rep(I, Tmax/dt)
		print('I vectorised')
	}
	N <- length(I) # number of data points
	V <- rep(v.init, length(I))
	v <- v.init
	u <- 0
	
	for (n in 1:N){
		dv.dt <- 0.14 * (v-v.rest)*(v-v.threshold)/tau - Rm * u / tau + Rm * I[n] / tau
		du.dt <- a * (b * (v-v.rest)-u)
		# cat('dv.dt', dv.dt, 'du.dt', du.dt, 'v', v, '\n')
		v <- v + dv.dt * dt
		u <- u + du.dt * dt

		V[n] <- v
		if (v > v.peak){
			v <- v.reset
			u <- u + d
		}
		n <- n + 1
	}
	V
}

# stim <- rep(0, 1000); stim[100:1000] <- 70
# v <- sim.Izhikevich(stim)
# plot(v, t="l")