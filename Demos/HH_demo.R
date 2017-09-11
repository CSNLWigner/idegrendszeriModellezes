#############################################################
## demo for simulating the Hodgkin-Huxley equations

source('HH_sim.R')

cat('simulating the Hodgkin-Huxley equations\n')
cat('there are two simulations in this demo:\n')
cat('1: constant current\n')
cat('2: current pulses\n')
id.sim <- readline('press [1] or [2] to choose between them, or press anything else to quit \n')
i.run <- 'Y'


###########################################
## simulation with constant input
cols <- c(2,3,4)

if (id.sim == '1'){
	cat('Constant current is injected in a neuron of 20 um diameter with Hodgkin-Huxley channels\n')
	cat('Try to find the AP threshold!\n\n')

	I0 <- 0.05/1000
	while (i.run  %in% c('Y', 'y')){

		# simulation
		times <-  seq(0,50, by=0.02) # in ms!
		input <- cbind(times, rep(I0, length(times)))
		I.ext <- approxfun(input[,1], input[,2], method = "linear", rule = 2)
		
		params <- c(gK=gK, gNa=gNa, gL=gL, cm=cm, E.Na=E.Na, E.K=E.K, E.L=E.L)
		state <- c(v=-65, m=.053, h=.596, n=.317)
		
		out <- ode(y = state, times = times, func = sim.HH, parms = params)
	
		# plotting		
		layout(matrix(c(1,2,3,4,6,5), 3), c(1,1), c(4, 2, 4))
		par(mar=c(1,4,1,1))
		plot(out[,1], out[,2], ylim=c(-80, 40), axes=F, xlab="", ylab="membrane potential (mV)", t="l"); axis(2, las=2)
		plot(out[,1], out[,6] * 1000, axes=F, xlab="", ylab="stimulus (nA)", t="l"); axis(2, las=2)

		par(mar=c(4,4,1,1))
		matplot(out[,1], out[,3:5], ylim=c(0, 1), axes=F, xlab="time (ms)", ylab="probability", t="l", lty=1, col=cols); axis(2, las=2)
		axis(1)
		legend('topright', c('m-gate', 'h-gate', 'n-gate'), col=cols, lty=1, bty="n")
		
		
		par(mar=c(4,4,4,4))
		plot(out[,2], out[,3], pch=15, cex=0.3, axes=F, xlab="voltage (mV)", ylab="m-gate")
		axis(1); axis(2, las=2)

		plot(out[,2], out[,5], pch=15, cex=0.3, axes=F, xlab="voltage (mV)", ylab="h-gate")
		axis(1); axis(2, las=2)

		# parameters
		i.run <- readline('do you want to rerun with different current [Y or N]?')
		if (i.run %in% c('Y', 'y')){
			I0.new <- as.numeric(readline(paste('set new value for I0 (default = ', I0.default*1000, 'nA, range: [-0.1, 1]):')))
			if (!(is.na(I0.new))) I0 <- I0.new / 1000
		}
	}
}


###########################################
## simulation with a train of excitatory pulses

if (id.sim == '2'){
	cat('Short (2 ms) current pulse is injected in a neuron of 20 um diameter with Hodgkin-Huxley Na and K channels\n')
	cat('Try to find the AP threshold!\n')
	cat('Try also negative currents, or trains of pulses!\n\n')

	I.stim <- 1/20 / 1000
	while (i.run  %in% c('Y', 'y')){

		# simulation
		times <-  seq(0,50, by=0.02) # in ms!
		input <- cbind(times, rep(I0, length(times)))
		n.pulses <- 1
		delay <- 10
		input <- set.input(delay, delay+2, I.stim, input)
						
		I.ext <- approxfun(input[,1], input[,2], method = "linear", rule = 2)
		params <- c(gK=gK, gNa=gNa, gL=gL, cm=cm, E.Na=E.Na, E.K=E.K, E.L=E.L)
		state <- c(v=-65, m=.053, h=.596, n=.317)

		out <- ode(y = state, times = times, func = sim.HH, parms = params)

		## plotting
		layout(matrix(c(1,2), 2), 1, c(5, 3))
		par(mar=c(1,4,1,1))
		plot(out[,1], out[,2], ylim=c(-80, 40), axes=F, xlab="", ylab="membrane potential (mV)", t="l"); axis(2, las=2)
		par(mar=c(4,4,1,1))
		plot(out[,1], out[,6] * 1000, axes=F, xlab="time (ms)", ylab="stim (nA)", t="l"); axis(2, las=2)
		axis(1)		

		# parameters
		i.run <- readline('do you want to rerun with different current [Y or N]?')
		if (i.run %in% c('Y', 'y')){
			I.stim.new <- as.numeric(readline(paste('set new value for I.stim (default = 1/20, nA, range: [-0.2, 1]):')))
			if (!(is.na(I.stim.new))) I.stim <- I.stim.new / 1000
		}
	}
}

		
###########################################
## simulation with short pulses
## 1 pulse, 2 ms duration: threshold is around -1/5 and 1/20 nA
## 3 E pulse: optimal delay is 17 ms, 1/22.4 nA
## I-E train: optimal delay is 9 ms, -1/10 nA followed by 1/25 nA

# amps <- c(-1/5, -1/10, -1/20, 1/100, 1/80, 1/60, 1/40, 1/30, 1/20, 1/10) / 1000
# resps <- matrix(NA, 2501, 10)

# for (i.amp in 1:10){
	# times <-  seq(0,50, by=0.02) # in ms!
	# I0 <- 0
	# input <- cbind(times, rep(I0, length(times)))
	# input <- set.input(10,12, amps[i.amp], input)
	
	# I.ext <- approxfun(input[,1], input[,2], method = "linear", rule = 2)
	# params <- c(gK=gK, gNa=gNa, gL=gL, cm=cm, E.Na=E.Na, E.K=E.K, E.L=E.L)
	# state <- c(v=-65, m=.053, h=.596, n=.317)
	
	# out <- ode(y = state, times = times, func = sim.HH, parms = params)
	# resps[,i.amp] <- out[,2]
# }
# matplot(resps, t="l", col=rainbow(10, end=0.7), lty=1, xlab="time (ms)", ylab="voltage (mV)")
