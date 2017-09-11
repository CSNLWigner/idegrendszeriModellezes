#############################################################
## demo for simulating the buildup of the equilibrium potential in a simple neuron

source('ions_sim.R')

cat('simulating the buildup of the equilibrium potential in a simple neuron\n')
cat('both electric and chemical changes are simulated \n')
cat('the goal is to compare the time course and magnitude of the changes under different conditions \n\n')
cat('there are two simulations in this demo:\n')
cat('1: Nernst potential\n')
cat('2: equilibrium potential\n')
id.sim <- readline('press [1] or [2] to choose between them, or press anything else to quit \n')
i.run <- 'Y'


#########################################
## Nernst equation

if (id.sim == '1'){
	cat('We will simulate the flux of K ions through a simple membrane permeable only to K ions\n')
	cat('initially there are no electric potential difference between the two sides of the membrane\n')
	cat('but there is a difference in the ionic concentrations\n')

	while (i.run  %in% c('Y', 'y')){
		
		## simulation
		params <- c(gK=gK, cm=cm, vi=vi, ve=ve)
		state <- c(v=0, C.Ki=C.Ki.init, C.Ke=C.Ke.init)
		
		times <-  seq(0,30, by=1/10)
		out <- ode(y = state, times = times, func = sim.Nernst, parms = params)
		# diagnostics(out)

		## plotting
		icol <- rainbow(24)[sample(1:24, 1)]
		layout(matrix(1:3, 3), 1, c(2,2,3))
		par(mar=c(1,5,3,1))
		plot(out[,1], out[,2], t="l", col=icol, axes=F, xlab="", ylab="", main="membrane potential (mV)"); axis(2, las=2)
		par(mar=c(1,5,3,1))
		plot(out[,1], out[,3], t="l", col=icol, axes=F, xlab="", ylab="", main="[K]-intra (mM)"); axis(2, las=2)
		par(mar=c(4,5,3,1))
		plot(out[,1], out[,4], t="l", col=icol, axes=F, xlab="time (ms)", ylab="", main="[K]-extra (mM)"); axis(2, las=2)
		axis(1)
	
		## parameters
		i.run <- readline('do you want to rerun it with different parameters [Y or N]?')
		if (i.run %in% c('Y', 'y')){
			gK.new <- as.numeric(readline(paste('set new value for gK (default = ', gK.default, '):')))
			if (!(is.na(gK.new))) gK <- gK.new else gK <- gK.default
			cm.new <- as.numeric(readline(paste('set new value for cm (default = ', cm.default, '):')))
			if (!(is.na(cm.new))) cm <- cm.new	else cm <- cm.default
			C.Ki.init.new <- as.numeric(readline(paste('set new value for C.Ki.init (default = ', C.Ki.default, '):')))
			if (!(is.na(C.Ki.init.new))) C.Ki.init <- C.Ki.init.new else C.Ki.init <- C.Ki.default
		}
	}
}


#########################################
## Equilibrium

if (id.sim == '2'){
	cat('We will simulate the buildup of the resting potantial of a neuron.\n')
	cat('Initially there are no electric or chemical gradients.\n')
	cat('At t=0, we switch on the Na/K exchanger, that starts transporting ions.\n')
	cat('Both Na and K ions are simulated, the permeability is assumed to be constant.\n')

	while (i.run  %in% c('Y', 'y')){

		# simulation
		params <- c(gK=gK, gNa=gNa, cm=cm, vi=vi, ve=ve, I.pump.K=I.pump.K, I.pump.Na=I.pump.Na)
		state <- c(v=0, C.Nai=70, C.Nae=70, C.Ki=65, C.Ke=65)
		
		times <-  seq(0, 5000000, by=1000) # in ms!
		out <- ode(y = state, times = times, func = sim.equilibrium, parms = params)


		# plotting
		layout(matrix(1:2, 2), 1, c(4,5))
		par(mar=c(1,5,3,1))
		plot(out[,1], out[,2], t="l", col=1, axes=F, xlab="", ylab="", main="membrane potential (mV)"); axis(2, las=2)
		par(mar=c(4,5,3,1))
		matplot(out[,1], out[,3:6], t="l", col=c(2,2,3,3), lty=c(1,2,1,2), axes=F, xlab="time (s)", ylab="", main="ion concentrations (mM)"); axis(2, las=2)
		axis(1, seq(0, 5000000, by= 1000000), seq(0, 5000, by= 1000))
		legend('topright', c('[Na] intra', '[Na] extra', '[K] intra', '[K] extra'), col=c(2,2,3,3), lty=c(1,2,1,2), bg=grey(1, alpha=0.75))


		## parameters
		i.run <- readline('do you want to rerun it with different parameters [Y or N]?')
		if (i.run %in% c('Y', 'y')){
			gK.new <- as.numeric(readline(paste('set new value for gK (default = ', gK.default, '):')))
			if (!(is.na(gK.new))) gK <- gK.new else gK <- gK.default
			
			gNa.new <- as.numeric(readline(paste('set new value for gNa (default = ', gNa.default, '):')))
			if (!(is.na(gNa.new))) gNa <- gNa.new else gNa <- gNa.default
		}
		
	}	
}



