

###############################################
cell.maps.runs <- function(spt, pos, i.runs, dx=0.05, MazeRange, cell.IDs=NULL){
	# spt: matrix of two columns: time, cell.ID
	# pos: matrix of two columns: time, x position
	# i.runs: row index of the start and end of run trials in the pos matrix
	# dx: x resolution of the output matrix
	# MazeRange: start and end coordinates for the maze
	# cell.IDs: cells to be included in the spike.map

	# output: an array with the spikes of each cell on each trial with spatial resolution dx
	# dimensions: cell x distance x trials
	# the last cell is not a true cell, but the time (in seconds) the rat spent at each location at each trial
	
	if (is.null(cell.IDs)) cell.IDs <- sort(unique(spt[,2]))
	N.cells <- length(cell.IDs)
	x.bins <- seq(MazeRange[1], MazeRange[2], by=dx)
	x.bin.centres <- round(x.bins[-1] - dx/2,4)
	N.xbins <- length(x.bins) -1
	N.runs <- nrow(i.runs)
	
	dt.pos <- mean(diff(pos[,1]))

	## cellIDs are encoded as recorded from the tetrodes. 
	## We need to put them in a matrix instead where row number indexes the cell
	# cellnames <- seq(1:N.cells); 
	# names(cellnames) <- cell.IDs

	## the output matrix
	spike.map <- array(0, dim=c(N.cells+1, N.xbins, N.runs), dimnames=list(c(cell.IDs, 'time (s)'), paste(x.bin.centres, 'm'), paste('run', 1:N.runs)))
	
	## we only analyse runs
	for (i.run in 1:N.runs){
		isp.run.start <- min(which(spt[,1] > pos[i.runs[i.run,1],1])) # first spike of the run
		isp.run.end <- max(which(spt[,1] < pos[i.runs[i.run,2],1])) # last spike of the run
		spt.run.all <- spt[isp.run.start: isp.run.end,] # run spikes collected

		## tiem spent in each bin
		xpos.i <- pos[i.runs[i.run,1]:i.runs[i.run,2],2]
		xpos.i[xpos.i < MazeRange[1]] <- MazeRange[1]
		xpos.i[xpos.i > MazeRange[2]] <- MazeRange[2]
		spike.map[N.cells+1, ,i.run] <- (hist(xpos.i, br=x.bins, plot=F)$counts * dt.pos)

		i.spt.IDs <- spt.run.all[,2] %in% cell.IDs
		spt.run <- spt.run.all[i.spt.IDs,] # only cells that are given in the cell.ID
		
		## collect the position of the spikes
		x.spt.run <- get.x.spikes(spt.run[,1], pos)

		## match the cell names with the corresponding matrix index
		name.cells <- spt.run[,2] # orignial names 
		index.cells <- rep(NA, length(spt.run[,2])) # indexes
		for (i.cell in 1:N.cells) {
			isp.cell <- which(name.cells == cell.IDs[i.cell])
			if (length(isp.cell) > 0) index.cells[isp.cell] <- i.cell
		}
				
		N.sp.run <- nrow(spt.run)
		for (i.sp in 1:N.sp.run){
			ind.cell <- index.cells[i.sp]

			x.sp <- x.spt.run[i.sp]
			ind.x.sp <- ((x.sp - MazeRange[1]) %/% dx) + 1
			if (ind.x.sp < 1) ind.x.sp <- 1
			if (ind.x.sp > N.xbins) ind.x.sp <- N.xbins

			spike.map[ind.cell, ind.x.sp, i.run] <- spike.map[ind.cell, ind.x.sp, i.run] + 1
		}
		
		
	}
	
	spike.map
}



get.x.spikes <- function(t.spt, pos){
	## we assume that t.spt contains order time points
	## pos contains ordered time and position data
	## go along t and find the matching position
	
	x.spt <- rep(NA, length(t.spt))

	ii.first <- which.min((t.spt[1] - pos[,1])^2)
	x.spt[1] <- pos[ii.first,2]

	L <- length(t.spt)	
	ii.last <- which.min((t.spt[L] - pos[,1])^2)
	x.spt[L] <- pos[ii.last,2]
	
	pos.sp <- pos[ii.first:ii.last,]
	for (i in 2:L){
		ii <- which.min((t.spt[i] - pos.sp[,1])^2)
		x.spt[i] <- pos.sp[ii,2]
	}
	
	x.spt	
}

###############################################
pop.act.t <- function(spt, pos, i.runs, dt=0.1, cell.IDs=NULL, verbose=T){
	# spt: matrix of two columns: time, cell.ID
	# pos: matrix of two columns: time, x position
	# i.runs: row index of the start and end of run trials in the pos matrix
	# dt: approximate temporal resolution for the population activity - the true resolution will be an integer multiple of the position data, for simplicity...
	# cell.IDs: cells to be included in the spike.map

	# output: an matrix with the spikes of each cell on each time frame
	# dimensions: cell x time frames
	
	if (is.null(cell.IDs)) cell.IDs <- sort(unique(spt[,2]))
	N.cells <- length(cell.IDs)
	N.runs <- nrow(i.runs)
	
	dt.pos <- mean(diff(pos[,1]))
	nbins.pos <- ceiling(dt / dt.pos)
	
	## we only analyse runs
	for (i.run in 1:N.runs){
		if (verbose) cat ('calculating trial ', i.run, ' ')
			
		i.first <- rat$iruns.up[i.run,1]
		i.last <- rat$iruns.up[i.run,2]
		nframes.pos.all <- (i.last - i.first) %/% nbins.pos
		
		popact.run <- matrix(0, N.cells, nframes.pos.all)
		pos.run <- rep(NA, nframes.pos.all)
		
		for (i.frame in 1:nframes.pos.all){
			ipos1 <- i.first + (i.frame-1) * 4
			ipos2 <- i.first + i.frame * 4
			pos.run[i.frame] <- mean(pos[ipos1:(ipos2-1),2])
	
			t1 <- pos[ipos1,1]
			t2 <- pos[ipos2+1,1]
			
			isp.run.start <- min(which(spt[,1] > t1)) # first spike of the bin
			isp.run.end <- max(which(spt[,1] < t2)) # last spike of the bin
			
			spt.frame.all <- spt[isp.run.start: isp.run.end,2]
			i.spt.IDs <- spt.frame.all %in% cell.IDs
			spt.frame <- spt.frame.all[i.spt.IDs] # only cells that are given in the cell.ID
			if (length(spt.frame > 0)){
				for (sp in spt.frame) 	popact.run[which(cell.IDs == sp), i.frame] <- popact.run[which(cell.IDs == sp), i.frame] + 1
			}		
		}
		
		if (i.run==1) {
			popact <- popact.run 
			pos.all <- pos.run
		} else {
			popact <- cbind(popact, popact.run)
			pos.all <- c(pos.all, pos.run)
		}
			
		if (verbose) cat (sum(popact.run), ' spikes found in', nframes.pos.all, ' frames \n')

	}
	
	data <- list(pos=pos.all, popact=popact)
	data
}


skaggs93.info <- function(r.x, Px){
	## I = \sum_x r(x)p(x) log(r(x)/r)
	## r.x: estimated firing rate as function of x (spikes/sec)
	## Px: P(x), unnormalised
	## output: information rate, bit/s
	## Px <- total.t; r.x <- rate.x[,43]
	Px <- Px / sum(Px)
	r <- sum(Px * r.x)
	KL <- r.x * Px * log(r.x/r, 2)
	KL[r.x==0] <- 0
	I <- sum(KL)
	I
}




