########################################################################
## read data from the matlab file
###########################################################################
## spike-file; clusters are made from an entire day (called a 2-3 h recording session in Karlsson 2009)
require('R.matlab')
require('colormap')

static <- readMat('data_v1_binned_static.mat')
moving <- readMat('data_v1_binned_static.mat')

i.ses <- 2
summary(static$data[[1]][[1]])


for (i.ses in 1:17){
	ses <- static$data[[i.ses]][[1]]
	names(ses) <- c('date', 'subject', 'conditions', 'contamination', 'tetrode', 'spikes', 'times')
	# ses <- moving$data[[i.ses]][[1]]
	# print(ses[[2]])
	# print(dim(ses[[3]]))

	condition <- ses[[3]][,1,]
	contamin <- ses[[4]]
	tet <- ses[[5]]
	spikes <- ses[[6]] # cells x conditions x time x repetitions
	times <- as.vector(ses[[7]])
	
	
	oris <- unlist(condition['orientation',])
	n.ori <- length(unique(oris))
	scaled.ori <- oris/180
	
	contrasts <- unlist(condition['contrast',])
	range.contrasts <- range(contrasts)
	scaled.contrasts <- (contrasts - min(contrasts)) / (max(contrasts)-min(contrasts)) * 2 + 1/2
	
	n.cells <- dim(spikes)[1]
	n.conditions <- dim(spikes)[2]
	L <- dim(spikes)[3]
	
	
	par(mfcol=c(5,6)); par(mar=c(2,2,1,1))
	for (i.cell in 1:n.cells){
		rates <- matrix(0, L, n.conditions)
		for (i.condition in 1:n.conditions){
			sp <- spikes[i.cell,i.condition,,]
			rates[,i.condition] <- apply(sp, 1, mean)*100 # Hz
		}
		matplot(times, rates, t="l", col=colormap(colormaps$jet,180)[oris], lty=1, lwd=scaled.contrasts, axes=F, ylim=c(0, 200))
		# readline(i.cell)
	}
	plot(apply(spikes*100, 1, mean), pch=16)
	readline(i.ses)
}


### selecting session 2, active cells with low conatmination

data <- static$data[[i.ses]][[1]]
names(data) <- c('date', 'subject', 'conditions', 'contamination', 'tetrode', 'spikes', 'times')
		
condition <- data[[3]][,1,]
contamin <- as.vector(data[[4]])
tet <- as.vector(data[[5]])
spikes <- data[[6]] # cells x conditions x time x repetitions
times <- as.vector(data[[7]])

rates <- apply(spikes*100, 1, mean)
keep.cells <- as.vector((rates>1) & (contamin<0.05))
plot(rates, col=keep.cells+1)

data$spikes <- data$spikes[keep.cells,,,]
data$contamination <- as.vector(data$contamination[keep.cells])
data$tetrode <- as.vector(data$tetrode[keep.cells])
data$conditions <- data$conditions[,1,]
data$times <- as.vector(data$times)
save(data, file='data_v1_binned_static_ses2.RData')


### plotting the rates 

n.cells <- dim(data$spikes)[1]
n.conditions <- dim(data$spikes)[2]
L <- dim(data$spikes)[3]


par(mfcol=c(3,3)); par(mar=c(2,2,1,1))
for (i.cell in 1:n.cells){
	rates <- matrix(0, L, n.conditions)
	for (i.condition in 1:n.conditions){
		sp <- data$spikes[i.cell,i.condition,,]
		rates[,i.condition] <- apply(sp, 1, mean)*100 # Hz
	}
	matplot(times, rates, t="l", col=colormap(colormaps$jet,180)[oris], lty=1, lwd=scaled.contrasts, axes=F, ylim=c(0, 200))
	abline(h=0); axis(2); axis(1)
	# readline(i.cell)
}

cell <- data$spikes[7,,,]
apply(cell, 1, mean)*100

par(mar=c(3,3,1,1))
plot.raster(cell[16,,])

plot.raster <- function(sp, t.sp=NULL){
	L <- nrow(sp)
	if (is.null(t.sp)) t.sp <- 1:L
	dt <- t.sp[2]-t.sp[1]	
	nrep <- ncol(sp)
	for (irep in 1:nrep){
		sp.rep <- sp[,irep]
		nsp <- unique(sp.rep)[unique(sp.rep)>0]
		spt <- 0
		for (ii in nsp) spt <- c(spt, rep(t.sp[which(sp.rep==ii)], each=ii))
		spt <- spt[-1]
		spt <- spt + runif(sum(sp.rep), -dt, 0)
		spt <- sort(spt)
		if (irep==1){
			plot(spt, rep(irep, length(spt)), pch=16, cex=0.3, ylim=c(0, nrep), xlim=c(0, max(t.sp)), axes=F, xlab='', ylab='')
			axis(1)
			axis(2, las=2)
			mtext('time', 1, 2)
			mtext('repetition', 2, 2)
		} else {
			points(spt, rep(irep, length(spt)), pch=16, cex=0.3)			
		}
	}
}
