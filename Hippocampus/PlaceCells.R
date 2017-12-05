########################################################################
## read data from the R datafile, generated from the matlab file

source("PlaceCellFunctions.R")
load('./Achilles.RData')
summary(rat)
require(viridis)
pos <- rat$pos
spt <- rat$spt


dx <- 0.05 # cm, resolution 
x.breaks <- seq(rat$MazeRange[1], rat$MazeRange[2], by=dx)
x.mids <- round(x.breaks[-1] - dx/2, 3)

act.runs <- cell.maps.runs(spt, pos, i.runs=rat$iruns.up, dx=0.05, MazeRange=rat$MazeRange, cell.IDs=rat$PyrIDs)



###############################################
## divide spike count with the occupancy time to get firing rates
ratemaps.t <- apply(act.runs[1:120,,], c(1,2), sum)
Tmap <- apply(act.runs[121,,], 1, sum)

ratemaps.all <- t(ratemaps.t) / Tmap
matplot(x.mids, ratemaps.all, t='l', lty=1, col=rainbow(120), xlab='x position (m)', ylab='firing rate (Hz)')

##############################################################
## firing rates estimated on individual runs - still need to divide with occupancy time!

image(x.mids, 1:42, act.runs[11,,] / act.runs[121,,], col=viridis(24), xlab='x position (m)', ylab='trial')

## index of active cells
## active: peak firing rate is larger than 5 Hz

i.cells.active <- which(apply(ratemaps.all, 2, max) > 5) 
par(mfcol=c(7,8)); par(mar=c(1,1,3,1))
for (i.cell in i.cells.active){
	image(x.mids, 1:42, act.runs[i.cell,,] / act.runs[121,,], col=viridis(24), xlab='', ylab='', main='', ylim=c(0, 42), axes=F)
	lines(x.mids, ratemaps.all[,i.cell], col=viridis(3)[3], lwd=2)
	info.cell <- skaggs93.info(ratemaps.all[,i.cell], Tmap) # bits/sec
	title(main=paste('info:', round(info.cell, 3), 'bit / s'))
	# readline(i.cell)
}


par(mfcol=c(1,2)); par(mar=c(4,4,4,2))
for (i.cell in c(70, 97)){
	image(x.mids, 1:42, act.runs[i.cell,,] / act.runs[121,,], col=viridis(24), xlab='x position (m)', ylab='trial', main='', ylim=c(0, 42))
	lines(x.mids, ratemaps.all[,i.cell], col=viridis(3)[3], lwd=2)
	info.cell <- skaggs93.info(ratemaps.all[,i.cell], Tmap) # bits/sec
	title(main=paste('info:', round(info.cell, 3), 'bit / s'))
	# readline(i.cell)
}

##############################################################
## ratemaps - sorted according to the position of the peak

ratemaps <- ratemaps.all[,i.cells.active]
ii.maxs <- apply(ratemaps, 2, which.max)
sort.peaks <- sort(ii.maxs, ind=T)$ix
par(mfcol=c(1,2))
matplot(ratemaps[,sort.peaks], t='l', lty=1, col=rainbow(60), xlab='x position (m)', ylab='firing rate (Hz)')
image(ratemaps[,sort.peaks], col=viridis(24), xlab='x position (m)', ylab='trial')



##############################################################
## we prepare population activity for decoding

all.data <- pop.act.t(spt, pos, i.runs=rat$iruns.up, dt=0.1, cell.IDs=rat$PyrIDs[i.cells.active])

L.all <- length(all.data$pos)
L.train <- 2000
L.test <- L.all - L.train

i.train <- sample (1:L.all, L.train)
i.test <- which(!(1:L.all %in% i.train))

data.train <- list(pos=all.data$pos[i.train], popact=all.data$popact[,i.train])
data.test <- list(pos=all.data$pos[i.test], popact=all.data$popact[,i.test])

##############################################################
## we will use a perceptron for decoding - alternative would be a Bayesian decoding...

patterns.train <- cbind(t(data.train$popact), rep(1, L.train))
patterns.test <- cbind(t(data.test$popact), rep(1, L.test))
target.train <-data.train$pos / 4 + 0.35
target.test <-data.test$pos / 4 + 0.35



n.obs <- ncol(patterns.train)
n.hidden <- 16 # number of hidden neurons
set.seed(114)

w.in <- matrix(rnorm(n.obs*n.hidden), n.obs)
w.out <- matrix(rnorm(n.hidden+1), n.hidden+1)

source('../SynapseLearning/perceptron.R', chdir = TRUE)
ww <- train.perceptron.2L(w.in, w.out, patterns.train, target.train, Tmax=200, graphics=F)
ww <- train.perceptron.2L(ww$w.in, ww$w.out, patterns.train, target.train, Tmax=200, graphics=F, e=10)

pred.pos.train <- (eval.perceptron.2L(ww$w.in, ww$w.out, patterns.train, graphics=F)-0.35)*4
pred.pos.test <- (eval.perceptron.2L(ww$w.in, ww$w.out, patterns.test, graphics=F)-0.35)*4


err.train <- sqrt(mean((data.train$pos-pred.pos.train)^2))
plot(data.train$pos, pred.pos.train, xlab='true position', ylab='decoded position', main=paste('training error:', round(err.train, 2), 'cm'))
abline(0, 1, col=2)

err.test <- sqrt(mean((data.test$pos-pred.pos.test)^2))
plot(data.test$pos, pred.pos.test, xlab='true position', ylab='decoded position', main=paste('test error:', round(err.test, 2), 'cm'))
abline(0, 1, col=2)
