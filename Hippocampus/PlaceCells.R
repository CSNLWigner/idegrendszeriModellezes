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

act.runs.left <- cell.maps.runs(spt, pos, i.runs=rat$iruns.up, dx=0.05, MazeRange=rat$MazeRange, cell.IDs=rat$PyrIDs)
act.runs.right <- cell.maps.runs(spt, pos, i.runs=rat$iruns.down, dx=0.05, MazeRange=rat$MazeRange, cell.IDs=rat$PyrIDs)



###############################################
## divide spike count with the occupancy time to get firing rates
ratemaps.left.t <- apply(act.runs.left[1:120,,], c(1,2), sum)
Tmap.left <- apply(act.runs.left[121,,], 1, sum)

ratemaps.right.t <- apply(act.runs.right[1:120,,], c(1,2), sum)
Tmap.right <- apply(act.runs.right[121,,], 1, sum)

ratemaps.left.all <- t(ratemaps.left.t) / Tmap.left
ratemaps.right.all <- t(ratemaps.right.t) / Tmap.right
matplot(x.mids, ratemaps.left.all, t='l', lty=1, col=rainbow(120), xlab='x position (m)', ylab='firing rate (Hz)')
matplot(x.mids, ratemaps.right.all, t='l', lty=1, col=rainbow(120), xlab='x position (m)', ylab='firing rate (Hz)')

##############################################################
## firing rates estimated on individual runs - still need to divide with occupancy time!

image(x.mids, 1:42, act.runs.left[11,,] / act.runs.left[121,,], col=viridis(24), xlab='x position (m)', ylab='trial')

## index of active cells
## active: peak firing rate is larger than 5 Hz

i.cells.active.left <- which(apply(ratemaps.left.all, 2, max) > 5) 
N.cells.active.left <- length(i.cells.active.left)

i.cells.active.right <- which(apply(ratemaps.right.all, 2, max) > 5) 
N.cells.active.right <- length(i.cells.active.right)


i.cells.active <- which((apply(ratemaps.left.all, 2, max) > 5) | (apply(ratemaps.right.all, 2, max) > 5)) 
N.cells.active <- length(i.cells.active)

par(mfcol=c(8,10)); par(mar=c(1,1,1,1))
for (i.cell in i.cells.active){
	image(x.mids, 1:42, act.runs.left[i.cell,,] / act.runs.left[121,,], col=viridis(24), xlab='', ylab='', main='', ylim=c(0, 42), axes=F, xlim=c(-0.3, 3))
	image(x.mids+1.65, 1:42, act.runs.right[i.cell,,] / act.runs.left[121,,], col=viridis(24, option='B'), xlab='', ylab='', main='', ylim=c(0, 42), axes=F, xlim=c(-0.3, 3), add=T)
	lines(x.mids, ratemaps.left.all[,i.cell], col=viridis(3)[3], lwd=2)
	lines(x.mids+1.65, ratemaps.right.all[,i.cell], col=viridis(3)[3], lwd=2)
	# info.cell <- skaggs93.info(ratemaps.all[,i.cell], Tmap) # bits/sec
	# title(main=paste('info:', round(info.cell, 3), 'bit / s'))
	# readline(i.cell)
}


par(mfcol=c(1,2)); par(mar=c(4,4,4,2))
for (i.cell in c(70, 97)){
	image(x.mids, 1:42, act.runs.left[i.cell,,] / act.runs.left[121,,], col=viridis(24), xlab='x position (m)', ylab='trial', main='', ylim=c(0, 42))
	lines(x.mids, ratemaps.left.all[,i.cell], col=viridis(3)[3], lwd=2)
	info.cell <- skaggs93.info(ratemaps.left.all[,i.cell], Tmap.left) # bits/sec
	title(main=paste('info:', round(info.cell, 3), 'bit / s'))
	# readline(i.cell)
}

##############################################################
## ratemaps - sorted according to the position of the peak

ratemaps.left <- ratemaps.left.all[,i.cells.active]
ii.maxs.left <- apply(ratemaps.left, 2, which.max)
sort.peaks.left <- sort(ii.maxs.left, ind=T)$ix

ratemaps.right <- ratemaps.right.all[,i.cells.active]
ii.maxs.right <- apply(ratemaps.right, 2, which.max)
sort.peaks.right <- sort(ii.maxs.right, ind=T)$ix


par(mfcol=c(1,2))
matplot(x.mids, ratemaps.left[,sort.peaks.left], t='l', lty=1, col=rainbow(60), xlab='x position (m)', ylab='firing rate (Hz)')
image(x.mids, 1:N.cells.active, ratemaps.left[,sort.peaks.left], col=viridis(24), xlab='x position (m)', ylab='cell')

par(mfcol=c(2,2))
image(x.mids, 1:N.cells.active, ratemaps.left[,sort.peaks.left], col=viridis(24), xlab='x position (m)', ylab='cell', main="left, sorted: left")
image(x.mids, 1:N.cells.active, ratemaps.right[,sort.peaks.left], col=viridis(24), xlab='x position (m)', ylab='cell', main="right, sorted: left")
image(x.mids, 1:N.cells.active, ratemaps.left[,sort.peaks.right], col=viridis(24), xlab='x position (m)', ylab='cell', main="left, sorted: right")
image(x.mids, 1:N.cells.active, ratemaps.right[,sort.peaks.right], col=viridis(24), xlab='x position (m)', ylab='cell', main="right, sorted: right")


## correlations between the place fields
rate.l <- ratemaps.left
rate.r <- ratemaps.right
cor.rl <- rep(0, N.cells.active)
cor.rl.rand <- rep(0, N.cells.active)
for (i in 1:N.cells.active) cor.rl[i] <- cor(ratemaps.left[,i], ratemaps.right[,i])
for (i in 1:N.cells.active) cor.rl.rand[i] <- cor(ratemaps.left[,i], ratemaps.right[,sort.peaks.right[i]])

br <- seq(-1,1,by=0.05)
hist(cor.rl.rand, br=br, col=rgb(.75, 0, 0, alpha=1), main='place field correlations', xlab='correlation')
hist(cor.rl, br=br, col=rgb(0, 0.35, .75, alpha=0.5), add=T)
abline(v=mean(cor.rl.rand), col=rgb(.75, 0, 0, alpha=1), lwd=2)
abline(v=mean(cor.rl), col=rgb(0, 0.35, 0.75, alpha=1), lwd=2)
legend('topleft', legend=c('true', 'shuffled'), pch=22, pt.bg=c(rgb(0, 0.35, 0.75, alpha=.5), rgb(.75, 0, 0, alpha=1)))

##############################################################
## deocoding population activity using a perceptron
##############################################################

##############################################################
## we prepare population activity for decoding

all.data <- pop.act.t(spt, pos, i.runs=rat$iruns.left, dt=0.1, cell.IDs=rat$PyrIDs[i.cells.active])

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


##############################################################
## preplay and replay analysis
##############################################################

## collect all spikes in a big population activity matrix

i1 <- 3001 # start a bit before the first run
i2 <- 80000 # end after the last run
isp.1 <- min(which(spt[,1] > pos[i1,1]))
isp.2 <- max(which(spt[,1] < pos[i2,1]))
spt.pop.all <- spt[isp.1:isp.2,]

ii.act.pop <- spt.pop.all[,2] %in% rat$PyrIDs[i.cells.active.left]
spt.pop.IDs <- spt.pop.all[ii.act.pop,]

spt.pop <- spt.pop.IDs
for (i.cell in 1:N.cells.active.left){
	ii <- rat$PyrIDs[i.cells.active.left[i.cell]]
	i.sp <- which(spt.pop.IDs[,2] == ii)
	spt.pop[i.sp,2] <- i.cell
}


## time and poition for the population activity...
tpop <- pos[i1:i2, 1]
xpop <- pos[i1:i2, 2]
xpop <- xpop - min(xpop)

dt.pos <- mean(diff(tpop))
popact <- Matrix(0, N.cells.active, i2-i1+1)

## for each active cell we find its spikes and add to the population activity
for (i.cell in 1:N.cells.active.left){
	t.sp <- spt.pop[which(spt.pop[,2] == i.cell),1]
	i.sp <- (t.sp - pos[i1,1]) %/% dt.pos
	for (jj in i.sp) popact[i.cell,jj] <- popact[i.cell,jj] + 1
	cat('cell', i.cell, 'done \n')
}

sum(popact)


## total population activity - to detect candidate replay events
poprate <- colSums(popact)
sum(poprate)
poprate.f <- filter(poprate, rep(1, 4))

## speed of the animal 
speed <- abs(diff(xpop))
speed <- c(speed[1], speed)


## candidate spw is where there are at least 40 spikes in 0.1 s and the rat is not moving
ind.spw <- which((speed < 0.0005) & (poprate.f > 40))

## some spw-s are detected more than once
ind.different.spw <- rep(T, length(ind.spw))
for (ii in 2:length(ind.spw)){
	if ((ind.spw[ii] - ind.spw[ii-1]) < 10) ind.different.spw[ii] <- F
}
ind.spw <- ind.spw[ind.different.spw]
t.spw <- tpop[ind.spw]
x.spw <- xpop[ind.spw]

ratemaps.left <- ratemaps.left.all[,i.cells.active.left]
ii.maxs.left <- apply(ratemaps.left, 2, which.max)
sort.peaks.left <- sort(ii.maxs.left, ind=T)$ix

png('ReplayEvents.png', 2500, 1800, pointsize=36)
par(mfcol=c(7,8)); par(mar=c(1,1,1,1))
for (i.spw in 1:length(ind.spw)){
	t.start <- t.spw[i.spw] - 0.15
	t.end <- t.spw[i.spw] + 0.35

	isp.1 <- min(which(spt.pop[,1] > t.start))
	isp.2 <- max(which(spt.pop[,1] < t.end))
	spw <- spt.pop[isp.1:isp.2,]
	
	cells <- spw[,2]
	cells.left <- cells
	# cells.right <- cells
	cols.cells <- rep(0, length(cells))
	for (i.sp in 1:length(cells)) {
		cells.left[i.sp] <- which(sort.peaks.left == cells[i.sp])
		# cells.right[i.sp] <- which(sort.peaks.right == cells[i.sp])
		cols.cells[i.sp] <- rainbow(N.cells.active.left, end=0.7)[which(sort.peaks.left == cells[i.sp])]
	}
	
	
	# plot(spw[,1], cells.left, pch=16, col=4)
	# abline(v=t.spw[i.spw], col=grey(0.75))
	title <- paste('t:', round(t.spw[i.spw], 1), 'x:', round(x.spw[i.spw], 2))
	plot(spw[,1], cells.left, pch=16, col=cols.cells, axes=F, main=title)
	abline(v=t.spw[i.spw], col=grey(0.75))
	box()
	# readline(i.spw)	
}
dev.off()



## look at the different runs - individually
N.runs <- nrow(rat$iruns.up)

# for (i in 1:N.runs){
i <- 30
is <- matrix(NA, 3,2)
is[1,] <- c(19211.35, 19211.6)
is[2,] <- c(rat$pos[rat$iruns.up[i,1],1] - 5, rat$pos[rat$iruns.up[i,2],1] + 5)
is[3,] <- c(19222.3, 19222.55)

pdf(file='ReplayPreplay.pdf', 10, 2.5, useDingbats=F)
layout(matrix(c(1,2,3), 1), c(1,3,1))
par(mar=c(2,3,2,1))
source('~/Programs/hGLMs/Rlib/Graphics/scalebar.R')

for (j in 1:3){

	t.start <- is[j,1]
	t.end <- is[j,2]

	isp.1 <- min(which(spt.pop[,1] > t.start))
	isp.2 <- max(which(spt.pop[,1] < t.end))
	spw <- spt.pop[isp.1:isp.2,]
	
	cells <- spw[,2]
	cells.left <- cells
	# cells.right <- cells
	cols.cells <- rep(0, length(cells))
	for (i.sp in 1:length(cells)) {
		cells.left[i.sp] <- which(sort.peaks.left == cells[i.sp])
		# cells.right[i.sp] <- which(sort.peaks.right == cells[i.sp])
		cols.cells[i.sp] <- rainbow(N.cells.active.left, end=0.7)[which(sort.peaks.left == cells[i.sp])]
	}
	
	
	# plot(spw[,1], cells.left, pch=16, col=4)
	# abline(v=t.spw[i.spw], col=grey(0.75))
	# title <- paste('t:', round(t.spw[i.spw], 1), 'x:', round(x.spw[i.spw], 2))
	plot(spw[,1], cells.left, pch=16, cex=0.6, col=cols.cells, axes=F, main='', xlab='', ylab='', ylim=c(0, 55), xlim=c(t.start, t.end))
	abline(v=t.spw[i.spw], col=grey(0.75))
	

	lines(c(is[1,1],is[1,1],is[1,2],is[1,2],is[1,1]), c(0,55, 55, 0, 0))
	lines(c(is[3,1],is[3,1],is[3,2],is[3,2],is[3,1]), c(0,55, 55, 0, 0))
	if (j==2){
		i1 <- rat$iruns.up[i,1] - 5/dt.pos
		i2 <- rat$iruns.up[i,2] + 5/dt.pos
		lines(rat$pos[i1:i2,1], (rat$pos[i1:i2,2]+0.5)*25, t='l')
		scalebar2(1, 0, '1 s', pos='topleft')
		axis(2, las=2)
		mtext('cells', 2, 2, cex=0.7)
	}
}
dev.off()	
	# readline(i)	
	
# }
