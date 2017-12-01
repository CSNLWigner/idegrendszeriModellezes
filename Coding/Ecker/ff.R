load(file="ecker_sess8_spikes.rda")
nBins <- dim(spikes)[3]
nCell <- dim(spikes)[1]
nConditions <- dim(spikes)[2]
binIdx <- seq(1,nBins,1)

evokedStartBin <- 30
evokedEndBin <- 70

lcIdx <- seq(1, nConditions, 2)
hcIdx <- lcIdx + 1
lcSpikes <- spikes[,lcIdx,,]
hcSpikes <- spikes[,hcIdx,,]

lcVars <- apply(lcSpikes, c(1,2,3), var)
lcMeans <- apply(lcSpikes, c(1,2,3), mean)
lcFF <- lcVars / lcMeans
lcFF[is.nan(lcFF)] <- 0 
lcMeanFF <- apply(lcFF, c(3), mean)
lcMeanRate <- apply(lcSpikes,c(3),mean)

hcVars <- apply(hcSpikes, c(1,2,3), var)
hcMeans <- apply(hcSpikes, c(1,2,3), mean)
hcFF <- hcVars / hcMeans
hcFF[is.nan(hcFF)] <- 0 
hcMeanFF <- apply(hcFF, c(3), mean)
hcMeanRate <- apply(hcSpikes,c(3),mean)

par(mfrow=c(1,2))
plot(binIdx, hcMeanRate, type="l", xlab="time bin", ylab="Firing rate")
lines(binIdx, lcMeanRate, lwd=2, col="grey")
plot(binIdx, hcMeanFF, type="l", xlab="time bin", ylab="FF")
lines(binIdx, lcMeanFF, lwd=2, col="grey")
