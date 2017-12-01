
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
			plot(spt, rep(irep, length(spt)), pch=16, cex=0.3, ylim=c(0, nrep), xlim=c(min(t.sp)-dt, max(t.sp)), axes=F, xlab='', ylab='')
			axis(1)
			axis(2, las=2)
			mtext('time', 1, 2)
			mtext('repetition', 2, 2)
		} else {
			points(spt, rep(irep, length(spt)), pch=16, cex=0.3)			
		}
	}
}
