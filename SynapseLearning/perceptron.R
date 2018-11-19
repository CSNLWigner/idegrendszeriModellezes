
# sigmoid nonlinearity
sigm <- function(x, a=1, b=1, c=0, d=0){
	# a / (1 + exp(-b(x-c)))
	#a - amplitude: limit in inf; b - steepness; c - inflection
	# d: the derivative of the sigmoid function at x. d=0: the function; d=1: first derivative; d=2: second derivative
	bc <- exp(-b*(x-c))
	if (d>2) return ("d must be smaller than 2!")
	if (d==0) {
		sx <- a / (1+bc)
		} else {
			if (d==1) {
				sx <- (a*b*bc)/(1+bc)^2
			} else {
				sx <- (a*b^2*bc*(bc-1))/(1+bc)^3
			}
		}
	sx
}

## colormap from blue to red
blueRed2 <- function(x){
	# x has to be between 0 and 1
	if (min(x) < 0){ 
		x <- x - min(x)
		x <- x / max(x)
	}
	if (max(x) > 1){ 
		x <- x - min(x)
		x <- x / max(x)
	}
	
	cols <- rep(NA, length(x))
	for (i in 1:length(x)){
		if (x[i]<0.5) cols[i] <- rgb(1, 2*x[i], 2*x[i]) else cols[i] <- rgb(2-2*x[i],2-2*x[i],1)
	}
	cols
}


## one layer perceptron
train.perceptron <- function(w.init, patterns, target, Tmax=10, learn_rate=10, graphics=T){
	n.patterns <- nrow(patterns)
	n.inputs <- ncol(patterns)
	if (graphics){
		par(mfcol=c(3,4))
		par(mar=c(1,1,1,1))
	}
	w <- w.init
	m <- ceiling(max(abs(patterns)))
  miss <- rep(0, Tmax)
  
	for (time in 1:Tmax){
		output <- rep(0, n.patterns)
		dw <- rep(0, n.inputs)
		for (i in 1:n.patterns){
		    u <- patterns[i,]
		    v.u <- sigm(u %*% w)
			output[i] <- v.u
		    dw <- dw + as.vector(v.u * (1-v.u) * (v.u - target[i])) * u # gradient learning rule
		}

		v.u <- sigm(patterns %*% w)
		vv <- round(sigm(v.u, 1, 1000, 0.5))
		miss[time] <- sum(xor(vv, target))
		# cat('misclassification rate: ', round(miss[time]/n.patterns*100), '% \n')
		
		if (graphics){
			if (time < 12)	{
				plot(patterns[,1], patterns[,2], pch=16, col=blueRed2(output), xlab="", ylab="", cex=0.7, axes=F, xlim=c(-m, m), ylim=c(-m, m))
				sep2 <- (-1) * (w[1] * x  + w[3]) / w[2]
				lines(x, sep2, col=3, lwd=2)
			}
			if (time == Tmax)	{
				plot(patterns[,1], patterns[,2], pch=16, col=blueRed2(output), xlab="", ylab="", cex=0.7, main='final', axes=F, xlim=c(-m, m), ylim=c(-m, m))
				# axis(1); axis(2, las=2)
				sep2 <- (-1) * (w[1] * x  + w[3]) / w[2]
				lines(x, sep2, col=3, lwd=2)
			}
		}
		w <- w - learn_rate/2 * dw / n.patterns # gradient learning rule
	}
	
	w.final <- list(w=w, miss=miss)
	w.final
	
}


## two layer perceptron
train.perceptron.2L <- function(w.in, w.out, patterns, target, Tmax=50, e=10, graphics=T, verbose=T){
  # w.in  - initial input weights
  # w.out - initial output weights
  # patterns - patterns to classify (rows: different points, columns: dimensionality)
  # target: the target classification
  # Tmax: number of steps for training
  # e: learning rate
  # graphics: T or F, whether or not to plot the results
  # verbose: whether the error   should be reported during training
  
  n.patterns <- nrow(patterns)
	n.inputs <- ncol(patterns)
	n.hidden <- length(w.out)-1
	if (graphics){
		par(mfcol=c(3,4))
		par(mar=c(1,1,1,1))
		nshow <- ceiling(Tmax/11)+1
	}

	m <- ceiling(max(abs(patterns)))
	
	for (time in 1:Tmax){
		output <- rep(0, n.patterns)
		dw.in <- matrix(0, n.inputs, n.hidden)
		dw.out <- rep(0, n.hidden+1)
		error <- 0
		for (i in 1:n.patterns){
		 	u <- patterns[i,]
		  	y.u <- c(sigm(u %*% w.in), 1)
			v.u <- sigm(y.u %*% w.out)
			output[i] <- v.u
			dE.dw.out <- v.u * (1-v.u) * (v.u - target[i])
			error <- error +  (v.u - target[i])^2
		  	dw.out <- dw.out + as.vector(dE.dw.out) * y.u # gradient learning rule
			dy.dw.in <- outer(u, y.u[1:(n.hidden)] * (1-y.u[1:(n.hidden)]))
		 	dw.in <- dw.in +  t(as.vector(as.numeric(dE.dw.out) * w.out[1:(n.hidden)]) * t(dy.dw.in)) # gradient learning rule
		}
		
		if (graphics){
			if ((time %% nshow) == 1)	{
				title <- paste('t=', time)
				plot(patterns[,1], patterns[,2], pch=16, col=blueRed2(output), xlab="", ylab="", main=title, cex=1.7, axes=F, xlim=c(-m, m), ylim=c(-m, m))
			}
			if (time == Tmax)	{
				plot(patterns[,1], patterns[,2], pch=16, col=blueRed2(output), xlab="", ylab="", cex=1.7, main='final', axes=F, xlim=c(-m, m), ylim=c(-m, m))
			}
		}
		
		if (verbose) cat('error', error, '\n')
		w.in <- w.in - e/2 * dw.in / n.patterns # gradient learning rule
		w.out <- w.out - e/2 * dw.out / n.patterns # gradient learning rule
	}
	
    y.u <- cbind(sigm(patterns %*% w.in), rep(1, n.patterns))
	v.u <- sigm(y.u %*% w.out)
	vv <- round(sigm(v.u, 1, 1000, 0.5))
	miss <- sum(xor(vv, target))
	# cat('misclassification rate: ', round(miss/n.patterns*100), '% \n')
	
	w.final <- list(w.in=w.in, w.out=w.out, miss=miss)
	w.final
}


## two layer perceptron
eval.perceptron.2L <- function(w.in, w.out, patterns, graphics=T){
  n.patterns <- nrow(patterns)
  n.inputs <- ncol(patterns)
  n.hidden <- length(w.out)-1
  
  m <- ceiling(max(abs(patterns)))

  y.u <- cbind(sigm(patterns %*% w.in), rep(1, n.patterns))
  v.u <- sigm(y.u %*% w.out)
  vv <- round(sigm(v.u, 1, 1000, 0.5))

  if (graphics){
    # par(mfcol=c(1,1))
    par(mar=c(4,4,1,1))
    plot(patterns[,1], patterns[,2], pch=16, col=blueRed2(v.u), xlab="", ylab="", cex=0.7, axes=F, xlim=c(-m, m), ylim=c(-m, m))
  }
  v.u
}
