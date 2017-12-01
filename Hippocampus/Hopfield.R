# Train a recurrent network using covariance rule
# Optimal recall dynamics is from Savin, Dayan and Lengyel, 2013
library(gplots)

f <- 0.15 # the sparseness of the representation
sig2.w <- f^2 - 2*f^3 + f^4 # the varinace of the weights after storing a single pattern
N <- 1024 # number of neurons
# plot(seq(0, N), dbinom(seq(0,N), N, f))

M <- 50 # number of patterns stored try with 100, 200 and 300
X <- matrix(rbinom(N*M, 1, f), M, N) # the patterns to be stored


library(png)
ims <- array(NA, dim=c(7, 32, 32))
imnames <- list.files('Figs32/')
for (i.image in 1:7){
	imname <- paste('Figs32/', imnames[i.image], sep='')
	im1 <- readPNG(imname)
	for (i in 1:32){
		for (j in 1:32){
			ims[i.image, i,j] <- round(1-im1[33-j, i,1])
		}
	}
	image(ims[i.image,,], col=grey(seq(0,1, by=0.1)))
	X[i.image,] <- as.vector(ims[i.image,,])
}

apply(X, 1, sum)[1:7]
mean(apply(X, 1, sum)[1:6])
# W <- matrix(0, N, N) # the weight matrix
# for (m in 1:M){
	# W <- W + outer(X[m,] - f, X[m,]-f) # additive, covariance rule
# }

W <- t(X-f) %*% (X-f)
diag(W) <- 0


### We assume a noise model parametrised by the probability of switching on: P(x~=1 | x=0) = r < f
### the probability of switching off  P(x~=0 | x=1) = (1-f)*r/f
### this way the noise does NOT change the sparseness of the input.
### if r = f/2 the noise is strong, but recall is still possible.
x <- X[1,]
r <- f/2; rho <- (1-f) * r / f
px <- rep(rho, N); px[!x] <- r 
cx <- rbinom(N, 1, px)
nx <- xor(x, cx)

image(matrix(x, 32), col=grey(seq(0,1, by=0.1)))
image(matrix(nx, 32), col=grey(seq(0,1, by=0.1)))

plot (X %*% x)
points (X %*% nx, col=2, pch=21, bg=2)

###########################
## check that patterns are stable - input is transient


stats.recall <- rep(NA, 9)
kk <- 0
mems <- matrix(NA, 4, 50)
rownames(mems) <- c("recalled vs target", "recalled vs best", "cue vs. target", "i.recalled")
quartz()

for (m in 1:7){ # for each pattern
	x <- X[m,] 		# we start from the stored pattern
	px <- rep(rho, N); px[!x] <- r 
	cx <- rbinom(N, 1, px)
	x.t <- xor(x, cx)
	nx <- x.t
	converged <- F
	k <- 0

	while (converged == F){
		dev.set(2)
		par(mfcol=c(1,1))
		par(mar=c(1,1,1,1))
		I <- 1/(sig2.w * (M-1)) * (W %*% x.t - (1-2*f)^2/2*sum(x.t) - f * rowSums(W) - f^2 * (1-2*f) / 2) + log(f/(1-f)) # optimal dynamics from Eq. 6. of Savin et al., 2013
		x.I <- rep(F, N)
		x.I [I>0] <- T
		change.x <- which(xor(x.t, x.I))
		if (length(change.x) == 0){ # convergence
			converged <- T
		} else {
			ch.x <- sample(change.x, 1)
			x.t[ch.x] <- 1 - x.t[ch.x]
		}
		# if ((k %% 20) == 1) image(matrix(x.t, 32), col=grey(seq(0,1, by=0.1)), axes=F)	
		k <- k+1
		if (k>1000) break
	}
	
	if (m < 8){
		dev.set(3)
		par(mfcol=c(1,3))
		par(mar=c(1,1,1,1))
		image(matrix(x, 32), col=grey(seq(0,1, by=0.1)), axes=F)
		image(matrix(nx, 32), col=grey(seq(0,1, by=0.1)), axes=F)
		image(matrix(x.t, 32), col=grey(seq(0,1, by=0.1)), axes=F)	
	}
	
	# plot (X %*% x)
	# points (X %*% nx, col=2, pch=21, bg=2, cex=0.7)
	mem <- X %*% x.t
	# points (mem, col=3, t="h")
	i.recalled <- which.max(mem)
	x.recalled <- X[i.recalled,]
	mems[1,m] <- 2 * x %*% x.t / (x.t %*% x.t + x %*% x)
	mems[2,m] <- 2 * x.recalled %*% x.t / (x.t %*% x.t + x.recalled %*% x.recalled)
	mems[3,m] <- 2 * nx %*% x / (x %*% x + nx %*% nx)
	mems[4,m] <- i.recalled

	cat("", m)
}
	
# stats.recall[1:2]  <- c(mean(mems[1,]), sd(mems[1,]))	
# stats.recall[3:4]  <- c(mean(mems[2,]), sd(mems[2,]))	
# stats.recall[5:6]  <- c(mean(mems[3,]), sd(mems[3,]))	
# stats.recall[7]  <- sum((!mems[4,]-seq(1,100))) 		
# stats.recall[8]  <- sum(mems[1,] == 1) 		
# stats.recall[9]  <- sum(mems[1,] > 0.8) 		

# stats.recall <- stats.recall[1:6,]

# plot(Ms, stats.recall[,1], ylim=c(0, 1), log="x", xlab="number of stored patterns", ylab="quality of recall", axes=F, t="n"); axis(1); axis(2, las=2)
# polygon(c(Ms, rev(Ms)), c(stats.recall[,1] + stats.recall[,2], rev(stats.recall[,1] - stats.recall[,2])), col=2, lty=2)
# lines(Ms, stats.recall[,1], lwd=2)
# polygon(c(Ms, rev(Ms)), c(stats.recall[,5] + stats.recall[,6], rev(stats.recall[,5] - stats.recall[,6])), col=3, lty=2)
# lines(Ms, stats.recall[,5], lwd=2)

# lines(Ms, stats.recall[,9]/100, col=4, lwd=1, pch=21, t="o", bg=4)
# legend(100, 0.5, bty="n", legend=c("cue", "recalled", "percent of > 0.8"), lty=1, col=c(3, 2, 4), lwd=c(3,3,1))

# abline(0,1)

