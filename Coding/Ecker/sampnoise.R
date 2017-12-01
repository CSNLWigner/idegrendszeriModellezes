x <- seq(-2, 2, length=100)
nSamps <- c(5,20,50,100,200)
trueVar = 1;
hx <- dnorm(x,0,sqrt(trueVar) / sqrt(nSamps[5]))
hgx <- dgamma(x,(nSamps[5] - 1) / 2, nSamps[5] / (2 * trueVar))

degf <- c(1, 3, 8, 30)
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c("5", "20", "50", "100", "200")

par(mfrow=c(1,2))

plot(x, hx, type="l", lty=2, xlab="measured value",
  ylab="Density", main="Mean, true value 0")

for (i in 1:4){
  lines(x, dnorm(x,0,sqrt(trueVar) / sqrt(nSamps[i])), lwd=2, col=colors[i])
}

plot(x, hgx, type="l", lty=2, xlab="measured value",
  ylab="Density", main="Variance, true value 1")

for (i in 1:4){
  lines(x, dgamma(x,(nSamps[i] - 1) / 2, nSamps[i] / (2 * trueVar)), lwd=2, col=colors[i])
}

legend("topleft", inset=.05, title="Number of samples",
  labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)