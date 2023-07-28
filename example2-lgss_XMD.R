##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a LGSS model.
#
# Johan Dahlin <uni (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################
setwd("/home/daixiongming/Documents/Code/pmh-tutorial-2.1/r_chebyshev")

source("helpers/dataGeneration.R")
source("helpers/stateEstimation.R")
source("helpers/parameterEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

# Save plot to file
savePlotToFile <- FALSE
savePlotToFile <- TRUE

##############################################################################
# Define the model and generate data
# x[t + 1] = phi * x[t] + sigmav * v[t],    v[t] ~ N(0, 1)
# y[t] = x[t] + sigmae * e[t],              e[t] ~ N(0, 1)
##############################################################################
phi <- 0.75
sigmav <- 1.00
sigmae <- 0.10
T <- 250
initialState <- 0

data <- generateData(c(phi, sigmav, sigmae), T, initialState)



####################################################################################################
######################################################################################################
###################################################################################################
####################################################################################################
p=1
k=4*p
f=function(x)
  dnorm(x,0,1)

lf=function(x) #log funtion of normal distribution
  -x^2/2

energy1=function(x)
{
  d=abs(D-x)+10^(-10)
  #val=-2*log(f(x))+log(sum(1/(fev^2*d^k)))
  val=-2*lf(x)+log(sum(exp(-2*lfev-k*log(d))))
  return(val)
}
#U=1
#U=-1
#U=-0.5
#U=0.5
#U=1.5
#U=-1.5
U=2
#U=-2
#U=-8  if U=8 the result is very bad
#W=2 it is normal but not good

#W=1.2 it is good
#W=-1.2 it is bad
#W=1.5 it is good not better than 1.2
#W=1.1  it is good
#W=0.6 good

#W=0.4 not enough good

#W=0.2 not enough good
W=1.6

# This test works well for different U values and W values, the interval of the W is less than 2, we can 

energy2=function(x)
{
  d=abs(D-x)+10^(-10)
  #val=-2*log(f(x))+log(sum(1/(fev^2*d^k)))
  val=-2*W*lf(x)+log(sum(exp(-2*W*lfev+U*d-k*log(d)))) #This is very bad!
  return(val)
}


ac=1 # 1 is also good!
k=80
m=40
W=1
energy=function(x)
{
  d=abs(D-x)+10^(-10)
  #d_kernel=kernel_CF(d)
  #d=1/d
  d_kernel=d
  #val=-2*log(f(x))+log(sum(1/(fev^2*d^k)))
  val=-(k/(2*p))*lf(x)+log(sum(exp(-(k/(2*p))*lfev-k*log(d_kernel))))
  A=k*log(d_kernel)
  B=ac*d_kernel
  C=(-1)*lfev*lf(x)*W
  D=C+B
  val_w=log(sum(exp(-(m)*D-A)))  #This negative is very important to be loaded in!
  return(val_w)
}



energy_XMD2=function(x)
{
  d=abs(D-x)+10^(-10)
  d_kernel=kernel_CF(d)
  #val=-2*log(f(x))+log(sum(1/(fev^2*d^k)))
  val=-(k/(2*p))*lf(x)+log(sum(exp(-(k/(2*p))*lfev+k*log(d_kernel))))
  A=k*log(d_kernel)
  B=ac*d_kernel
  C=(-1)*lfev*lf(x)
  D=C+B
  val_w=log(sum(exp(-(m)*D+A)))
  return(val_w)
}



N=200
D=0
n=length(D)
plot(cbind(D,0), xlim=c(-3,3), type="n")
text(cbind(D,0), labels=1:n)
one=rep(1,n)
fev=f(D)
lfev=lf(D)
S=matrix(seq(-3,3,length=10*N),ncol=1)
S1=matrix(seq(-3,3,length=10*95),ncol=1)
S2=matrix(seq(-3,3,length=10*15),ncol=1)
S3=matrix(seq(-3,3,length=10*180),ncol=1)
S4=matrix(seq(-3,3,length=10*55),ncol=1)



for(j in 2:N)
{
  val=apply(S,1,energy)
  val1=apply(S1,1,energy)
  val2=apply(S2,1,energy)
  val3=apply(S3,1,energy)
  val4=apply(S4,1,energy)
  new=S[which.min(val)]
  new1=S1[which.min(val1)]
  
  new2=S2[which.min(val2)]
  
  new3=S3[which.min(val3)]
  
  new4=S4[which.min(val4)]
  
  new=max(new,new1,new2,new3,new4)
  points(new,runif(1,-.1,.1), type="n")
  text(new,runif(1,-.1,.1),labels=j, cex=.5)
  D=c(D,new)
  n=length(D)
  #fev=c(fev,f(new))
  lfev=c(lfev,lf(new))
}
hist(D,breaks=15, xlim=c(-3,3), prob=T, main="Normal (0,1)")
curve(dnorm(x,0,1), from=-3,to=3,col=2, add=T)

pdf(file="normal_paper.pdf")
hist(D,breaks=15, xlim=c(-3,3), prob=T, main="Normal (0,1)" ,xlab="x", ylab="density", cex.lab=1.25)
curve(dnorm(x,0,1), from=-3,to=3,col=2,lty=2, lwd=2, add=T)
dev.off()

qqnorm(D)

pdf(file="QQnorm.pdf")
plot(qnorm(((1:N)-.5)/N),sort(D),xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="Normal(0,1)")
abline(0,1)
dev.off()
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################






##############################################################################
# PMH
##############################################################################
#initialPhi <- 0.50
initialPhi <- 0.40
noParticles <- 100
noBurnInIterations <- 1000
#noIterations <- 5000
noIterations <- 2200


if (loadSavedWorkspace) {
  load("savedWorkspaces/example2-lgss.RData")
} else {
  res1 <- particleMetropolisHastingsD(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      #stepSize = 0.01,
      stepSize = 0.05,
      D
    )
  res2 <- particleMetropolisHastingsD(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      stepSize = 0.10,
      D
    )
  res3 <- particleMetropolisHastingsD(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      stepSize = 0.50,
      D
    )
}

##############################################################################
# Plot the results
##############################################################################
resTh1 <- res1[noBurnInIterations:noIterations,]
resTh2 <- res2[noBurnInIterations:noIterations,]
resTh3 <- res3[noBurnInIterations:noIterations,]

# Estimate the KDE of the marginal posteriors
kde1  <- density(resTh1,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)
kde2  <- density(resTh2,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)
kde3  <- density(resTh3,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("figures/example2-lgss.pdf",
            height = 6,  #10, 8 is the original
            width = 12)
}

layout(matrix(1:9, 3, 3, byrow = TRUE))
par   (mar = c(4, 5, 0, 0))

# Plot the parameter posterior estimate
hist(
  resTh1,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde1, lwd = 2, col = "#7570B3")
abline(v = mean(resTh1),
       lwd = 1,
       lty = "dotted")

hist(
  resTh2,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde2, lwd = 2, col = "#E7298A")
abline(v = mean(resTh2),
       lwd = 1,
       lty = "dotted")

hist(
  resTh3,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde3, lwd = 2, col = "#66A61E")
abline(v = mean(resTh3),
       lwd = 1,
       lty = "dotted")

# Plot the trace of the Markov chain during 1000 iterations after the burn-in
grid <- seq(noBurnInIterations, noBurnInIterations + 1000 - 1, 1)
#grid <- seq(noBurnInIterations, noBurnInIterations + 200 - 1, 1)


plot(
  grid,
  resTh1[1:1000],
  col = '#7570B3',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh1),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh1[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)

plot(
  grid,
  resTh2[1:1000],
  col = '#E7298A',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh2),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh2[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25)
)

plot(
  grid,
  resTh3[1:1000],
  col = '#66A61E',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh3),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh3[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)

# Plot the ACF of the Markov chain

res1ACF <- acf(resTh1, plot = FALSE, lag.max = 60)
plot(
  res1ACF$lag,
  res1ACF$acf,
  col = '#7570B3',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res1ACF$lag, rev(res1ACF$lag)),
  c(res1ACF$acf, rep(0, length(res1ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

res2ACF <- acf(resTh2, plot = FALSE, lag.max = 60)
plot(
  res2ACF$lag,
  res2ACF$acf,
  col = '#E7298A',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res2ACF$lag, rev(res2ACF$lag)),
  c(res2ACF$acf, rep(0, length(res2ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

res3ACF <- acf(resTh3, plot = FALSE, lag.max = 60)
plot(
  res3ACF$lag,
  res3ACF$acf,
  col = '#66A61E',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res3ACF$lag, rev(res3ACF$lag)),
  c(res3ACF$acf, rep(0, length(res3ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}

# Estimate the parameter posterior mean
mean(res1[grid])
mean(res2[grid])
mean(res3[grid])

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example2-lgss.RData")
}

