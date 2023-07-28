##############################################################################
# State estimation in a LGSS model using particle and Kalman filters
#
# Johan Dahlin <uni (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################
setwd("/home/daixiongming/Documents/Code/pmh-tutorial-2.1/r_chebyshev")
source("helpers/dataGeneration.R")
source("helpers/stateEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

# Save plot to file
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
x <- data$x
y <- data$y

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("figures/example1-lgss.pdf",
            height = 10,
            width = 8)
}

# Plot the latent state and observations
layout(matrix(c(1, 1, 2, 2, 3, 4), 3, 2, byrow = TRUE))
par   (mar = c(4, 5, 0, 0))

grid <- seq(0, T)

plot(
  grid,
  y,
  col = "#1B9E77",
  type = "l",
  xlab = "time",
  ylab = "observation",
  ylim = c(-6, 6),
  bty = "n"
)
polygon(c(grid, rev(grid)),
        c(y, rep(-6, T + 1)),
        border = NA,
        col = rgb(t(col2rgb("#1B9E77")) / 256, alpha = 0.25))



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


ac=1 # 1 is also good! than 2.
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
D200=D
pdf(file="normal_paper.pdf")
hist(D,breaks=15, xlim=c(-3,3), prob=T, main="Normal (0,1)" ,xlab="x", ylab="density", cex.lab=1.25)
curve(dnorm(x,0,1), from=-3,to=3,col=2,lty=2, lwd=2, add=T)
dev.off()

qqnorm(D)

pdf(file="QQnorm200.pdf")
plot(qnorm(((1:N)-.5)/N),sort(D),xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="200 Chebyshev Particles",col="#9933FF")
abline(0,1,col="#CC0066",lty=2)
dev.off()
###############################################################################################
###############################################################################################
###############################################################################################
##############################################
noParticles <- 200
EnergyParticles<- seq(1:noParticles) 
#################################################
###############################################################################################
###############################################################################################





##############################################################################
# State estimation using the particle filter and Kalman filter
##############################################################################
if (loadSavedWorkspace) {
  load("savedWorkspaces/example1-lgss.RData")
} else {
  # Using noParticles = 20 particles and plot the estimate of the latent state
  #noParticles <- 20
  outputPF <-
    particleFilterD(y, c(phi, sigmav, sigmae), noParticles, initialState,D)
  outputKF <-
    kalmanFilter(y, c(phi, sigmav, sigmae), initialState, 0.01)
  difference <-
    outputPF$xHatFiltered - outputKF$xHatFiltered[-(T + 1)]
}

grid <- seq(0, T - 1)
plot(
  grid,
  difference,
  col = "#7570B3",
  type = "l",
  xlab = "time",
  ylab = "error in state estimate",
  ylim = c(-0.1, 0.1),
  bty = "n"
)
polygon(
  c(grid, rev(grid)),
  c(difference, rep(-0.1, T)),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)

# Compute bias and MSE
logBiasMSE <- matrix(0, nrow = 7, ncol = 2)
gridN <- c(10, 20, 50, 100, 200, 500, 1000)

for (ii in 1:length(gridN)) {
  pfEstimate <-
    particleFilter(y, c(phi, sigmav, sigmae), gridN[ii], initialState)
  pfEstimate <- pfEstimate$xHatFiltered
  kfEstimate <- outputKF$xHatFiltered[-(T + 1)]

  logBiasMSE[ii, 1] <- log(mean(abs(pfEstimate - kfEstimate)))
  logBiasMSE[ii, 2] <- log(mean((pfEstimate - kfEstimate) ^ 2))
}

##############################################################################
# Plot the bias and MSE for comparison
##############################################################################


plot(
  gridN,
  logBiasMSE[, 1],
  #col = "#E7298A",
  col="#009999",   #http://www.sthda.com/english/wiki/colors-in-r
  type = "l",
  xlab = "no. particles (N)",
  ylab = "log-bias",
  ylim = c(-7,-3),
  bty = "n"
)
polygon(
  c(gridN, rev(gridN)),
  c(logBiasMSE[, 1], rep(-7, length(gridN))),
  border = NA,
  col = rgb(t(col2rgb("#009999")) / 256, alpha = 0.25)
)
#points(gridN, logBiasMSE[, 1], col = "#009999", pch = 19)
points(gridN, logBiasMSE[, 1], col = "#009999", pch = 21)
#https://www.google.com/search?q=How+do+you+fill+PCH+in+R?&sxsrf=AOaemvKrDUUAK3jbZAqv83x44nUrIyInLw:1635209574408&tbm=isch&source=iu&ictx=1&fir=J2JAWpPtTvK7JM%252C-_hlpsv9txrjuM%252C_&vet=1&usg=AI4_-kRQJYD6mnkCePhOxSCKIgZZAjhRTg&sa=X&ved=2ahUKEwiLxcKu7ubzAhXSk2oFHapiDZUQ9QF6BAgUEAE&biw=2987&bih=723&dpr=1#imgrc=J2JAWpPtTvK7JM
plot(
  gridN,
  logBiasMSE[, 2],
 # col = "#66A61E",
  col="#0000FF",
  lwd = 1.5,
  type = "l",
  xlab = "no. particles (N)",
  ylab = "log-MSE",
  ylim = c(-12,-6),
  bty = "n"
)
polygon(
  c(gridN, rev(gridN)),
  c(logBiasMSE[, 2], rep(-12, length(gridN))),
  border = NA,
  col = rgb(t(col2rgb("#0000FF")) / 256, alpha = 0.25)
)
#points(gridN, logBiasMSE[, 2], col = "#0000FF", pch = 19)
points(gridN, logBiasMSE[, 2], col = "#0000FF", pch = 22)

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}

# Print a table (no. particles, log-bias, log-mse)
print(t(rbind(gridN, t(logBiasMSE))))

# gridN
# [1,]    10 -3.696997  -6.938594
# [2,]    20 -3.964671  -7.493297
# [3,]    50 -4.567552  -8.718346
# [4,]   100 -4.850363  -9.294468
# [5,]   200 -5.192173  -9.905719
# [6,]   500 -5.668407 -10.866745
# [7,]  1000 -6.077648 -11.671646

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example1-lgss.RData")
}

