##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a SV model
#
# Johan Dahlin <uni (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################
setwd("/home/daixiongming/Documents/Code/pmh-tutorial-2.1/r_chebyshev")

#setwd("~/Documents/Code/pmh-tutorial-2.1/r")
library("Quandl")
library("mvtnorm")
source("helpers/stateEstimation.R")
source("helpers/parameterEstimation.R")
source("helpers/plotting.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

loadSavedWorkspace <- FALSE


# Save plot to file
#savePlotToFile <- FALSE
savePlotToFile <- TRUE

#nPlot <- 2500
nPlot <- 1200

#Quandl.api_key("RAznaYRym2wbozfQTddG")
#Quandl.api_key("CyzgihyVzJC3pYHR16Jt")  #https://www.quandl.com/tools/r 
Quandl.api_key("kHAzgtSa1bSyAC-norNn")  #https://www.quandl.com/tools/r 

##############################################################################
# Load data
##############################################################################
d <-
  Quandl(
    "NASDAQOMX/OMXS30",
    start_date = "2012-02-04", #2012-2014
    end_date = "2014-02-04",
    type = "zoo"
  )

#d<-Quandl("NASDAQOMX/OMXS30", api_key="RAznaYRym2wbozfQTddG")
y <- as.numeric(100 * diff(log(d$"Index Value")))

##############################################################################
# PMH
##############################################################################
initialTheta <- c(0, 0.9, 0.2)
#initialTheta <- c(0, 0.8, 0.1)

noParticles <- 500
noBurnInIterations <- 2500
#noBurnInIterations <- 100

noIterations <- 7500  #original
#noIterations <- 2000
#noIterations <- 3200

#stepSize <- diag(c(0.10, 0.01, 0.05) ^ 2)
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
p=1
k=4*p
k=80
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


ac=1# 1 is also good!
k=80
m=40  #k must the 2d times than m!!! please keep in mind that!
W=1
energy=function(x)
{
  d=abs(D-x)+10^(-10)
  # d_kernel=kernel_CF(d)
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
  #d_kernel=kernel_CF(d)
  #val=-2*log(f(x))+log(sum(1/(fev^2*d^k)))
  val=-(k/(2*p))*lf(x)+log(sum(exp(-(k/(2*p))*lfev+k*log(d_kernel))))
  A=k*log(d_kernel)
  B=ac*d_kernel
  C=(-1)*lfev*lf(x)
  D=C+B
  val_w=log(sum(exp(-(m)*D+A)))
  return(val_w)
}



#N=100
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


################################################################
###############################################################
#############################################################

# elif D == 3:
#   # set gaussian ceters and covariances in 3D
#   means = np.array([[2, 2, 2],
#                     [2, 2, -2],
#                     [-2, 2, -2],
#                     [-2, 2, 6],
#                     [-6, 2, 6],
#                     [-6, 2, 2],
#                     [6, 2, 2],
#                     [6, 2, -2],
#                     [6, 2, 6],
#                     [6, 2, -6],
#                     [-6, 2, -2],
#                     [-6, 2, -6],
#                     [-2, 2, -6],
#                     [-2, 2, 2],
#                     [2, 2, 6],
#                     [2, 2, -6]])
# covs = np.array([np.diag([0.63, 0.74, 1]),
#                  np.diag([0.68, 0.61, 1]),
#                  np.diag([0.72, 0.81, 0.91]),
#                  np.diag([0.81, 0.84, 0.81]),
#                  np.diag([0.97, 1, 0.65]),
#                  np.diag([1, 0.85, 0.61]),
#                  np.diag([1, 0.91, 0.71]),
#                  np.diag([0.61, 0.75, 0.81]),
#                  np.diag([0.88, 0.81, 0.91]),
#                  np.diag([0.71, 0.95, 1.0]),
#                  np.diag([0.68, 0.88, 0.68]),
#                  np.diag([0.91, 0.75, 0.69]),
#                  np.diag([0.98, 0.88, 0.77]),
#                  np.diag([1.0, 0.75, 0.91]),
#                  np.diag([0.68, 0.82, 0.71]),
#                  np.diag([0.92, 0.95, 0.81]),
#                  np.diag([1.0, 0.95, 0.88]),
#                  np.diag([1.0, 1.0, 0.68])])

##############################################################
##############################################################
##############################################################
##############################################################
Sigama_mu=1;
mu_mu=0;

Sigama_fa=0.05;
mu_fa=0.95;

Sigama_vi=0.03;
mu_vi=0.2;


Mu_1=D*Sigama_mu+mu_mu;
Fa_2=D*Sigama_fa+mu_fa;
Vi_3=D*Sigama_vi+mu_vi;

# r sample multiple times without replacement
#sample (c(1:10), size=3, replace =F)

# r sample with replacement from vector
#sample (c(1:10), size=3, replace=T)
NumberOfSize=1600



X_1=sample(Mu_1,size=NumberOfSize,replace=T)
X_2=sample(Fa_2,size=NumberOfSize,replace=T)
X_3=sample(Vi_3,size=NumberOfSize,replace=T)


X_12=cbind(X_1,X_2)
X_123=cbind(X_12,X_3)



####################################################################################################################
#########################################################################################################################################################################################################################
####################################################################################################################
####################################################################################################################


#stepSize<-Sigmainv/100
#stepSize<-Sigmainv/10
stepSize <- diag(c(0.10, 0.01, 0.05) ^ 2)
#stepSize <- diag(c(1, 1,1) ^ 2)


if (loadSavedWorkspace) {
  load("savedWorkspaces/example3-sv_XMD.RData")
} else {
  res <- particleMetropolisHastingsSVmodel(y, initialTheta, noParticles, noIterations, stepSize,X_123)
}
savePlotToFile <- TRUE # why I need to put it here so that it can save?
##############################################################################
# Plot the results
##############################################################################
if (savePlotToFile) {
  cairo_pdf("figures/example3-sv_XMD.pdf",
            height = 10,
            width = 8)
}

iact <- makePlotsParticleMetropolisHastingsSVModel(y, res, noBurnInIterations, noIterations, nPlot)

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}

# Print the estimate of the posterior mean and standard deviation
resTh <- res$theta[noBurnInIterations:noIterations, ]
thhat   <- colMeans(resTh)
thhatSD <- apply(resTh, 2, sd)

print(thhat)
print(thhatSD)

#[1] -0.2337134  0.9708399  0.1498914
#[1] 0.37048000 0.02191359 0.05595271

# Compute an estimate of the IACT using the first 100 ACF coefficients
print(iact)
# [1] 135.19084  85.98935  65.80120

# Estimate the covariance of the posterior to tune the proposal
estCov <- var(resTh)
#               [,1]          [,2]          [,3]
# [1,]  0.137255431 -0.0016258103  0.0015047492
# [2,] -0.001625810  0.0004802053 -0.0009973058
# [3,]  0.001504749 -0.0009973058  0.0031307062

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example3-sv_XMD.RData")
}

################ test program ##########################

sigma <- matrix(c(3,2,2,6), 2, 2)
mu <- c(5,10)
x <- rmvnorm(1000, mean = mu, sigma = sigma)
head(x)
summary(x)
plot(x[,1], x[,2])


sigma <- matrix(c(3,2,1,2,6,2,1,2,1), 3, 3)
mu <- c(5,10,1)
x <- rmvnorm(1000, mean = mu, sigma = sigma)
head(x)
summary(x)
plot(x[,1], x[,2])
pairs(x)

###########This function should be simulated#########################
stepSize <- diag(c(0.10, 0.01, 0.05) ^ 2)
mu2=c(0,0.9,0.1)
x <- rmvnorm(1000, mean = mu2, sigma = stepSize) # square already have!
head(x)
summary(x)
plot(x[,1], x[,2])
pairs(x)


########### True Test ########################
#Sigma
x <- rmvnorm(1000, mean = mu2, sigma = Sigma)
head(x)
summary(x)
plot(x[,1], x[,2])
pairs(x)

