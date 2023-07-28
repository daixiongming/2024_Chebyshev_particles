setwd("/home/daixiongming/Documents/Code/pmh-tutorial-2.1/r_chebyshev")

load("savedWorkspaces/example1-lgss_40_chebyshev.RData")


load("savedWorkspaces/example1-lgss_120_chebyshev.RData")



load("savedWorkspaces/example1-lgss_200_chebyshev.RData")


pdf(file="QQnormtotal.pdf")
par(mfrow=c(3,1))
N_40=40
N_120=120
N_200=200
plot(qnorm(((1:N_40)-.5)/N_40),sort(D40),xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="40 Chebyshev Particles",col="#9933FF")
abline(0,1,col="#CC0066",lty=2)

plot(qnorm(((1:N_120)-.5)/N_120),sort(D120),xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="120 Chebyshev Particles",col="#9933FF")
abline(0,1,col="#CC0066",lty=2)

plot(qnorm(((1:N_200)-.5)/N_200),sort(D200),xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="200 Chebyshev Particles",col="#9933FF")
abline(0,1,col="#CC0066",lty=2)



dev.off()



