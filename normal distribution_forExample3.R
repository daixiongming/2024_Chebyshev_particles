#Normal distribution
setwd("~/Documents/Code/Minimal_Energy_Way/1323273/MED Rcodes/CameraPlacement")
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
N=60

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

Sigama_vi=3.5;
mu_vi=2;


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



D_matrix_3 <- matrix(0, nrow = NumberOfSize, ncol = 3)






# 
# 
# 
# Sigama1=0.98
# Sigama2=0.89
# Sigama3=0.98
# Sigama4=0.7
# mu1=2.5
# mu2=-2.5
# mu3=3
# mu4=-3
# 
# 
# Sigama12=0.87
# mu12=2
#   
#   
# Sigama11=0.96
# Sigama22=0.95
# Sigama33=0.99
# Sigama44=0.91
# mu11=2.5
# mu22=-2.5
# mu33=-3
# mu44=3
#   
# 
# X_1=D*Sigama1+mu1
# X_2=D*Sigama2+mu2
# X_3=D*Sigama3+mu3
# X_4=D*Sigama4+mu4
# 
# X12=append(X_1,X_2)
# X123=append(X12,X_3)
# X1234=append(X123,X_4)
# 
# 
# 
# Y_1=D*Sigama12+mu12
# 
# 
# Z_11=D*Sigama11+mu11
# Z_22=D*Sigama22+mu22
# Z_33=D*Sigama33+mu33
# Z_44=D*Sigama44+mu44
# 
# Z12=append(Z_11,Z_22)
# Z123=append(Z12,Z_33)
# Z1234=append(Z123,Z_44)
# 
# library("scatterplot3d")
# 
# NumberOfSize=200
# 
# # r sample multiple times without replacement
# #sample (c(1:10), size=3, replace =F)
# 
# # r sample with replacement from vector
# #sample (c(1:10), size=3, replace=T)
# 
# X_1=sample(X_1,size=NumberOfSize)
# X_2=sample(X_2,size=NumberOfSize)
# X_3=sample(X_3,size=NumberOfSize)
# X_4=sample(X_4,size=NumberOfSize)
# 
# X_12=append(X_1,X_2)
# X_123=append(X_12,X_3)
# X_1234=append(X_123,X_4)
# X=X_1234
# 
# Y_1=sample(Y_1,size=NumberOfSize)
# Y_12=append(Y_1,Y_1)
# Y_123=append(Y_12,Y_1)
# Y_1234=append(Y_123,Y_1)
# Y=Y_1234
# 
# Z_11=sample(Z_11,size=NumberOfSize)
# Z_22=sample(Z_22,size=NumberOfSize)
# Z_33=sample(Z_33,size=NumberOfSize)
# Z_44=sample(Z_44,size=NumberOfSize)
# 
# Z_1122=append(Z_11,Z_22)
# Z_112233=append(Z_1122,Z_33)
# Z_11223344=append(Z_112233,Z_44)
# Z=Z_11223344
# 
# scatterplot3d(X,Y,Z, color=par("col"),pch = 16, grid=FALSE, box=FALSE)
# 
# 
# #save(X,file="x_r_output.txt")
# 
# sink("x_r.txt")
# cat(X)
# 
# write.table(X,"X_table_r.txt",sep="\n",row.names = FALSE)
# 
# sink("y_r.txt")
# cat(Y)
# write.table(Y,"Y_table_r.txt",sep="\n",row.names = FALSE)
# 
# 
# 
# sink("z_r.txt")
# cat(Z)
# write.table(Z,"Z_table_r.txt",sep="\n",row.names = FALSE)
# 
# 

#c->file("x_r.txt","w")
#writeLines(X,c)
#write.table(X,"/home/daixiongming/Documents/Code/Minimal_Energy_Way/1323273/MED Rcodes/CameraPlacement")

#scatterplot3d(X,Z,Y, color=par("col"),pch = 16, grid=FALSE, box=FALSE)
#scatterplot3d(X,Y,Z, pch = 16, grid=c("xy","xz","yz"))
