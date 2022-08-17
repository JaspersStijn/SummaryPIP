# install.packages("LPS")
library(LPS)
require(mvnfast)
beta11 = -10

TRUE_PIPs_faster = function(beta01,beta11,sigma,prop1,sampsize){
  N = sampsize
  PIP_theoretical_true = pnorm(abs(beta11)*prop1/(2*sigma))* (1-prop1) + pnorm(abs(beta11)*(1-prop1)/(2*sigma))* (prop1)
  
  #Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(prop1*(1-prop1))/(sigma^2)/N),2,2)
  Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N),2,2)
  
  E1 <- matrix(nrow=10000000/2, ncol = 2)
  class(E1) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,-0.5*beta11), Sigma, A=E1)
  E2 <- matrix(nrow=10000000/2, ncol = 2)
  class(E2) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,0.5*beta11), Sigma, A=E2)
  
  PIP_expected_true =  (sum(E1[,1]^2< E1[,2]^2) + sum(E2[,1]^2< E2[,2]^2))/10000000
  
  N = 100000000
  #Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(prop1*(1-prop1))/(sigma^2)/N),2,2)
  Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N),2,2)
  
  E1 <- matrix(nrow=10000000/2, ncol = 2)
  class(E1) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,-0.5*beta11), Sigma, A=E1)
  E2 <- matrix(nrow=10000000/2, ncol = 2)
  class(E2) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,0.5*beta11), Sigma, A=E2)
  PIP_expected_true_limit =  (sum(E1[,1]^2< E1[,2]^2) + sum(E2[,1]^2< E2[,2]^2))/10000000
  
  return(list("PIP_theor" = PIP_theoretical_true, "PIP_exp" = PIP_expected_true,"PIP_exp_limit"=PIP_expected_true_limit ))
}


f = function(beta11){
sampsize=4000
set.seed(1988)
beta01 = 10
sigma = 2
prop1=0.5
# True value for m0 intercept

beta00 = beta01 + beta11*0.5

# Sample dataset with specified sample size
N = sampsize
X = c(rep(0,(1-prop1)*sampsize),rep(1,prop1*sampsize))
U <- rnorm(sampsize, sd = sigma)
Y = beta01 + beta11*X + U

training=data.frame(X, Y)

mod = lm(Y~X,data=training)
means = predict(mod,data.frame(X=c(0,1)))

y = seq(c(means + c(6,-6))[2],c(means + c(6,-6))[1],length=100)

f1 = dnorm(y,means[1],sigma*sqrt(1+2/sampsize))
f2 = dnorm(y,means[2],sigma*sqrt(1+2/sampsize))

plot(y,f2,type="l")
lines(y,f1,col="red")

PIPs = TRUE_PIPs_faster(beta01,beta11,sigma,prop1,sampsize)
Exp = PIPs$PIP_exp
The = PIPs$PIP_theor

title(paste( paste(paste("Theor", round(The,4),sep="="),paste("Exp", round(Exp,4),sep="="),sep="\n")   ,paste("OVL" , round(1-OVL(means, rep(sigma*sqrt(1+2/sampsize),2), cutoff=1e-4, n=1e6),4),sep="="),sep="\n"))
}





par(mfrow=c(1,3))
f(0)
f(-1)
f(-4)


