# install.packages("LPS")
library(LPS)
require(mvnfast)
beta11 = 0

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


f = function(i,beta11){
  bill=c()
  bill2=c()
  bill3=c()
  the = c()
  exp = c()
  for(i in 1:1){
  set.seed(i)
  sampsize=20
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
  
  # plot(y,f2,type="l",lwd=2)
  # lines(y,f1,col="red",lwd=2,lty=2)
  # title(paste0("Effect Size= ",beta11))
  #PIPs = TRUE_PIPs_faster(beta01,beta11,sigma,prop1,sampsize)
  #Exp = PIPs$PIP_exp
  The = pnorm(abs(coef(mod)[2])/(4*sigma))
  Exp = 0
  
  bill = c(bill,1-OVL(means, rep(sigma*sqrt(1+2/sampsize),2), cutoff=1e-4, n=1e6))
  the = c(the,The)
  exp = c(exp,Exp)
  bill2 = c(bill2,1-pnorm(0.5*(means[1]+means[2]),min(means),sigma*sqrt(1+2/sampsize))+pnorm(0.5*(means[2]+means[1]),max(means),sigma*sqrt(1+2/sampsize)))

  OVL_new = 2*(1-pnorm((abs(means[2]-means[1]))/(sqrt(2*(sigma^2*(1+2/sampsize))))))

  
  bill3 = c(bill3,OVL_new)
  
  #title(paste( paste(paste("Theor", round(The,4),sep="="),paste("Exp", round(Exp,4),sep="="),sep="\n")   ,paste("OVL" , round(1-OVL(means, rep(sigma*sqrt(1+2/sampsize),2), cutoff=1e-4, n=1e6),4),sep="="),sep="\n"))
  }
  
  return(list("The" = mean(the),"Exp"=mean(exp),"Ovl1" =mean(bill),"Ovl2" =mean(bill2),"Ovl3" =mean(bill3) ))
}





par(mfrow=c(1,3))
f(1988,0)
f(1988,-1)
f(1988,-4)





require(doSNOW)
require(ggplot2)

par(mfrow=c(1,1))
cl <- makeCluster(7)
registerDoSNOW(cl)
iterations <- 100
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
output <- foreach(beta11=seq(0,15,length=100),.packages=c("MASS","Matrix","mvnfast","LPS"),
                  .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                    result <- f(1988,beta11)
                    pip_The =  result$The
                    pip_Exp =  result$Exp
                    Bill_Ovl1 =  result$Ovl1
                    Bill_Ovl2 =  result$Ovl2
                    Bill_Ovl3 =  result$Ovl3
                    return(cbind(pip_The,pip_Exp,Bill_Ovl1,Bill_Ovl2,Bill_Ovl3))
                  }
close(pb)
stopCluster(cl)

par(mfrow=c(1,2))
plot(output[,"pip_The"],output[,"Bill_Ovl2"],type="l",lwd=2,xlim=c(0.45,1),xlab="plugin Theoretical PIP",ylab="OVL1")  
lines(output[,"pip_Exp"],output[,"Bill_Ovl2"],col="red",lty=2,lwd=2)    
sampsize=20
rel = function(pip){
  return(2*pnorm(-2/sqrt(1+2/sampsize)*qnorm(pip)))
}
lines(seq(0.5,0.999,length=100),sapply(seq(0.5,0.999,length=100),rel),col="darkgreen",lwd=2,lty=3)  
legend("topright",c("Simulation","Theoretical"),lty=c(1,2),lwd=2,col=c("black","darkgreen"))

plot(output[,"pip_The"],output[,"Bill_Ovl3"],type="l",lwd=2,xlim=c(0.45,1),xlab="plugin Theoretical PIP",ylab="OVL2")  
lines(output[,"pip_Exp"],output[,"Bill_Ovl3"],col="red",lty=2,lwd=2)    

rel2 = function(pip){
  return(2*(1-pnorm(4/sqrt(2*(1+2/sampsize))*qnorm(pip))))
}

lines(seq(0.5,0.999,length=100),sapply(seq(0.5,0.999,length=100),rel2),col="darkgreen",lwd=2,lty=3)  

legend("topright",c("Simulation","Theoretical"),lty=c(1,2),lwd=2,col=c("black","darkgreen"))





pnorm(-0.5*qnorm(0.5*0.5))


x = seq(-10,10,length=1000)
y = seq(-10,10,length=1000)

par(mfrow=c(1,1))
plot(x,dnorm(x,0,1),type="l")
lines(y,dnorm(y,0.164,1),col="red")
abline(v=(0+0.164)/2,lty=2)

lines(x,dnorm(x,0,1)*dnorm(y,0.164,1),col="darkgreen")
lines(x,dnorm(y,0,1)*dnorm(x,0.164,1),col="darkblue")


f = function(x,y){return(min(dnorm(x)*dnorm(y,0.164,1),dnorm(y)*dnorm(x,0.164,1)))}





integrate(function(y) { 
  sapply(y, function(y) {
    integrate(function(x) f(x,y), -10, 10)$value
  })
}, -10, 10)





fx0 = function(x){return(f(x,0))}
plot(x,sapply(x,fx0),type="l")
fx1 = function(x){return(f(x,1))}
lines(x,sapply(x,fx1),col="red")
fx2 = function(x){return(f(x,2))}
lines(x,sapply(x,fx2),col="blue")






means = c(10,9.5)
sigma=2
test = function(x,y){return(abs(dnorm(x,means[1],sigma)*dnorm(y,means[2],sigma)-dnorm(y,means[1],sigma)*dnorm(x,means[2],sigma)))}

OVL_new = 1 - 0.5*(integrate(function(y) { 
  sapply(y, function(y) {
    integrate(function(x) test(x,y), -200, 200)$value
  })
}, -200, 200)$value)











