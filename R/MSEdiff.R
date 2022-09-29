# MSE diff vs PIP

PIP_K_rep_cv = function(data,K,reps,seed=1988,alpha=0.05,plot=FALSE){
  set.seed(seed)
  
  # Create K equally sized folds after randomly shuffling the data
  PIP_cv = c()
  mse0_CV = c()
  mse1_CV = c()
  for(l in 1:reps){
    yourData<-data[sample(nrow(data)),]
    cvIndex <- createFolds(data$y, K, returnTrain = T)
  
    #Perform K fold cross validation
    
    pip_cv = c()
    mse0_CV_sub=c()
    mse1_CV_sub=c()
    
    for(j in names(cvIndex)){
      trainData = yourData[cvIndex[[j]],]
      testData = yourData[-cvIndex[[j]],]
      
      
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
      
      pred0 = predict(mod0,testData)
      pred1 = predict(mod1,testData)
      
      pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )
      
      mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
      mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
    }
    if(plot==TRUE){
      dat_plot = data.frame("pip"=pip_cv,"mse_diff"=mse1_CV_sub-mse0_CV_sub)
      dat_plot = dat_plot[order(dat_plot$pip),]
    if(l==1){plot(dat_plot$pip,dat_plot$mse_diff,col=l,type="l",ylim=c(-1,2),xlim=c(0,1))}
    if(l>1){lines(dat_plot$pip,dat_plot$mse_diff,col=l)}}
    PIP_cv = c(PIP_cv,mean(pip_cv))
    
    mse0_CV = mean(mse0_CV_sub)
    mse1_CV = mean(mse1_CV_sub)
  }
  return(list("PIP_cv"=mean(PIP_cv),"PIP_cv_lower"=quantile(PIP_cv,alpha),"PIP_cv_upper" = quantile(PIP_cv,1-alpha),"mse0"=mean(mse0_CV),"mse1"=mean(mse1_CV)))
}


cv_out = PIP_K_rep_cv(dat_in,5,50,plot=TRUE)

f = function(i,beta11,sampsize){
  set.seed(i)
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
  mod0 = lm(Y~1,data=training)
  
  The = pnorm(abs(coef(mod)[2])/(4*sigma))
  
  X_star = c(rep(0,(1-prop1)*1000000),rep(1,prop1*1000000))
  Y_star = beta01 + beta11*X_star + rnorm(1000000, sd = sigma)
  dat_star = data.frame("X"=X_star,"y"=Y_star)
  
  pred1 = predict(mod,dat_star)
  pred0 = predict(mod0,dat_star)
  true_mse_diff = mean((dat_star$y-pred1)^2) - mean((dat_star$y-pred0)^2) 
  
  # dat_star=data.frame("X"=X_star, Y_star)
  
  mse_diff = summary(mod)$sigma^2 -summary(mod0)$sigma^2 
  pval = coef(summary(mod))[2,4]
  
  The_emp = mean(((dat_star$y-pred1)^2) < ((dat_star$y-pred0)^2)) + 0.5*mean(((dat_star$y-pred1)^2) == ((dat_star$y-pred0)^2)) 
  
  dat_in = training
  names(dat_in) = c("x","y")
  cv_out = PIP_K_rep_cv(dat_in,5,50)

  return(list("The" = The,"The_emp" = The_emp,"MSE_diff"=mse_diff,"MSE_diff_emp"=true_mse_diff,"pval"=pval,"beta11_hat"=coef(mod)[2],"PIP_CV" = cv_out$PIP_cv,"MSE_diff_cv" = cv_out$mse1-cv_out$mse0 ))
}

require(doSNOW)
require(ggplot2)

# par(mfrow=c(1,1))
# cl <- makeCluster(7)
# 
# registerDoSNOW(cl)
# iterations <- 1000
# pb <- txtProgressBar(max = iterations, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# output <- foreach(beta11=seq(0,15,length=100),.packages=c("MASS","Matrix","mvnfast","LPS"),
#                   .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
#                     result <- f(1988,beta11)
#                     pip_The =  result$The
#                     mse_diff = result$MSE_diff
#                     pval = result$pval
#                     return(cbind(pip_The,mse_diff,pval))
#                   }
# close(pb)
# stopCluster(cl)


par(mfrow=c(1,1))
cl <- makeCluster(7)

registerDoSNOW(cl)
iterations <- 1000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
output <- foreach(i=1:iterations,.packages=c("MASS","Matrix","mvnfast","LPS","caret"),
                  .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                    result <- f(i,-1,40)
                    pip_The =  result$The
                    mse_diff = result$MSE_diff
                    pval = result$pval
                    pip_The_emp =  result$The_emp
                    mse_diff_emp = result$MSE_diff_emp
                    beta1_est = result$beta11_hat
                    pip_cv = result$PIP_CV
                    mse_diff_cv = result$MSE_diff_cv
                    return(cbind(pip_The,mse_diff,pval,pip_The_emp,mse_diff_emp,beta1_est,pip_cv,mse_diff_cv))
                  }
close(pb)
stopCluster(cl)

output_keep=output

save(output,file=paste0(paste(paste0("/Volumes/GoogleDrive/My Drive/SummaryPIP/R/Output_Sims/investigate_emp",1),paste0("sampsize",40),sep="_"),".R"))


rel_mse = function(x){
  return(-4*4*qnorm(x)^2)
}

par(mfrow=c(1,2))
plot(output[,"pip_The"],output[,"mse_diff"],col=c("darkgreen","red")[(output[,"beta1_est"]>-1)+1],xlab="pip_est_plugin",ylab="MSE_diff")
lines(seq(0.5,0.999,length=1000),sapply(seq(0.5,0.999,length=1000),rel_mse),col="red")
plot(output[,"pip_The_emp"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]>-1)+1])
plot(output[,"pip_The"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]>-1)+1])





plot(output[,"pip_The"],output[,"pip_The_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]>-1)+1])


rel_mse = function(x){
  return(-4*4*qnorm(x)^2)
}

rel_mse2 = function(x,beta){
  return(-4*4*qnorm(x)^2+4^2/40)
}
sampsize=200000
lines(seq(0.5,0.999,length=1000),sapply(seq(0.5,0.999,length=1000),rel_mse),col="red")
lines(seq(0.5,0.999,length=1000),sapply(seq(0.5,0.999,length=1000),rel_mse2),col="darkgreen")


plot(output[,"pval"],output[,"pip_The"])

pvals = output[,"pval"]
f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-2))
  )
}
lines(seq(0,max(max(pvals),4e-10),length=1000),sapply(seq(0,max(max(pvals),4e-10),length=1000),f_pip,n=sampsize),col="red",lwd=2)


plot(output[,"pval"],output[,"mse_diff"])




# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
# unregister_dopar()


plot(seq(0.5,0.9999,length=1000),sapply(seq(0.5,0.9999,length=1000),rel_mse2,beta=beta11))


output[which((output[,"mse_diff"]>0.4)&(output[,"pip_The"]>0.6)),]

i=17
summary(mod0)









output[output[,"pip_The"]<0.51,]




pred1 = beta01+beta11*dat_star$X
pred0 = beta00
The_emp = mean(((dat_star$y-pred1)^2) < ((dat_star$y-pred0)^2)) + 0.5*mean(((dat_star$y-pred1)^2) == ((dat_star$y-pred0)^2)) 
The_emp;pnorm(abs(beta11)/(4*sigma))






beta=0.25
sampsize=200000
boxplot(log(output[,"pval"])/sampsize)
abline(h=-0.5*log(1+beta^2/(4*2^2)),col="red",lwd=2)
title(paste(paste0("Beta: ",beta),paste0("Sampsize: ",sampsize),sep="\n"))
legend("bottomright","-0.5log(1+beta^2/(4*sigma^2))",lwd=2,col="red")





10 - coef(mod)[1] ; 10 - coef(mod0)

9 - (coef(mod)[1]+coef(mod)[2]); 9 - coef(mod0)




PIP_K_rep_cv = function(data,K,reps,seed=1988,alpha=0.05){
  set.seed(seed)
  
  # Create K equally sized folds after randomly shuffling the data
  PIP_cv = c()
  mse0_CV = c()
  mse1_CV = c()
  for(l in 1:reps){
    yourData<-data[sample(nrow(data)),]
    cvIndex <- createFolds(data$y, K, returnTrain = T)
    
    #Perform K fold cross validation
    
    pip_cv = c()
    mse0_CV_sub=c()
    mse1_CV_sub=c()
    
    for(j in names(cvIndex)){
      trainData = yourData[cvIndex[[j]],]
      testData = yourData[-cvIndex[[j]],]
      
      
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
      
      pred0 = predict(mod0,testData)
      pred1 = predict(mod1,testData)
      
      pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )
      
      mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
      mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
    plot()
    }
    
    PIP_cv = c(PIP_cv,mean(pip_cv))
    
    mse0_CV = mean(mse0_CV_sub)
    mse1_CV = mean(mse1_CV_sub)
  }
  return(list("PIP_cv"=mean(PIP_cv),"PIP_cv_lower"=quantile(PIP_cv,alpha),"PIP_cv_upper" = quantile(PIP_cv,1-alpha),"mse0"=mean(mse0_CV),"mse1"=mean(mse1_CV)))
}

library(caret)
dat_in = training
names(dat_in) = c("x","y")
PIP_K_rep_cv(dat_in,5,50)


output = output_keep


#-------------
# Seed 17, sampsize 40, beta=-1 shows strange results.
#-------------
i = 17
f(i,-1,40)  # Significant p-value, PIP is 0.63, and below the CV PIP is also 0.67

set.seed(i)
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
dat_in = training
names(dat_in) = c("x","y")
PIP_K_rep_cv(dat_in,5,50)

# But: based on dat_star, the amount with which model m1 overshoots when x = 0 and undershoots when x=1 is higher as compared o model m0:
# Hence the empirical PIP and empirical MSE diff show that m0 is better
mod = lm(Y~X,data=training)
mod0 = lm(Y~1,data=training)
10 - coef(mod)[1] ; 10 - coef(mod0)
9 - (coef(mod)[1]+coef(mod)[2]); 9 - coef(mod0)


plot(training$X,training$Y)


colnames(output)
par(mfrow=c(2,2))
plot(output[,"beta1_est"],output[,"pval"])
plot(output[,"beta1_est"],output[,"pip_The"])
plot(output[,"beta1_est"],output[,"pip_The_emp"])

par(mfrow=c(2,2))
plot(output[,"beta1_est"],output[,"mse_diff_emp"])
plot(output[,"beta1_est"],output[,"mse_diff"])
plot(output[,"beta1_est"],output[,"mse_diff_cv"])



par(mfrow=c(1,1))

plot(output[,"mse_diff"],output[,"mse_diff_cv"],col=c("darkgreen","red")[(output[,"mse_diff_emp"]<0)+1],pch=19)
abline(a=0,b=1,col="red",lwd=2)

plot(output[,"mse_diff"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
plot(output[,"mse_diff_cv"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"mse_diff_emp"]<0)+1],pch=19)


plot(output[,"pip_The"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
abline(h=0)

plot(output[,"pip_cv"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
abline(h=0)

plot(output[,"pip_cv"],output[,"pip_The_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
abline(v=pnorm(1/8))


lim_new = function(PIP){
  return(-0.5*log(1+4*qnorm(PIP)^2))
}



plot(seq(-10,10,length=10000),pnorm(abs(seq(-10,10,length=10000))/8))



PIP = seq(0,1,length=10000)
plot(PIP,sapply(PIP,lim_new),ylab="Limit for n goining to infinity",type="l")

lim_new(0.9994995)
plot(sapply(PIP,lim_new),PIP,type="l")



plot(output[,"pip_The_emp"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)

###
par(mfrow=c(2,2))
plot(output[,"mse_diff"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
plot(output[,"mse_diff_cv"],output[,"mse_diff_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
plot(output[,"pip_The"],output[,"pip_The_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)
plot(output[,"pip_cv"],output[,"pip_The_emp"],col=c("darkgreen","red")[(output[,"beta1_est"]<0)+1],pch=19)

##
par(mfrow=c(2,2))
plot(output[,"beta1_est"],output[,"pval"])
plot(output[,"beta1_est"],output[,"pip_The"])
plot(output[,"beta1_est"],output[,"pip_The_emp"])

par(mfrow=c(2,2))
plot(output[,"beta1_est"],output[,"mse_diff_emp"])
plot(output[,"beta1_est"],output[,"mse_diff"])
plot(output[,"beta1_est"],output[,"mse_diff_cv"])



plot(density(output[output[,"pip_The_emp"]<0.5,"beta1_est"]))

min(output[output[,"pip_The_emp"]>0.5,"beta1_est"])
max(output[output[,"pip_The_emp"]>0.5,"beta1_est"])




load(paste0(paste(paste0("/Volumes/GoogleDrive/My Drive/SummaryPIP/R/Output_Sims/investigate_emp",1),paste0("sampsize",40),sep="_"),".R"))


