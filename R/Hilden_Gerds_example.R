# Example from Hilden and Gerds

require(mvnfast)
require(Matrix)
require(MASS)
require(caret)


PIP_K_rep_cv = function(data,K,type,alpha=0.05,reps,seed=1988){
  set.seed(seed)
  # Check if K is smaller than sample size
  if(K>nrow(data)){print("error: K should be less than or equal to n")}
  
  # Create K equally sized folds after randomly shuffling the data
  PIP_cv = c()
  mse0_CV = c()
  mse1_CV = c()
  for(l in 1:reps){
    yourData<-data[sample(nrow(data)),]
    cvIndex <- createFolds(factor(data$x1), K, returnTrain = T)
    
    #Perform K fold cross validation
    
    pip_cv = c()
    mse0_CV_sub=c()
    mse1_CV_sub=c()
    
    for(j in names(cvIndex)){
      trainData = yourData[cvIndex[[j]],]
      testData = yourData[-cvIndex[[j]],]
    
      if(type=="binomial"){
        mod1 = glm(y ~ x1+x2+x3, family="binomial", data=trainData, maxit = 5000)
        mod0 = glm(y ~ x1+x2, family="binomial", data=trainData, maxit = 5000)
      }
      
         pred0 = predict(mod0,testData,type="response")
         pred1 = predict(mod1,testData,type="response")
        

      pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )
      #pip_cv = c(pip_cv,mean((-testData$y*log(pred1)- (1-testData$y)*log(1-pred1)) < (-testData$y*log(pred0)- (1-testData$y)*log(1-pred0))) + 0.5*mean((-testData$y*log(pred1)- (1-testData$y)*log(1-pred1)) == (-testData$y*log(pred0)- (1-testData$y)*log(1-pred0))) )
      
      mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
      mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
    }
    PIP_cv = c(PIP_cv,mean(pip_cv))
    
    mse0_CV = mean(mse0_CV_sub)
    mse1_CV = mean(mse1_CV_sub)
  }
  return(list("PIP_cv"=mean(PIP_cv),"PIP_cv_lower"=quantile(PIP_cv,alpha),"PIP_cv_upper" = quantile(PIP_cv,1-alpha),"mse0"=mean(mse0_CV),"mse1"=mean(mse1_CV)))
}


expit = function(x){return(exp(x)/(1+exp(x)))}

do_sim = function(i){
  set.seed(i)
  sampsize=30000


  x1 = rbinom(30000,1,0.5)
  x2 = rbinom(30000,1,0.3)
  x3 = rnorm(sampsize)

  linear = -3+x1+x2+0.4*x3

  prob = expit(linear)

  y = rbinom(sampsize,1,prob)
  dat = data.frame(y,x1,x2,x3)
  pip_out = PIP_K_rep_cv(dat,5,"binomial",0.05,1,1988)
  return(list("pip"=pip_out$PIP_cv))
}


cl <- makeCluster(7)
registerDoSNOW(cl)
iterations <- 10000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
output <- foreach(i=1:iterations,.packages=c("MASS","Matrix","mvnfast","caret"),
                  .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                    result <- do_sim(i)
                    pip = result$pip
                    return(cbind(pip))
                  }
close(pb)
stopCluster(cl)

boxplot(output)
mean(output)


output_sq_loss = output
boxplot(output_sq_loss)
mean(output_sq_loss)



plot(output,output_sq_loss)
