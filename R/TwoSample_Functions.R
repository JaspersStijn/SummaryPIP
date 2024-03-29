# Function to calculate the K-fold cross-validation PIP
#

require(mvnfast)
require(Matrix)
require(MASS)
require(caret)

# Helper functions

logit = function(x){
  return(log(x/(1-x)))
}
expit = function(x){
  return((exp(x)/(1+exp(x))))
}

#---------------
# Main functions
#---------------

# Theoretical PIP


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


## Empirical PIP

Emp_PIP = function(dataset,total_n,prop1,beta00,beta01,beta11,sigma){

  xstar = c(rep(0,(1-prop1)*total_n),rep(1,prop1*total_n))
  ustar <- rnorm(total_n, sd = sigma)
  ystar = beta01 + beta11*xstar  + ustar
  dat_star = data.frame('X'=xstar,'y'= ystar)


  # Empirical version

  sub = dataset
  N = nrow(sub)
  # Fit models m0 and m1

  mod1 = lm(Y ~ X, data = sub)
  mod0 = lm(Y ~ 1, data = sub)

  ### plug-in

  pred0 = predict(mod0,dat_star)
  pred1 = predict(mod1,dat_star)

  pip_emp = mean((pred1-dat_star$y)^2 < (pred0-dat_star$y)^2) + 0.5*mean((pred1-dat_star$y)^2 == (pred0-dat_star$y)^2)

  ### expected

  coef_mod1 = coef(mod1)
  coef_mod0 = coef(mod0)

  var_covar0 = vcov(mod0)
  var_covar1 = vcov(mod1)
  var_b01 = var_covar1[1,1]
  var_b11 = var_covar1[2,2]
  covar_mod1 = var_covar1[1,2]

  sigma_pars = matrix(c(sigma^2/sampsize,sigma^2/sampsize,0,sigma^2/sampsize,2*sigma^2/sampsize,-2*sigma^2/sampsize,0,-2*sigma^2/sampsize,4*sigma^2/sampsize),nrow=3,ncol=3)

  p=tryCatch({samples = mvrnorm(n = total_n, c(beta00,beta01,beta11), sigma_pars, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)}, error=function(e){})
  if(is.null(p)){
    sigma_pars = nearPD(sigma_pars)$mat
    samples = mvrnorm(n = total_n, c(beta00,beta01,beta11), sigma_pars, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  }
  else{samples = p}
  dat_star$pred0 =  samples[,1]
  dat_star$pred1 = samples[,2] +  samples[,3]*dat_star$X

  pip_exp_emp = mean((dat_star$pred1-dat_star$y)^2 < (dat_star$pred0-dat_star$y)^2) +0.5*mean((dat_star$pred1-dat_star$y)^2 == (dat_star$pred0-dat_star$y)^2)

  return(list("Cond_Emp" = pip_emp,"Exp_Emp" = pip_exp_emp))
}

## Conditional PIPs

estimated_conditional_PIP = function(dataset){

  sub = dataset

  # Fit models m0 and m1

  mod1 = lm(Y ~ X, data = sub)
  mod0 = lm(Y ~ 1, data = sub)

  # Get predicted values

  sub$pred0 = mod0$fitted.values
  sub$pred1 = mod1$fitted.values

  # Check when m0 provides larger values than m1 (i.e. conditions C1 or C2 in slides ISNPS2022)
  sub$Ind =  sub$pred0>sub$pred1

  # Input for the cdf of Y* and proportion of both conditions
  sub$input = (sub$pred0+sub$pred1)/2
  input_1 = mean(subset(sub,pred0>pred1)$input)
  input_2 = mean(subset(sub,pred0<pred1)$input)
  prop_input_1 = nrow(subset(sub,pred0>pred1))/nrow(dataset)
  prop_input_2 = nrow(subset(sub,pred0<pred1))/nrow(dataset)

  # For version V2, we use the empirical distribution of Y* in both C1 and C2 condition groups
  df_input1 = ecdf(subset(sub,pred0>pred1)$Y)
  df_input2 = ecdf(subset(sub,pred0<pred1)$Y)

  # For version V1, we assume Y* to be distributed according to m1

  df_input1_norm = function(x){return(pnorm(x,coef(mod1)[1]+coef(mod1)[2]*subset(sub,pred0>pred1)$X[1],sd=sigma(mod1)))}
  df_input2_norm = function(x){return(pnorm(x,coef(mod1)[1]+coef(mod1)[2]*subset(sub,pred0<pred1)$X[1],sd=sigma(mod1)))}

  # We calculate versions V1 and V2

  PIP_theor_est_1=df_input1_norm(input_1)*prop_input_1+(1-df_input2_norm(input_2))*prop_input_2
  PIP_theor_est_2=df_input1(input_1)*prop_input_1+(1-df_input2(input_2))*prop_input_2

  # pvalue of m1
  pval = summary(mod1)$coefficients[2,4]

  return(list("V1" = PIP_theor_est_1,"V2" = PIP_theor_est_2,"pval"=pval))
}

## Expected PIP

estimated_expected_PIP = function(dataset,sim_size){
  require(mvnfast)
  require(MASS)
  require(Matrix)

  sub = dataset
  prop1 = mean(sub$X==1)
  N = nrow(sub)
  # Fit models m0 and m1

  mod1 = lm(Y ~ X, data = sub)
  mod0 = lm(Y ~ 1, data = sub)

  # Get parameter estimates

  beta01_est = summary(mod1)$coefficients[1,1]
  beta11_est = summary(mod1)$coefficients[2,1]
  sigma_est = sigma(mod1)

  beta00_est = summary(mod0)$coefficients[1,1]


  #Sigma_est <- sigma_est^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(beta11_est^2*prop1*(1-prop1))/(sigma_est^2)/N),2,2)
  Sigma_est <- sigma_est^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N),2,2)

  E1 <- matrix(nrow=sim_size/2, ncol = 2)
  class(E1) <- "numeric"
  rmvn(sim_size/2, c(0,-prop1*beta11_est), Sigma_est, A=E1)
  E2 <- matrix(nrow=sim_size/2, ncol = 2)
  class(E2) <- "numeric"
  rmvn(sim_size/2, c(0,prop1*beta11_est), Sigma_est, A=E2)

  PIP_expected_estimate =  (sum(E1[,1]^2< E1[,2]^2) + sum(E2[,1]^2< E2[,2]^2))/sim_size

  return(list("Exp"=PIP_expected_estimate))
}

## Cross validation

# Inputs: dataset, number of folds, type of models to be used
# Output: PIP comparing model m0 to model m1; MSE of models m0 and m1


PIP_K_cv = function(data,K,type,alpha=0.05){
  set.seed(1988)
  # Check if K is smaller than sample size
  if(K>nrow(data)){print("error: K should be less than or equal to n")}

  # Create K equally sized folds after randomly shuffling the data
  yourData<-data
  if(!(type%in%c("gbm","rf"))) {cvIndex <- createFolds(factor(data$x), K, returnTrain = T)}
  if((type%in%c("gbm","rf"))) {cvIndex <- createFolds(data$y, K, returnTrain = T)}
  
  #Perform K fold cross validation

  pip_cv = c()
  mse0_CV_sub=c()
  mse1_CV_sub=c()

  for(j in names(cvIndex)){
    trainData = yourData[cvIndex[[j]],]
    testData = yourData[-cvIndex[[j]],]

    if(type=="gaussian"){
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
    }

    if(type=="poisson"){
      mod1 = glm(y ~ x, poisson(link = "log"), data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, poisson(link = "log"), data=trainData, maxit = 5000)
    }

    if(type=="binomial"){
      mod1 = glm(y ~ x, family="binomial", data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, family="binomial", data=trainData, maxit = 5000)
    }
    
    
    if(type=="gbm"){
      ctrl <- trainControl(method = "none")
       grid <-  expand.grid(interaction.depth = 2, 
                       n.trees = 50, 
                       shrinkage = 0.1,
                      n.minobsinnode = 2)
      # trainControl <- trainControl(method="cv", number=10)
      # metric <- "RMSE"
      set.seed(99)
      # mod1 <- train(y ~ .
      #               , data=trainData
      #               , distribution="gaussian"
      #               , method="gbm"
      #               , trControl=trainControl
      #               , verbose=FALSE
      #               , tuneGrid=grid
      #               , metric=metric
      #               , bag.fraction=0.75
      # )      
      # 
      # mod0 <- train(y ~ x1+x3
      #               , data=trainData
      #               , distribution="gaussian"
      #               , method="gbm"
      #               , trControl=trainControl
      #               , verbose=FALSE
      #               , tuneGrid=grid
      #               , metric=metric
      #               , bag.fraction=0.75
      # )   
      mod1 <- train(y ~ ., data = trainData,
                    method = "gbm",
                    trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  grid)   
      
      mod0 <- train(y ~ x1+x2+x3, data = trainData,
                    method = "gbm",
                    trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  grid)  
      
      
      pred0 = predict(mod0,testData,type="raw")
      pred1 = predict(mod1,testData,type="raw")
    }
    
    if(type=="rf"){
      ctrl <- trainControl(method = "none")
      
      mod1 <- train(y ~ ., data = trainData,
                    method = "rf",
                    ntree = 100,
                    trControl = ctrl,
                    tuneGrid = data.frame(mtry = 2))    
      mod0 <- train(y ~ x1+x3, data = trainData,
                    method = "rf",
                    ntree = 100,
                    trControl = ctrl,
                    tuneGrid = data.frame(mtry = 2))    
      
      pred0 = predict(mod0,testData,type="raw")
      pred1 = predict(mod1,testData,type="raw")
    }
    
    if(!(type%in%c("gbm","rf"))){
    pred0 = predict(mod0,testData,type="response")
    pred1 = predict(mod1,testData,type="response")
    }
    if(type=="poisson") {pip_cv = c(pip_cv,mean((pred1-testData$y*log(pred1)) < (pred0-testData$y*log(pred0)))+0.5*mean((pred1-testData$y*log(pred1)) == (pred0-testData$y*log(pred0))))}
    else{    pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )}

    mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
    mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
  }

  PIP_cv = mean(pip_cv)

    mse0_CV = mean(mse0_CV_sub)
  mse1_CV = mean(mse1_CV_sub)

  return(list("PIP_cv"=PIP_cv,"mse0"=mse0_CV,"mse1"=mse1_CV))
}



PIP_K_rep_cv = function(data,K,type,alpha=0.05,reps,seed=1988){
  set.seed(seed)
  grid <-  expand.grid(interaction.depth = 2, 
                       n.trees = 50, 
                       shrinkage = 0.1,
                       n.minobsinnode = 2)
  metric <- "RMSE"
  # Check if K is smaller than sample size
  if(K>nrow(data)){print("error: K should be less than or equal to n")}
  
  # Create K equally sized folds after randomly shuffling the data
  PIP_cv = c()
  mse0_CV = c()
  mse1_CV = c()
  for(l in 1:reps){
    yourData<-data[sample(nrow(data)),]
    if(!(type%in%c("gbm","rf"))) {cvIndex <- createFolds(factor(data$x), K, returnTrain = T)}
    if((type%in%c("gbm","rf"))) {cvIndex <- createFolds(data$y, K, returnTrain = T)}
    
  #Perform K fold cross validation
  
  pip_cv = c()
  mse0_CV_sub=c()
  mse1_CV_sub=c()
  
  for(j in names(cvIndex)){
    trainData = yourData[cvIndex[[j]],]
    testData = yourData[-cvIndex[[j]],]
    
    if(type=="gaussian"){
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
    }
    
    if(type=="poisson"){
      mod1 = glm(y ~ x, poisson(link = "log"), data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, poisson(link = "log"), data=trainData, maxit = 5000)
    }
    
    if(type=="binomial"){
      mod1 = glm(y ~ x, family="binomial", data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, family="binomial", data=trainData, maxit = 5000)
    }
    
    if(type=="gbm"){
      set.seed(99)
      ctrl <- trainControl(method = "none")
      
      # mod1 <- train(y ~ .
      #               , data=trainData
      #               , distribution="gaussian"
      #               , method="gbm"
      #               , verbose=FALSE
      #               , tuneGrid=grid
      #               , metric=metric
      #               , bag.fraction=0.75
      # )      
      # 
      # mod0 <- train(y ~ x1+x3
      #               , data=trainData
      #               , distribution="gaussian"
      #               , method="gbm"
      #               , verbose=FALSE
      #               , tuneGrid=grid
      #               , metric=metric
      #               , bag.fraction=0.75
      # )   
      
      mod1 <- train(y ~ ., data = trainData,
                    method = "gbm",
                    trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  grid)   
      
      mod0 <- train(y ~ x1+x2+x3, data = trainData,
                    method = "gbm",
                    trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  grid)  
      
      pred0 = predict(mod0,testData,type="raw")
      pred1 = predict(mod1,testData,type="raw")
    }
    
    if(type=="rf"){
      ctrl <- trainControl(method = "none")
      
      mod1 <- train(y ~ ., data = trainData,
                    method = "rf",
                    ntree = 100,
                    trControl = ctrl,
                    tuneGrid = data.frame(mtry = 2))    
      mod0 <- train(y ~ x1+x3, data = trainData,
                    method = "rf",
                    ntree = 100,
                    trControl = ctrl,
                    tuneGrid = data.frame(mtry = 2))    
      
      pred0 = predict(mod0,testData,type="raw")
      pred1 = predict(mod1,testData,type="raw")
    }
    
    if(!(type%in%c("gbm","rf"))){
      pred0 = predict(mod0,testData,type="response")
      pred1 = predict(mod1,testData,type="response")
    }
    
    if(type=="poisson") {pip_cv = c(pip_cv,mean((pred1-testData$y*log(pred1)) < (pred0-testData$y*log(pred0)))+0.5*mean((pred1-testData$y*log(pred1)) == (pred0-testData$y*log(pred0))))}
    else{pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )}
    
    mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
    mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
  }
  PIP_cv = c(PIP_cv,mean(pip_cv))
  
  mse0_CV = mean(mse0_CV_sub)
  mse1_CV = mean(mse1_CV_sub)
  }
    return(list("PIP_cv"=mean(PIP_cv),"PIP_cv_lower"=quantile(PIP_cv,alpha),"PIP_cv_upper" = quantile(PIP_cv,1-alpha),"mse0"=mean(mse0_CV),"mse1"=mean(mse1_CV)))
}


## Split sample

PIP_SS = function(data){
  set.seed(1988)
  trainIndex <- createDataPartition(data$X, p = .5,
                                    list = FALSE,
                                    times = 100)

  PIP = c()

  for(j in 1:ncol(trainIndex)){
    training = data[ trainIndex[,j],]
    testing  = data[-trainIndex[,j],]

    sub=training
    mod1 = lm(Y ~ X, data = sub)
    mod0 = lm(Y ~ 1, data = sub)

    PIP = c(PIP,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2) +0.5* mean((testing$Y-predict(mod1,testing))^2==(testing$Y-predict(mod0,testing))^2))
  }

  pip_split_sample_rep100 = mean(PIP)
  pip_split_sample_rep50 = mean(PIP[1:50])
  pip_split_sample = PIP[1]

  return(list("pip_SS_rep100" = pip_split_sample_rep100,"pip_SS_rep50" = pip_split_sample_rep50,"pip_SS" = pip_split_sample))
}


## Simulation function

do_SIM = function(i,beta01,beta11,sigma,prop1,sampsize){
  set.seed(i)

  # True value for m0 intercept

  beta00 = beta01 + beta11*prop1

  # Sample dataset with specified sample size
  N = sampsize
  X = c(rep(0,(1-prop1)*sampsize),rep(1,prop1*sampsize))
  U <- rnorm(sampsize, sd = sigma)
  Y = beta01 + beta11*X + U

  training=data.frame(X, Y)

  # Conditional PIP
  cond_pip = estimated_conditional_PIP(training)

  PIP_theoretical_plugin = cond_pip$V1
  pval = cond_pip$pval

  pip_full = cond_pip$V2
  pip_cond_check = cond_pip$V1

  # Expected PIP

  PIP_expected_estimate =  estimated_expected_PIP(training, 1000000)$Exp

  # Empirical versions

  emp_out = Emp_PIP(training,1000000,prop1,beta00,beta01,beta11,sigma)

  pip_emp = emp_out$Cond_Emp
  pip_exp_emp = emp_out$Exp_Emp

  #Non-parametric approaches

  # Split-sample

  pip_split_sample = PIP_SS(training)

  # Leave-one-out CV
  CV_input=training
  names(CV_input) = c("x","y")
  PIP_cv = PIP_K_cv(CV_input,nrow(CV_input),type="gaussian",alpha=0.05)$PIP_cv

  # 5-fold CV
  out_CV_5 = PIP_K_cv(CV_input,5,type="gaussian",alpha=0.05)

  PIP_cv_5 = out_CV_5$PIP_cv
  MSE0_CV = out_CV_5$mse0
  MSE1_CV = out_CV_5$mse1

  
  # Repeated 5-fold CV
  
  out_rep_CV_5 = PIP_K_rep_cv(CV_input,5,type="gaussian",alpha=0.05,reps=100)
  
  PIP_rep_cv_5 = out_rep_CV_5$PIP_cv
  MSE0_rep_CV = out_rep_CV_5$mse0
  MSE1_rep_CV = out_rep_CV_5$mse1
  
  
  
  return(list("PIP_cond" = PIP_theoretical_plugin,
              "PIP_exp" = PIP_expected_estimate,
              "Emp_Cond"  = pip_emp,
              "Emp_Exp"  = pip_exp_emp,
              "PIP_SS" = pip_split_sample$pip_SS,
              "PIP_SS_rep50" = pip_split_sample$pip_SS_rep50,
              "PIP_SS_rep100" = pip_split_sample$pip_SS_rep100,
              "PIP_LOO" = PIP_cv,
              "pval_mod1" = pval,
              "pip_full" =    pip_full,
              "pip_cond_check" = pip_cond_check,
              "PIP_CV5" = PIP_cv_5,
              "MSE0_CV" = MSE0_CV,
              "MSE1_CV" = MSE1_CV,
              "PIP_rep_CV5" = PIP_rep_cv_5,
              "MSE0_rep_CV" = MSE0_rep_CV,
              "MSE1_rep_CV" = MSE1_rep_CV
              
  ))
}


# Relation between pvalue and PIP_C1

f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-2))
  )
}



# Assess coverage of repeated CV

do_SIM_coverage = function(i,beta01,beta11,sigma,prop1,sampsize){
  set.seed(i)
  # True value for m0 intercept
  beta00 = beta01 + beta11*prop1
  # Sample dataset with specified sample size
  N = sampsize
  X = c(rep(0,(1-prop1)*sampsize),rep(1,prop1*sampsize))
  U <- rnorm(sampsize, sd = sigma)
  Y = beta01 + beta11*X + U
  training=data.frame(X, Y)
  CV_input=training
  names(CV_input) = c("x","y")
  # Repeated 5-fold CV
  
  out_rep_CV_5 = PIP_K_rep_cv(CV_input,5,type="gaussian",alpha=0.05,reps=100)
  
  PIP = out_rep_CV_5$PIP_cv
  PIP_lower = out_rep_CV_5$PIP_cv_lower
  PIP_upper = out_rep_CV_5$PIP_cv_upper
  
  return(list("PIP" = PIP,
              "PIP_lower" = PIP_lower,
              "PIP_upper"  = PIP_upper
  ))
}



do_SIM_GBM = function(i, sampsize){
  set.seed(i)

  trainControl <- trainControl(method="cv", number=10)
  metric <- "RMSE"
  
  gbmGrid <-  expand.grid(interaction.depth = 2, 
                          n.trees = 50, 
                          shrinkage = 0.1,
                          n.minobsinnode = 2)
  
  x1=rnorm(sampsize,5,0.9)
  x2=rnorm(sampsize,10,2.5)
  x3=rnorm(sampsize,20,6.4)

    
  error=rnorm(sampsize,0,1.6)
  y=15-(4*x1)+(2.5*x2)+(5*x3)+error
  
  dat = data.frame(cbind(y,x1,x2,x3))
  
  set.seed(99)
  mod1 <- train(y ~ .
                , data=dat
                , distribution="gaussian"
                , method="gbm"
                , verbose=FALSE
                , tuneGrid=gbmGrid
                , metric=metric
                , bag.fraction=0.75
  )                  
  
  mod0 <- train(y ~ x1+x3
                , data=dat
                , distribution="gaussian"
                , method="gbm"
                , verbose=FALSE
                , tuneGrid=gbmGrid
                , metric=metric
                , bag.fraction=0.75
  )                  
  

  
  # Empirical
  total_n = 1000000
  x1_star=rnorm(total_n,5,0.9)
  x2_star=rnorm(total_n,10,2.5)
  x3_star=rnorm(total_n,20,6.4)

  error_star=rnorm(total_n,0,1.6)
  y_star=15-(4*x1_star)+(2.5*x2_star)+(5*x3_star)+error_star
  
  dat_star = data.frame(cbind("x1"=x1_star,"x2"=x2_star,"x3"=x3_star))
  
  pred0_star = predict(mod0, newdata=dat_star, type="raw")
  pred1_star = predict(mod1, newdata=dat_star, type="raw")
  
  pip_emp = mean((pred1_star-y_star)^2 < (pred0_star-y_star)^2) + 0.5*mean((pred1_star-y_star)^2 == (pred0_star-y_star)^2)

  CV5 = PIP_K_cv(dat,5,"gbm",alpha=0.05)
  
  PIP = CV5$PIP_cv 
  mse1_CV = CV5$MSE1 
  mse0_CV = CV5$MSE0 

  
  rep_CV5 = PIP_K_rep_cv(dat,5,"gbm",alpha=0.05,reps=20)

  #Split sample
  select = sample(nrow(dat),nrow(dat)/2)
  dat_train = dat[select,]
  dat_test = dat[-select,]
  
  set.seed(99)
  mod1 <- train(y ~ .
                , data=dat_train
                , distribution="gaussian"
                , method="gbm"
                , verbose=FALSE
                , tuneGrid=gbmGrid
                , metric=metric
                , bag.fraction=0.75
  )                  
  
  mod0 <- train(y ~ x1+x3
                , data=dat_train
                , distribution="gaussian"
                , method="gbm"
                , verbose=FALSE
                , tuneGrid=gbmGrid
                , metric=metric
                , bag.fraction=0.75
  )   
  
  pred0 = predict(mod0, newdata=dat_test, type="raw")
  pred1 = predict(mod1, newdata=dat_test, type="raw")
  
  pip_SS = mean((pred1-dat_test$y)^2 < (pred0-dat_test$y)^2) + 0.5*mean((pred1-dat_test$y)^2 == (pred0-dat_test$y)^2)
  
  
  return(list("PIP_SS"=pip_SS,"PIP_CV5" = PIP,"PIP_rep_CV5" = rep_CV5$PIP_cv,
              "PIP_emp" = pip_emp))
}






do_SIM_GBM_new = function(i, sampsize){
  set.seed(i)

  x1=round(runif(sampsize,0,6))
  x2=round(rnorm(sampsize))
  x3=rbinom(sampsize,1,prob=0.5)
  x4=round(runif(sampsize,0,1),1)
  x5=round(runif(sampsize,1,2),1)
  
  error=rnorm(sampsize,0,1.6)
  y=(4*abs(x1))^(3*x4)+(5*x2)+(2*x3)^x5+error

  dat = data.frame(cbind(y,x1,x2,x3,x4,x5))

  set.seed(99)
  ctrl <- trainControl(method = "none")
  gbmGrid <-  expand.grid(interaction.depth = 2, 
                          n.trees = 50, 
                          shrinkage = 0.1,
                          n.minobsinnode = 2)
  
  mod1 <- train(y ~ ., data = dat,
                method = "gbm",
                trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  gbmGrid )   
  
  mod0 <- train(y ~ x1+x2+x3, data = dat,
                method = "gbm",
                trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  gbmGrid)   

  
  # Empirical
  total_n = 1000000
  
  x1_star=round(runif(total_n,0,6))
  x2_star=round(rnorm(total_n))
  x3_star=rbinom(total_n,1,prob=0.5)
  x4_star=round(runif(total_n,0,1),1)
  x5_star=round(runif(total_n,1,2),1)

  error_star=rnorm(total_n,0,1.6)
  y_star=(4*abs(x1_star))^(3*x4_star)+(5*x2_star)+(2*x3_star)^x5_star+error_star
  

  dat_star = data.frame(cbind("x1"=x1_star,"x2"=x2_star,"x3"=x3_star,"x4"=x4_star,"x5"=x5_star))
  
  pred0_star = predict(mod0, newdata=dat_star, type="raw")
  pred1_star = predict(mod1, newdata=dat_star, type="raw")
  
  pip_emp = mean((pred1_star-y_star)^2 < (pred0_star-y_star)^2) + 0.5*mean((pred1_star-y_star)^2 == (pred0_star-y_star)^2)
  
  CV5 = PIP_K_cv(dat,5,"gbm",alpha=0.05)
  
  PIP = CV5$PIP_cv 
  mse1_CV = CV5$mse1 
  mse0_CV = CV5$mse0 
  

  rep_CV5 = PIP_K_rep_cv(dat,5,"gbm",alpha=0.05,reps=10)
  
  #Split sample
  select = sample(nrow(dat),nrow(dat)/2)
  dat_train = dat[select,]
  dat_test = dat[-select,]
  
  set.seed(99)
  mod1 <- train(y ~ ., data = dat_train,
                method = "gbm",
                trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  gbmGrid)   
  
  mod0 <- train(y ~ x1+x2+x3, data = dat_train,
                method = "gbm",
                trControl = ctrl,verbose=F, distribution="gaussian",tuneGrid=  gbmGrid)   
  
  
  

  
  
  pred0 = predict(mod0, newdata=dat_test, type="raw")
  pred1 = predict(mod1, newdata=dat_test, type="raw")
  
  pip_SS = mean((pred1-dat_test$y)^2 < (pred0-dat_test$y)^2) + 0.5*mean((pred1-dat_test$y)^2 == (pred0-dat_test$y)^2)
  
  
  return(list("PIP_SS"=pip_SS,"PIP_CV5" = PIP,"PIP_rep_CV5" = rep_CV5$PIP_cv,
              "PIP_emp" = pip_emp,"mse1_CV"=mse1_CV, "mse0_CV"=mse0_CV,
              "mse1_rep_CV"=rep_CV5$mse1, "mse0_rep_CV"=rep_CV5$mse0))
}













do_SIM_RF_fasterSS = function(i, sampsize){
  set.seed(i)
  
  x1=rnorm(sampsize,5,0.9)
  x2=rnorm(sampsize,10,2.5)
  x3=rnorm(sampsize,20,6.4)
  x4=rnorm(sampsize,15,3)
  x5=rnorm(sampsize,5,1)
  
  
  error=rnorm(sampsize,0,1.6)
  y=15-(4*x1)+(2.5*x2)+(5*x3)+(6*x4)+(10*x5)+error
  
  dat = data.frame(cbind(y,x1,x2,x3,x4,x5))
  
  set.seed(99)
  ctrl <- trainControl(method = "none")
  
  mod11 <- train(y ~ ., data = dat,
                method = "rf",
                ntree = 10,
                trControl = ctrl,
                tuneGrid = data.frame(mtry = 2))    
  mod01 <- train(y ~ x1+x3, data = dat,
                method = "rf",
                ntree = 10,
                trControl = ctrl,
                tuneGrid = data.frame(mtry = 2))    
  
  
  #CV5 = PIP_K_cv(dat,5,"rf",alpha=0.05)
  
  #PIP = CV5$PIP_cv 
  #mse1_CV = CV5$MSE1 
  #mse0_CV = CV5$MSE0 
  
  
  #rep_CV5 = PIP_K_rep_cv(dat,5,"rf",alpha=0.05,reps=20)
  
  #Split sample
  # select = sample(nrow(dat),nrow(dat)/2)
  # dat_train = dat[select,]
  # dat_test = dat[-select,]
  # 
  # set.seed(99)
  # mod1 <- train(y ~ ., data = dat_train,
  #               method = "rf",
  #               ntree = 10,
  #               trControl = ctrl,
  #               tuneGrid = data.frame(mtry = 2))         
  # 
  # mod0 <- train(y ~ x1+x3, data = dat_train,
  #               method = "rf",
  #               ntree = 10,
  #               trControl = ctrl,
  #               tuneGrid = data.frame(mtry = 2))    
  # 
  # pred0 = predict(mod0, newdata=dat_test, type="raw")
  # pred1 = predict(mod1, newdata=dat_test, type="raw")
  
  #pip_SS = mean((pred1-dat_test$y)^2 < (pred0-dat_test$y)^2) + 0.5*mean((pred1-dat_test$y)^2 == (pred0-dat_test$y)^2)
  
  
  # Empirical
  set.seed(1988)
  total_n = 1000000
  x1_star=rnorm(total_n,5,0.9)
  x2_star=rnorm(total_n,10,2.5)
  x3_star=rnorm(total_n,20,6.4)
  x4_star=rnorm(total_n,15,3)
  x5_star=rnorm(total_n,5,1)
  
  error_star=rnorm(total_n,0,1.6)
  y_star=15-(4*x1_star)+(2.5*x2_star)+(5*x3_star)+(6*x4_star)+(10*x5_star)+error_star
  
  dat_star = data.frame(cbind("x1"=x1_star,"x2"=x2_star,"x3"=x3_star,"x4"=x4_star,"x5"=x5_star))
  
  pred0_star = predict(mod01, newdata=dat_star, type="raw")
  pred1_star = predict(mod11, newdata=dat_star, type="raw")
  
  pip_emp = mean((pred1_star-y_star)^2 < (pred0_star-y_star)^2) + 0.5*mean((pred1_star-y_star)^2 == (pred0_star-y_star)^2)
  
  
  return(list("PIP_SS"=0.5,"PIP_CV5" = 0.5,"PIP_rep_CV5" = 0.5,
              "PIP_emp" = pip_emp))
}


