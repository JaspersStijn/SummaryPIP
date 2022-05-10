# Function to calculate the K-fold cross-validation PIP
#
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

  Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(prop1*(1-prop1))/(sigma^2)/N),2,2)
  E1 <- matrix(nrow=10000000/2, ncol = 2)
  class(E1) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,-0.5*beta11), Sigma, A=E1)
  E2 <- matrix(nrow=10000000/2, ncol = 2)
  class(E2) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
  rmvn(10000000/2, c(0,0.5*beta11), Sigma, A=E2)

  PIP_expected_true =  (sum(E1[,1]^2< E1[,2]^2) + sum(E2[,1]^2< E2[,2]^2))/10000000

  N = 100000000
  Sigma <- sigma^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(prop1*(1-prop1))/(sigma^2)/N),2,2)
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

Emp_PIP = function(dataset,total_n,prop1,beta01,beta11,sigma){
  require(mvnfast)
  require(Matrix)
  require(MASS)
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

  sigma_pars = matrix(c(var_covar0,sigma(mod1)^2/N,0,sigma(mod1)^2/N,var_b01,covar_mod1,0,covar_mod1,var_b11),nrow=3,ncol=3)

  p=tryCatch({samples = mvrnorm(n = total_n, c(coef_mod0[1],coef_mod1[1],coef_mod1[2]), sigma_pars, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)}, error=function(e){})
  if(is.null(p)){
    sigma_pars = nearPD(sigma_pars)$mat
    samples = mvrnorm(n = total_n, c(coef_mod0[1],coef_mod1[1],coef_mod1[2]), sigma_pars, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
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


  Sigma_est <- sigma_est^2*matrix(c(1+2/N,1+1/N,1+1/N,1+1/N+(beta11_est^2*prop1*(1-prop1))/(sigma_est^2)/N),2,2)

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
  set.seed(123)
  # Check if K is smaller than sample size
  if(K>nrow(data)){print("error: K should be less than or equal to n")}

  # Create K equally sized folds after randomly shuffling the data
  yourData<-data[sample(nrow(data)),]
  folds <- cut(seq(1,nrow(yourData)),breaks=K,labels=FALSE)

  #Perform K fold cross validation

  pip_cv = c()
  mse0_CV_sub=c()
  mse1_CV_sub=c()

  for(i in 1:K){

    #Segment your data by fold using the which() function

    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourData[testIndexes, ]
    trainData <- yourData[-testIndexes, ]

    #Use the test and train data partitions based on model type

    # if(type=="nonlinear"){
    #   mod1 = nls(y~p1+p2*x**p3, data=trainData,
    #              start = list(p1=3,
    #                           p2= 3,
    #                           p3= 1))
    #   mod0 = nls(y~p1+p2*x, data=trainData,
    #              start = list(p1=2,
    #                           p2= 1))
    # }


    if(type=="gaussian"){
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
    }

    if(type=="poisson"){
      mod1 = glm(y ~ x, poisson(link = "log"), data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, poisson(link = "log"), data=trainData, maxit = 5000)
    }

    # if(type=="gamma"){
    #   mod1 = glm(y ~ x, Gamma(link = "log"), data=trainData, maxit = 5000,start=c(beta0+0.05,beta1-0.5))
    #   mod0 = glm(y ~ 1, Gamma(link = "log"), data=trainData, maxit = 5000)
    # }

    if(type=="binomial"){
      mod1 = glm(y ~ x, family="binomial", data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, family="binomial", data=trainData, maxit = 5000)
    }

    pred0 = predict(mod0,testData,type="response")
    pred1 = predict(mod1,testData,type="response")

    if(type=="poisson") {pip_cv = c(pip_cv,mean((pred1-testData$y*log(pred1)) < (pred0-testData$y*log(pred0)))+0.5*mean((pred1-testData$y*log(pred1)) == (pred0-testData$y*log(pred0))))}
    else{    pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2) )}

    mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
    mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
  }

  PIP_cv = mean(pip_cv)

  # PIP_cv_lower = PIP_cv - qnorm(1-alpha/2)*sqrt(PIP_cv*(1-PIP_cv)/K)
  # PIP_cv_upper = PIP_cv + qnorm(1-alpha/2)*sqrt(PIP_cv*(1-PIP_cv)/K)

  # PIP_cv_mean = mean(pip_cv)
  # PIP_cv_var = var(pip_cv)

  # PIP_cv_lower2 = expit(logit( PIP_cv_mean) - qnorm(1-alpha/2)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  # PIP_cv_upper2 = expit(logit( PIP_cv_mean) + qnorm(1-alpha/2)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))

  # PIP_cv_lower2 = expit(logit( PIP_cv_mean) - qt(1-alpha/2,K-1)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  # PIP_cv_upper2 = expit(logit( PIP_cv_mean) + qt(1-alpha/2,K-1)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))

  mse0_CV = mean(mse0_CV_sub)
  mse1_CV = mean(mse1_CV_sub)

  return(list("PIP_cv"=PIP_cv,"mse0"=mse0_CV,"mse1"=mse1_CV))
}


## Split sample

PIP_SS = function(data){
  N=nrow(data)
  smp_size <- floor(0.5*N)
  train_ind <- sample(seq_len(N), size = smp_size)
  training_sub <-data[train_ind, ]
  testing <- data[-train_ind, ]

  mod1 = lm(Y ~ X, data = training_sub)
  mod0 = lm(Y ~ 1, data = training_sub)

  pred0 = predict(mod0,testing)
  pred1 = predict(mod1,testing)

  pip_split_sample = mean((pred1-testing$Y)^2 < (pred0-testing$Y)^2) +0.5*mean((pred1-testing$Y)^2 == (pred0-testing$Y)^2)

  return(pip_split_sample)
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

  emp_out = Emp_PIP(training,1000000,prop1,beta01,beta11,sigma)

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

  return(list("PIP_cond" = PIP_theoretical_plugin,
              "PIP_exp" = PIP_expected_estimate,
              "Emp_Cond"  = pip_emp,
              "Emp_Exp"  = pip_exp_emp,
              "PIP_SS" = pip_split_sample,
              "PIP_LOO" = PIP_cv,
              "pval_mod1" = pval,
              "pip_full" =    pip_full,
              "pip_cond_check" = pip_cond_check,
              "PIP_CV5" = PIP_cv_5,
              "MSE0_CV" = MSE0_CV,
              "MSE1_CV" = MSE1_CV
  ))
}

