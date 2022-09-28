# Function used

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
  
  return(list("The" = The,"The_emp" = The_emp,"MSE_diff"=mse_diff,"MSE_diff_emp"=true_mse_diff,"pval"=pval,"beta11_hat"=coef(mod)[2] ))
}


## Focus on out-sample empirical PIP below
####  
# clear difference between in-sample/out-sample for seed = 17 
i=17
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

# Fit two models
mod = lm(Y~X,data=training)
mod0 = lm(Y~1,data=training)

# New observations
X_star = c(rep(0,(1-prop1)*1000000),rep(1,prop1*1000000))
Y_star = beta01 + beta11*X_star + rnorm(1000000, sd = sigma)
dat_star = data.frame("X"=X_star,"y"=Y_star)

pred1 = predict(mod,dat_star)
pred0 = predict(mod0,dat_star)
true_mse_diff = mean((dat_star$y-pred1)^2) - mean((dat_star$y-pred0)^2) 
pval = coef(summary(mod))[2,4]
The_emp = mean(((dat_star$y-pred1)^2) < ((dat_star$y-pred0)^2)) + 0.5*mean(((dat_star$y-pred1)^2) == ((dat_star$y-pred0)^2)) 


