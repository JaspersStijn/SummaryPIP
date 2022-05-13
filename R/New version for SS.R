
library(caret)



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

  # For version V3, we assume Y* to be distributed according to m0

  df_input1_norm3 = function(x){return(pnorm(x,10,sd=2))}
  df_input2_norm3 = function(x){return(pnorm(x,10,sd=2))}


  # We calculate versions V1, V2 and V3

  PIP_theor_est_1=df_input1_norm(input_1)*prop_input_1+(1-df_input2_norm(input_2))*prop_input_2
  PIP_theor_est_2=df_input1(input_1)*prop_input_1+(1-df_input2(input_2))*prop_input_2
  PIP_theor_est_3=df_input1_norm3(input_1)*prop_input_1+(1-df_input2_norm3(input_2))*prop_input_2

  # pvalue of m1
  pval = summary(mod1)$coefficients[2,4]

  return(list("V1" = PIP_theor_est_1,"V2" = PIP_theor_est_2,"pval"=pval,"V3" =PIP_theor_est_3 ))
}

do_SIM = function(i,beta11,sampsize){
  set.seed(i)
  beta01=10

  sigma=2
  prop1 = 0.5

  # True value for m0 intercept

  beta00 = beta01 + beta11*prop1

  # Sample dataset with specified sample size
  N = sampsize
  X = c(rep(0,(1-prop1)*sampsize),rep(1,prop1*sampsize))
  U <- rnorm(sampsize, sd = sigma)
  Y = beta01 + beta11*X + U

  dat = data.frame(X, Y)

  # Conditional PIPs
  conds = estimated_conditional_PIP(dat)

  set.seed(1988)
  trainIndex <- createDataPartition(dat$X, p = .5,
                                    list = FALSE,
                                    times = 100)

  PIP = c()

  for(j in 1:ncol(trainIndex)){
      training = dat[ trainIndex[,j],]
      testing  = dat[-trainIndex[,j],]

      sub=training
      mod1 = lm(Y ~ X, data = sub)
      mod0 = lm(Y ~ 1, data = sub)

      PIP = c(PIP,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2))
  }

    cvIndex <- createFolds(factor(dat$X), 10, returnTrain = T)
    PIP_sub10 =c()
  for(j in names(cvIndex)){
    training = dat[cvIndex[[j]],]
    testing = dat[-cvIndex[[j]],]
    sub=training
    mod1 = lm(Y ~ X, data = sub)
    mod0 = lm(Y ~ 1, data = sub)
    PIP_sub10 = c(PIP_sub10,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2))
  }

    cvIndex <- createFolds(factor(dat$X), 5, returnTrain = T)
    PIP_sub5 =c()
    for(j in names(cvIndex)){
      training = dat[cvIndex[[j]],]
      testing = dat[-cvIndex[[j]],]
      sub=training
      mod1 = lm(Y ~ X, data = sub)
      mod0 = lm(Y ~ 1, data = sub)
      PIP_sub5 = c(PIP_sub5,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2))
    }

    PIP_sub5 =c()
    for(j in names(cvIndex)){
      training = dat[cvIndex[[j]],]
      testing = dat[-cvIndex[[j]],]
      sub=training
      mod1 = lm(Y ~ X, data = sub)
      mod0 = lm(Y ~ 1, data = sub)
      PIP_sub5 = c(PIP_sub5,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2))
    }

  return(list("V1" = conds$V1,"V2" = conds$V2,"V3" = conds$V3,"PIP_SS_rep100" = mean(PIP),"PIP_SS_rep50" = mean(PIP[1:50]),"PIP_SS" = PIP[sample(100,1)],"PIP_CV5"=  mean(PIP_sub5),"PIP_CV10"=  mean(PIP_sub10),"pval_mod1"= summary(mod1)$coefficients[2,4]))
}

require(doSNOW)

for(beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
cl <- makeCluster(7)
registerDoSNOW(cl)
iterations <- 10000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
output <- foreach(i=1:iterations,.packages=c("caret"),
                  .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                    result <- do_SIM(i,beta11,sampsize)
                    pip_SS_rep100 = result$PIP_SS_rep100
                    pip_SS_rep50 = result$PIP_SS_rep50
                    pip_SS_corrected = result$PIP_SS
                    pip_CV10 = result$PIP_CV10
                    pip_CV5 = result$PIP_CV5
                    pip_C1=result$V1
                    pip_C2=result$V2
                    pip_C3=result$V3
                    pval = result$pval_mod1
                    return(cbind(pip_C1,pip_C2,pip_C3,pip_SS_rep100,pip_SS_rep50,pip_SS_corrected,pip_CV5,pip_CV10,pval))
                  }
close(pb)
stopCluster(cl)

save(output,file=paste0(paste(paste0("R/Output_Sims/sim_updatedSS_effect",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))

new_output = output

use = load(paste0(paste(paste0("R/Output_Sims/sim_effect",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
output = get(use)
colnames(output)[1] = "pip_C1"
colnames(output)[8] = "pip_C2"
colnames(output)
if(beta11==0){comp = cbind(output[,c("emp_cond")],new_output[,-ncol(new_output)],output[,c("pip_LOO")])}
if(beta11!=0){comp = cbind(output[,c("emp_cond")],new_output[,-c(3,ncol(new_output))],output[,c("pip_LOO")])}
colnames(comp)[ncol(comp)] = "pip_LOO"
colnames(comp)[1] = "emp_cond"

boxplot(comp)
abline(h=mean(output[,c("emp_cond")]),col="red",lwd=2)
abline(h=mean(output[,c("emp_exp")]),col="blue",lwd=2)
means = apply(comp,2,mean)
points(means,pch=20)

require(mvnfast)
true_vals = TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)

abline(h=true_vals$PIP_theor,lwd=2,lty=2)
abline(h=true_vals$PIP_exp,lwd=2,lty=3)


legend("bottomleft",c("Emp_cond","Emp_Exp","Theoretical","True_Exp"),col=c("red","blue","black","black"),lty=c(1,1,2,3),lwd=2)
title(paste(paste0("Beta11= ",beta11),paste0("Sampsize= ",sampsize),sep="   "))
}
}


plot(output[,"pval_mod1"],output[,"pip_SS"])
points(new_output[,"pval"],new_output[,"pip_SS_rep100"],col="red")
points(output[,"pval_mod1"],output[,"pip_C1"],col="blue")

pvals = new_output[,"pval"]

f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-1))
  )
}
lines(seq(0,max(max(pvals),4e-10),length=1000),sapply(seq(0,max(max(pvals),4e-10),length=1000),f_pip,n=sampsize),col="blue",lwd=2)
points(output[,"pval_mod1"],output[,"pip_exp"],col="darkgreen")

plot(new_output[,"pval"],new_output[,"pip_CV5"],col="black")














