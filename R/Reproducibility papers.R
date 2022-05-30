# Papers discussed in "Evaluating the replicability of social science
# experiments in Nature and Science between
# 2010 and 2015", Camerer et al. (2018)
source("R/TwoSample_Functions.R")


summary_pips = data.frame()

# Paper 1: Ackerman et al. (2010)

set.seed(62)
n1 = 26; mu1= 5.80; sd1 = 0.76
n2 = 28; mu2= 5.38; sd2 = 0.79

x= c(rep(0,n1),rep(1,n2))
y = c(rnorm(n1,mu1,sd1),rnorm(n2,mu2,sd2))

dat1 = data.frame(cbind(x,y))

mod = lm(y~x,data=dat1)
summary(mod)
summary(mod)$coef[2,4]

# pip_out = PIP_K_cv(dat1,5,"gaussian")
# pip_out$PIP_cv
# pip_out$PIP_cv_lower
# pip_out$PIP_cv_upper
# 
# summary_pips = rbind(summary_pips,data.frame(Study = "Ackerman et al. (2010)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = pip_out$PIP_cv_lower,PIP = pip_out$PIP_cv,PIP_upper = pip_out$PIP_cv_upper,Replicated="No",PredictionMarket = 0.15, Survey = 0.13))


pip_sub = c()
pip_sub_lower = c()
pip_sub_upper = c()

for(i in 1:10000){
  fit = PIP_K_rep_cv(dat1,5,"gaussian",0.05,100,seed=i)
  pip_sub = c(pip_sub,fit$PIP_cv)
  pip_sub_lower = c(pip_sub_lower,fit$PIP_cv_lower)
  pip_sub_upper = c(pip_sub_upper,fit$PIP_cv_upper)
  
  print(i)
}

summary_pips = rbind(summary_pips,data.frame(Study = "Ackerman et al. (2010)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = quantile(pip_sub,0.025),PIP =mean(pip_sub),PIP_upper = quantile(pip_sub,0.975),Replicated="No",PredictionMarket = 0.15, Survey = 0.13))



# Paper 2: Wilson et al. (2014)

set.seed(62)
n1 = 15; mu1= 3.20; sd1 = 2.23
n2 = 15; mu2= 6.87; sd2 = 1.91

x= c(rep(0,n1),rep(1,n2))
y = c(rnorm(n1,mu1,sd1),rnorm(n2,mu2,sd2))

dat1 = data.frame(cbind(x,y))

mod = lm(y~x,data=dat1)
summary(mod)

# pip_out = PIP_K_cv(dat1,5,"gaussian")
# pip_out$PIP_cv
# pip_out$PIP_cv_lower
# pip_out$PIP_cv_upper
# 
# summary_pips = rbind(summary_pips,data.frame(Study = "Wilson et al. (2014)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = pip_out$PIP_cv_lower,PIP = pip_out$PIP_cv,PIP_upper = pip_out$PIP_cv_upper,Replicated="Yes",PredictionMarket = 0.46, Survey = 0.52))

pip_sub = c()
pip_sub_lower = c()
pip_sub_upper = c()

for(i in 1:10000){
  fit = PIP_K_rep_cv(dat1,5,"gaussian",0.05,100,seed=i)
  pip_sub = c(pip_sub,fit$PIP_cv)
  pip_sub_lower = c(pip_sub_lower,fit$PIP_cv_lower)
  pip_sub_upper = c(pip_sub_upper,fit$PIP_cv_upper)
  
  print(i)
}

summary_pips = rbind(summary_pips,data.frame(Study = "Wilson et al. (2014)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = quantile(pip_sub,0.025),PIP =mean(pip_sub),PIP_upper = quantile(pip_sub,0.975),Replicated="Yes",PredictionMarket = 0.46, Survey = 0.52))


# Paper 3: Gervais et al. (2012)
set.seed(457)
n1 = 31; mu1= 61.55; sd1 = 35.68
n2 = 26; mu2= 41.42; sd2 = 31.47

x= c(rep(0,n1),rep(1,n2))
y = c(rnorm(n1,mu1,sd1),rnorm(n2,mu2,sd2))

dat1 = data.frame(cbind(x,y))

mod = lm(y~x,data=dat1)
summary(mod)

# pip_out = PIP_K_cv(dat1,5,"gaussian")
# pip_out$PIP_cv
# pip_out$PIP_cv_lower
# pip_out$PIP_cv_upper

# summary_pips = rbind(summary_pips,data.frame(Study = "Gervais et al. (2012)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = pip_out$PIP_cv_lower,PIP = pip_out$PIP_cv,PIP_upper = pip_out$PIP_cv_upper,Replicated="No",PredictionMarket = 0.17, Survey = 0.20))

pip_sub = c()
pip_sub_lower = c()
pip_sub_upper = c()

for(i in 1:10000){
  set.seed(i)
  fit = PIP_K_rep_cv(dat1,5,"gaussian",0.05,100,seed=i)
  pip_sub = c(pip_sub,fit$PIP_cv)
  pip_sub_lower = c(pip_sub_lower,fit$PIP_cv_lower)
  pip_sub_upper = c(pip_sub_upper,fit$PIP_cv_upper)
  
  print(i)
}

summary_pips = rbind(summary_pips,data.frame(Study = "Gervais et al. (2012)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = quantile(pip_sub,0.025),PIP =mean(pip_sub),PIP_upper = quantile(pip_sub,0.975),Replicated="No",PredictionMarket = 0.17, Survey = 0.20))

# Paper 4: Balafoutas et al. (2012)

set.seed(62)
n1 = 36; prop1= 0.306
n2 = 36; prop2= 0.583

x= c(rep(0,n1),rep(1,n2))
y = c(rep(1,11),rep(0,25),rep(1,21),rep(0,15))
cor(x,y)
dat1 = data.frame(cbind(factor(x),factor(y)))
chisq <- chisq.test(dat1$X1,dat1$X2,correct=FALSE)
chisq$observed
chisq$expected
chisq$statistic

dat1 = data.frame(cbind(x,y))

mod = glm(y~x,data=dat1,family = "binomial")
summary(mod)

# pip_out = PIP_K_cv(dat1,5,"binomial")
# pip_out$PIP_cv
# pip_out$PIP_cv_lower
# pip_out$PIP_cv_upper
# 
# summary_pips = rbind(summary_pips,data.frame(Study = "Balafoutas et al.  (2012)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = pip_out$PIP_cv_lower,PIP = pip_out$PIP_cv,PIP_upper = pip_out$PIP_cv_upper,Replicated="Yes",PredictionMarket = 0.75, Survey = 0.43))

pip_sub = c()
pip_sub_lower = c()
pip_sub_upper = c()

for(i in 1:10000){
  fit = PIP_K_rep_cv(dat1,5,"binomial",0.05,100,seed=i)
  pip_sub = c(pip_sub,fit$PIP_cv)
  pip_sub_lower = c(pip_sub_lower,fit$PIP_cv_lower)
  pip_sub_upper = c(pip_sub_upper,fit$PIP_cv_upper)
  
  print(i)
}

summary_pips = rbind(summary_pips,data.frame(Study = "Balafoutas et al.  (2012)",pval = round(summary(mod)$coef[2,4],4),PIP_lower = quantile(pip_sub,0.025),PIP =mean(pip_sub),PIP_upper = quantile(pip_sub,0.975),Replicated="Yes",PredictionMarket = 0.75, Survey = 0.43))


summary_pips[,c("pval","PIP_lower","PIP","PIP_upper")] = round(summary_pips[,c("pval","PIP_lower","PIP","PIP_upper")],4)

summary_pips[order(summary_pips[,c("PIP")] ),c("Study","pval","Replicated","PredictionMarket","Survey","PIP_lower","PIP","PIP_upper")]

save.image("R/Reproducibility.Rdata")



