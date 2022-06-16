# Application of the PIP for the two-sample case

source("R/TwoSample_Functions.R")

require(doSNOW)
require(ggplot2)

for (beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
    cl <- makeCluster(7)
    registerDoSNOW(cl)
    iterations <- 10000
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    output <- foreach(i=1:iterations,.packages=c("MASS","Matrix","mvnfast","caret"),
                      .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                        result <- do_SIM(i,10,beta11,2,0.5,sampsize)
                        pip_cond = result$PIP_cond
                        pip_exp = result$PIP_exp
                        emp_cond<-result$Emp_Cond
                        emp_exp<-result$Emp_Exp
                        pip_SS = result$PIP_SS
                        pip_SS_rep50 = result$PIP_SS_rep50
                        pip_SS_rep100 = result$PIP_SS_rep100
                        pip_LOO = result$PIP_LOO
                        pval_mod1 = result$pval_mod1
                        pip_full = result$pip_full
                        pip_cond_check=result$pip_cond_check
                        pip_CV5=result$PIP_CV5
                        mse0_CV = result$MSE0_CV
                        mse1_CV = result$MSE1_CV
                        pip_rep_CV5=result$PIP_rep_CV5
                        mse0_rep_CV = result$MSE0_rep_CV
                        mse1_rep_CV = result$MSE1_rep_CV
                        return(cbind(pip_cond,pip_exp,emp_cond,emp_exp,
                                     pip_SS,pip_SS_rep50,pip_SS_rep100,pip_LOO,pval_mod1,pip_full,pip_cond_check,
                                     pip_CV5,mse0_CV,mse1_CV,pip_rep_CV5,mse0_rep_CV,mse1_rep_CV))
                      }
    close(pb)
    stopCluster(cl)

    save(output,file=paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
  }

}

boxplot(output[,c("pip_CV5","pip_rep_CV5")])


PIP = c()
pval = c()
MSE_diff = c()
EffectSize = c()
Samplesize = c()

for(beta11 in c(0,-1,-4)){
  layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(4, 4,1))
  par(mai=rep(0.5, 4))
  for( sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "C1"
    colnames(output)[2] ="Exp"
    colnames(output)[10] = "C2"
    colnames(output)[5] ="SS"
    colnames(output)[7] ="SS_rep100"
    colnames(output)[8] ="LOO"
    colnames(output)[12] ="CV5"
    colnames(output)[15] ="rep_CV5"
  
    means <- apply(output,2,mean)
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    #boxplot(output[,c("pip_cond1","pip_cond2","pip_exp")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_cond1","pip_cond2","pip_exp","pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_C1","pip_C2","emp_cond","pip_exp","emp_exp" ,"pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    boxplot(output[,c("C1","C2","Exp","LOO","CV5","rep_CV5","SS","SS_rep100")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    abline(h=true_vals$PIP_theor,col="black",lty=2)
    abline(h=true_vals$PIP_exp ,col='black',lty=3)
    #abline(h=true_vals$PIP_exp_limit ,col='darkgreen',lty=2,lwd=2)
    #points(means[c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_SS","pip_CV5")],pch=25)
    points(means[c("C1","C2","Exp","LOO","CV5","rep_CV5","SS","SS_rep100")],pch=20)

    print(means)
    PIP = c(PIP,output[,c("rep_CV5")])
    pval = c(pval,output[,c("pval_mod1")])
    MSE_diff = c(MSE_diff,output[,c("mse1_CV")]-output[,c("mse0_CV")])
    EffectSize = c(EffectSize,rep(paste0('Effect Size: ',beta11),nrow(output)))
    Samplesize = c(Samplesize,rep(paste0('Sample Size: ',sampsize),nrow(output)))
  }
  plot.new()
  par(mai=c(0,0,0,0))
  plot.new()
  legend(x="center", ncol=2,legend=c('Theoretical PIP',paste0('Expected PIP')),
         col=c("black"),lty=c(2,3),lwd=2)

}


par(mfrow=c(1,1))

Gaussian_out = data.frame(cbind(PIP,pval,MSE_diff))
Gaussian_out$Ind_MSE = Gaussian_out$MSE_diff<0
Gaussian_out$MSE = Gaussian_out$Ind_MSE
Gaussian_out$MSE=replace(Gaussian_out$MSE,Gaussian_out$Ind_MSE==FALSE,"MSE1>MSE0")
Gaussian_out$MSE=replace(Gaussian_out$MSE,Gaussian_out$Ind_MSE==TRUE,"MSE1<MSE0")

Gaussian_out$Ind_PIP = Gaussian_out$PIP>0.5
Gaussian_out$PIP_Ind = Gaussian_out$Ind_PIP
Gaussian_out$PIP_Ind=replace(Gaussian_out$PIP_Ind,Gaussian_out$Ind_PIP==FALSE,"PIP<0.5")
Gaussian_out$PIP_Ind=replace(Gaussian_out$PIP_Ind,Gaussian_out$Ind_PIP==TRUE,"PIP>0.5")

Gaussian_out$pval_Ind = Gaussian_out$pval<0.05
Gaussian_out$pval_Ind=replace(Gaussian_out$pval_Ind,Gaussian_out$pval_Ind==FALSE,"pval>0.05")
Gaussian_out$pval_Ind=replace(Gaussian_out$pval_Ind,Gaussian_out$pval_Ind==TRUE,"pval<0.05")


Gaussian_out$EffectSize = EffectSize
Gaussian_out$SampleSize=Samplesize

#Gaussian_out$SampleSize= factor(Gaussian_out$SampleSize, levels=c('Sample Size: 20','Sample Size: 40','Sample Size: 60'))
Gaussian_out$SampleSize= factor(Gaussian_out$SampleSize, levels=c('Sample Size: 20','Sample Size: 40','Sample Size: 60','Sample Size: 100','Sample Size: 400'))

# head(Gaussian_out)
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) +  geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# 
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
#   stat_function(fun = f_pip, args = list(n = 100,p=pval), colour = "red")+
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 
# 
# 
# 
# 
# ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval))+
#  stat_function(fun = function(x) do.call(f_pip, c(n=400, list(x))), color = "blue") +
#   facet_wrap(~ EffectSize + SampleSize, scales = "free")
# 


f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-2))
  )
}

for(beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
    par(mfrow=c(1,1))
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "pip_cond1"
    colnames(output)[8] = "pip_cond2"
    pval = output[,"pval_mod1"]
    pip_CV5 =  output[,"pip_CV5"]
    pip_LOO =  output[,"pip_LOO"]
    plot(pval,pip_CV5)
    title(paste(beta11,sampsize,sep="   "),outer=TRUE,line=-2)
  }}

    mse_diff = output[,"mse1_CV"]-output[,"mse0_CV"]
    plot(pval,pip)
    lines(seq(0,max(max(pvals),4e-10),length=1000),sapply(seq(0,max(max(pvals),4e-10),length=1000),f_pip,n=sampsize),col="red",lwd=2)
    plot(mse_diff,pip)
    abline(h=0.5,col="red",lwd=2)
    abline(v=0,col="red",lwd=2)
    title(paste(beta11,sampsize,sep="   "),outer=TRUE,line=-2)
    par(mfrow=c(1,2))
    plot(output[,"pip_cond2"],output[,"pip_cond1"])
    plot(output[,"pip_cond2"],output[,"pip_CV5"])
    par(mfrow=c(1,2))
    plot(density(pip))
    plot(density(pval))
    title(paste(beta11,sampsize,sep="   "),outer=TRUE,line=-2)
  }}



set.seed(8)
sampsize=20
prop1=0.5
sigma=2
# Sample dataset with specified sample size
N = sampsize
X = c(rep(0,(1-prop1)*sampsize),rep(1,prop1*sampsize))
U <- rnorm(sampsize, sd = sigma)
Y = 10 + -1*X + U

training=data.frame("x"=X,"y"= Y)
K=nrow(training)
yourData<-training[sample(nrow(training)),]
folds <- cut(seq(1,nrow(yourData)),breaks=K,labels=FALSE)



pip_cv = c()
for(i in 1:nrow(yourData)){

testIndexes <- which(folds==i,arr.ind=TRUE)
testData <- yourData[testIndexes, ]
trainData <- yourData[-testIndexes, ]

mod1 = glm(y ~ x, family="gaussian", data=trainData)
mod0 = glm(y ~ 1, family="gaussian", data=trainData)

pred0 = predict(mod0,testData,type="response")
pred1 = predict(mod1,testData,type="response")
pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2) + 0.5*mean((pred1-testData$y)^2 == (pred0-testData$y)^2))

plot(trainData$x,trainData$y)
abline(h=coef(mod0),lwd=2,lty=2)
abline(a=coef(mod1)[1],b=coef(mod1)[2],lwd=2,lty=2,col="red")
points(testData$x,testData$y,col="darkgreen",pch=19)
title(paste(c((pred1-testData$y)^2,  (pred0-testData$y)^2)))
}

PIP_cv = mean(pip_cv)

mod1 = glm(y ~ x, family="gaussian", data=training)
mod0 = glm(y ~ 1, family="gaussian", data=training)
summary(mod1)

plot(training$x,training$y)
abline(h=coef(mod0),lwd=2,lty=1)
abline(a=coef(mod1)[1],b=coef(mod1)[2],lwd=2,lty=2,col="red")
}

output[8,c('pip_full','pip_CV5','pip_SS','pip_LOO')]

legend('bottom',c("m0","m1"),lty=c(1,2),col=c('black','red'),lwd=c(2,2))






# Presentation ISNPS2022

for(beta11 in c(0,-1,-4)){
  layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4,1))
  par(mai=rep(0.5, 4))
  for( sampsize in c(40,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "C1"
    colnames(output)[2] ="Exp"
    colnames(output)[10] = "C2"
    colnames(output)[5] ="SS"
    colnames(output)[7] ="SS_rep100"
    colnames(output)[8] ="LOO"
    colnames(output)[12] ="CV5"
    colnames(output)[15] ="rep_CV5"
    means <- apply(output,2,mean)
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    #boxplot(output[,c("C1","C2","Exp","SS","LOO","CV5","rep_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    boxplot(output[,c("C1","C2","Exp")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    abline(h=true_vals$PIP_theor,col="red",lwd=2)
    abline(h=true_vals$PIP_exp ,col='blue',lwd=2)
    #points(means[c("C1","C2","Exp","SS","LOO","CV5","rep_CV5")],pch=20)
    points(means[c("C1","C2","Exp")],pch=20)
    
  }
  plot.new()
  par(mai=c(0,0,0,0))
  legend(x="center", ncol=2,legend=c('Theoretical PIP',paste0('Expected PIP')),
         col=c("red","blue"),lty=c(1,1),lwd=2)

}


for(beta11 in c(0,-1,-4)){
  layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4,1))
  par(mai=rep(0.5, 4))
  for( sampsize in c(40,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "C1"
    colnames(output)[2] ="Exp"
    colnames(output)[10] = "C2"
    colnames(output)[5] ="SS"
    colnames(output)[7] ="SS_rep100"
    colnames(output)[8] ="LOO"
    colnames(output)[12] ="CV5"
    colnames(output)[15] ="rep_CV5"
    means <- apply(output,2,mean)
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    #boxplot(output[,c("C1","C2","Exp","SS","LOO","CV5","rep_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    boxplot(output[,c("C1","C2","Exp","emp_cond")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    abline(h=true_vals$PIP_theor,col="red",lwd=2)
    abline(h=true_vals$PIP_exp ,col='blue',lwd=2)
    #points(means[c("C1","C2","Exp","SS","LOO","CV5","rep_CV5")],pch=20)
    points(means[c("C1","C2","Exp","emp_cond")],pch=20)
    
  }
  plot.new()
  par(mai=c(0,0,0,0))
  legend(x="center", ncol=2,legend=c('Theoretical PIP',paste0('Expected PIP')),
         col=c("red","blue"),lty=c(1,1),lwd=2)
  
}







# Compare with empirical conditional PIP

for(beta11 in c(0,-1,-4)){
  layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(4, 4,1))
  par(mai=rep(0.5, 4))
  for( sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "C1"
    colnames(output)[2] ="Exp"
    colnames(output)[10] = "C2"
    colnames(output)[5] ="SS"
    colnames(output)[7] ="SS_rep100"
    colnames(output)[8] ="LOO"
    colnames(output)[12] ="CV5"
    colnames(output)[15] ="rep_CV5"
    
    colnames(output)
    use_diff  = output[,c("C1","C2","LOO","CV5","rep_CV5","SS","SS_rep100")] - output[,"emp_cond"]
    means = apply(use_diff,2,mean)
    boxplot(use_diff,main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    abline(h=0,lwd=2,lty=2)
    points(means,pch=20)
  }
  plot.new()
}





# Check consistency
par(mfrow=c(2,2))
for(beta11 in c(0,-1,-4)){
  for( sampsize in c(400,100,60,40,20)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "pip_C1"
    colnames(output)[10] = "pip_C2"
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    means <- apply(output,2,mean)
    print(means)
    if(sampsize==400){plot(density(output[,"pip_C1"]))}
    else{lines(density(output[,"pip_C1"]))}
  }
  abline(v=true_vals$PIP_theor,lwd=2,lty=2)
}






relation = c()
EffectSize = c()
Samplesize = c()
x = c()
for(beta11 in c(0,-1,-4)){
  par(mfrow=c(2,3))
  for(sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("/Volumes/GoogleDrive/My Drive/SummaryPIP/R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "pip_C1"
    colnames(output)[10] = "pip_C2"

    f_pip = function(n,p){
      return(
        pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-2))
      )
    }
    
    pvals = output[,"pval_mod1"]
    relation = c(relation,sapply(seq(0,max(max(pvals),0.05),length=1000),f_pip,n=sampsize))
    x = c(x,seq(0,max(max(pvals),0.05),length=1000))
    EffectSize = c(EffectSize,rep(paste0('Effect Size: ',beta11),1000))
    Samplesize = c(Samplesize,rep(paste0('Sample Size: ',sampsize),1000))
  }}

df_relation = data.frame(cbind("EffectSize"=EffectSize,"SampleSize"=Samplesize,"Relation"=relation,"x"=x))
df_relation$Relation = as.numeric(df_relation$Relation)
df_relation$x = as.numeric(df_relation$x)

df_relation$SampleSize= factor(df_relation$SampleSize, levels=c('Sample Size: 20','Sample Size: 40','Sample Size: 60','Sample Size: 100','Sample Size: 400'))
Gaussian_out$SampleSize= factor(Gaussian_out$SampleSize, levels=c('Sample Size: 20','Sample Size: 40','Sample Size: 60','Sample Size: 100','Sample Size: 400'))




ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval, color = MSE),alpha = 0.5) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  geom_line(aes(y = Relation,x=x), data = subset( df_relation,EffectSize=="Effect Size: 0"), colour = "black",linetype=1) +
  facet_wrap(~ EffectSize + SampleSize, scales = "free")+ scale_colour_grey()

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = pval, color = MSE),alpha = 0.5) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  geom_line(aes(y = Relation,x=x), data = subset( df_relation,EffectSize=="Effect Size: -1"), colour = "black",linetype=1) +
  facet_wrap(~ EffectSize + SampleSize, scales = "free")+ scale_colour_grey()

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = PIP, x = pval, color = MSE),alpha = 0.5) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  geom_line(aes(y = Relation,x=x), data = subset( df_relation,EffectSize=="Effect Size: -4"), colour = "black",linetype=1) +
  facet_wrap(~ EffectSize + SampleSize, scales = "free")+ scale_colour_grey()


# Summary Table

prop = c()
effect = c()
samplesize = c()
measure = c()
for(beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
colnames(output)
  if(beta11==0){
    correct_pval = 100*mean(output[,"pval_mod1"]>0.05)
    correct_LOO = 100*mean(output[,"pip_LOO"]<0.5)
    correct_MSE = 100*mean(output[,"mse1_rep_CV"]>output[,"mse0_rep_CV"])
    correct_repCV5 = 100*mean(output[,"pip_rep_CV5"]<0.5)
    correct_CV5 = 100*mean(output[,"pip_CV5"]<0.5)
  }
  else{
    correct_pval = 100*mean(output[,"pval_mod1"]<0.05)
    correct_LOO = 100*mean(output[,"pip_LOO"]>0.5)
    correct_MSE = 100*mean(output[,"mse1_rep_CV"]<output[,"mse0_rep_CV"])
    correct_repCV5 = 100*mean(output[,"pip_rep_CV5"]>0.5)
    correct_CV5 = 100*mean(output[,"pip_CV5"]>0.5)
  }
  prop = c(prop,c(correct_pval,correct_repCV5,correct_MSE))
  effect = c(effect,rep(beta11,3))
  samplesize = c(samplesize,rep(sampsize,3))
  measure = c(measure,c("p-value","PIP","MSE"))
  cat(beta11,"&",sampsize,"&",format(round(correct_pval, digits=2), nsmall = 2),"&",format(round(correct_MSE, digits=2), nsmall = 2),"&",format(round(correct_LOO, digits=2), nsmall = 2),"&",format(round(correct_CV5, digits=2), nsmall = 2),"&",format(round(correct_repCV5, digits=2), nsmall = 2),paste0("\\","\\"),"\n")
}
}
library(ggplot2)
out = data.frame(cbind(prop,effect,samplesize,measure))
out$effect = as.factor(out$effect)
out$prop = as.numeric(out$prop)
out$samplesize = as.numeric(out$samplesize)
ggplot(data=out)+geom_line(aes(x=samplesize,y=prop,linetype=effect,color= measure))






# Coverage

for (beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
    cl <- makeCluster(7)
    registerDoSNOW(cl)
    iterations <- 10000
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    output <- foreach(i=1:iterations,.packages=c("MASS","Matrix","mvnfast","caret"),
                      .options.snow = opts, .combine = rbind,.verbose = T) %dopar% { #, .errorhandling="remove"
                        result <- do_SIM_coverage(i,10,beta11,2,0.5,sampsize)
                        PIP = result$PIP
                        PIP_lower = result$PIP_lower
                        PIP_upper = result$PIP_upper
                        
                        return(cbind(PIP,PIP_lower,PIP_upper))
                      }
    close(pb)
    stopCluster(cl)
    
    save(output,file=paste0(paste(paste0("R/Output_Sims/sim_coverage",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
  }
  
}


library(glue)

for(beta11 in c(0,-1,-4)){
  for(sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_coverage",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    if(beta11 == 0){coverage = mean(output[,"PIP_lower"]<=0.5)}
    if(beta11 != 0) {coverage = mean(output[,"PIP_lower"]>0.5)}
    string="Effect size= {beta11}; Sample size= {sampsize}; Coverage= {coverage}."
    print(glue(string))
}}








