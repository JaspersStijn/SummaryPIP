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
                        return(cbind(pip_cond,pip_exp,emp_cond,emp_exp,
                                     pip_SS,pip_SS_rep50,pip_SS_rep100,pip_LOO,pval_mod1,pip_full,pip_cond_check,pip_CV5,mse0_CV,mse1_CV))
                      }
    close(pb)
    stopCluster(cl)

    save(output,file=paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
  }

}




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
    colnames(output)[1] = "pip_C1"
    colnames(output)[10] = "pip_C2"
    means <- apply(output,2,mean)
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    #boxplot(output[,c("pip_cond1","pip_cond2","pip_exp")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_cond1","pip_cond2","pip_exp","pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_C1","pip_C2","emp_cond","pip_exp","emp_exp" ,"pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    #boxplot(output[,c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_SS","pip_CV5")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    boxplot(output[,c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_CV5","pip_SS","pip_SS_rep50","pip_SS_rep100")],main=paste(paste0("Sample size: ",sampsize),paste0("Beta1: ",beta11),sep= "\n"))
    abline(h=true_vals$PIP_theor,col="red",lwd=2)
    abline(h=true_vals$PIP_exp ,col='blue',lwd=2)
    abline(h=true_vals$PIP_exp_limit ,col='darkgreen',lty=2,lwd=2)
    #points(means[c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_SS","pip_CV5")],pch=25)
    points(means[c("pip_C1","pip_C2","pip_exp","pip_LOO","pip_CV5","pip_SS","pip_SS_rep50","pip_SS_rep100")],pch=20)

    print(means)
    PIP = c(PIP,output[,c("pip_CV5")])
    pval = c(pval,output[,c("pval_mod1")])
    MSE_diff = c(MSE_diff,output[,c("mse1_CV")]-output[,c("mse0_CV")])
    EffectSize = c(EffectSize,rep(paste0('Effect Size: ',beta11),nrow(output)))
    Samplesize = c(Samplesize,rep(paste0('Sample Size: ',sampsize),nrow(output)))
  }
  plot.new()
  par(mai=c(0,0,0,0))
  plot.new()
  legend(x="center", ncol=3,legend=c('Theoretical PIP',paste0('Expected PIP for specified sampsize'),'Expected PIP for n=100000000' ),
         col=c("red","blue","darkgreen"),lty=c(1,1,2),lwd=2)

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

head(Gaussian_out)
ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = PIP, x = pval, color = MSE)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")


ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) +  geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")

ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = PIP, x = MSE_diff, color = pval_Ind)) + geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")




ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = pval, x =MSE_diff, color = PIP_Ind)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")
ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -2"), aes(y = pval, x = MSE_diff, color = PIP_Ind)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")
ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -4"), aes(y = pval, x = MSE_diff, color = PIP_Ind)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")




ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: 0"), aes(y = PIP, x = pval)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")
ggplot() + geom_point(data = subset(Gaussian_out,EffectSize=="Effect Size: -1"), aes(y = PIP, x = pval, color = pval_Ind)) + #geom_hline(yintercept = 0.5,linetype="dashed")+ geom_vline(xintercept = 0.05,linetype="dashed")+
  facet_wrap(~ EffectSize + SampleSize, scales = "free")





f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qt(1-0.5*p,df=n-1))
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
    colnames(output)[1] = "pip_C1"
    colnames(output)[10] = "pip_C2"
    means <- apply(output,2,mean)
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    boxplot(output[,c("pip_C1","pip_C2","pip_exp")])
    abline(h=true_vals$PIP_theor,col="red",lwd=2)
    abline(h=true_vals$PIP_exp ,col='blue',lwd=2)
    abline(h=true_vals$PIP_exp_limit ,col='darkgreen',lty=2,lwd=2)
    points(means[c("pip_C1","pip_C2","pip_exp")],pch=20)
  }
  plot.new()
  par(mai=c(0,0,0,0))
  legend(x="center", ncol=3,legend=c('Theoretical PIP',paste0('Expected PIP for specified sampsize'),'Expected PIP for n=100000000' ),
         col=c("red","blue","darkgreen"),lty=c(1,1,2),lwd=2)

}








# Compare with empirical conditional/expected PIP
for(beta11 in c(0,-1,-4)){
  layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(4, 4,1))
  par(mai=rep(0.5, 4))
  for( sampsize in c(20,40,60,100,400)){
    use = load(paste0(paste(paste0("R/Output_Sims/sim_effect_new_exp",abs(beta11)),paste0("sampsize",sampsize),sep="_"),".R"))
    output = get(use)
    colnames(output)[1] = "pip_C1"
    colnames(output)[10] = "pip_C2"
    true_vals=TRUE_PIPs_faster(10,beta11,2,0.5,sampsize)
    means <- apply(output,2,mean)
    print(means)
    }
  plot.new()
  par(mai=c(0,0,0,0))
  plot.new()
  legend(x="center", ncol=3,legend=c('Theoretical PIP',paste0('Expected PIP for specified sampsize'),'Expected PIP for n=100000000' ),
         col=c("red","blue","darkgreen"),lty=c(1,1,2),lwd=2)

}






