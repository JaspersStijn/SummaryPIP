pvals = c()
pip = c()
for(i in 1:1000){
sampsize = 4000
beta11 = 0
set.seed(i)
beta01=0

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



set.seed(1988)
trainIndex <- createDataPartition(dat$X, p = .5,
                                  list = FALSE,
                                  times = 50)

PIP = c()
for(j in 1:ncol(trainIndex)){
  training = dat[ trainIndex[,j],]
  testing  = dat[-trainIndex[,j],]

  sub=training
  mod1 = lm(Y ~ X, data = sub)
  mod0 = lm(Y ~ 1, data = sub)

  PIP = c(PIP,mean((testing$Y-predict(mod1,testing))^2<(testing$Y-predict(mod0,testing))^2))
}

pip = c(pip,mean(PIP))

#fhat1 = (testing$Y - predict(mod0,testing))^2 - (testing$Y - predict(mod1,testing))^2 + predict(mod1,testing)^2 - predict(mod0,testing)^2
fhat = (testing$Y)^2 - (testing$Y - predict(mod1,testing))^2 + predict(mod1,testing)^2

fbar = mean(fhat)
fvar = mean((fhat-fbar)^2)

fstat = sqrt(nrow(testing)) * fbar / sqrt(fvar)
pvals = c(pvals,1-pnorm(abs(fstat)))
}


plot(pvals,pip)
abline(h=0.5,lwd=2,lty=2,col="red")
abline(v=0.05,lwd=2,lty=2,col="red")

mean(pvals>0.05)


