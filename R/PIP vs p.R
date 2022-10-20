# PIP vs. pvalue vs n

PIP = seq(0.5,1,length=1000)
pvals = seq(0,1,length=1000)
sample_sizes = seq(10,300,length=1000)

f_p = function(n,pip){
  return(2*(1-pnorm(2*sqrt(n)*qnorm(pip))))
}

f_pip = function(n,p){
  return(
    pnorm(1/(2*sqrt(n))*qnorm(1-0.5*p))
  )
}

par(mfrow=c(2,2))

par(mfrow=c(1,1))
plot(sample_sizes,f_p(sample_sizes,0.5),lwd=2,lty=1,type='l',col=1,ylim=c(0,1),ylab="p-value",xlab="Sample size",xlim=c(0,450),main="p-value vs n")
lines(sample_sizes,f_p(sample_sizes,0.52),lwd=2,col=2)
lines(sample_sizes,f_p(sample_sizes,0.54),lwd=2,col=3)
lines(sample_sizes,f_p(sample_sizes,0.56),lwd=2,col=4)
lines(sample_sizes,f_p(sample_sizes,0.58),lwd=2,col=5)
lines(sample_sizes,f_p(sample_sizes,1),lwd=2,col=6)
abline(h=0.05,lwd=2,lty=2)
legend("topright",c(paste0("PIP= ",c(seq(0.5,0.58,0.02),1)),"5% cut-off"),lwd=2,lty=c(rep(1,6),2),col=c(1:6,1))

plot(sample_sizes,f_pip(sample_sizes,0.00001),lwd=2,lty=1,type='l',col=1,ylim=c(0.5,0.7),ylab="PIP",xlab="Sample size",xlim=c(0,450),main="PIP vs n")
lines(sample_sizes,f_pip(sample_sizes,0.05),lwd=2,col=2)
lines(sample_sizes,f_pip(sample_sizes,0.1),lwd=2,col=3)
lines(sample_sizes,f_pip(sample_sizes,0.5),lwd=2,col=4)
lines(sample_sizes,f_pip(sample_sizes,0.7),lwd=2,col=5)
lines(sample_sizes,f_pip(sample_sizes,1),lwd=2,col=6)
legend("topright",c(paste0("p= ",c(0.00001,0.05,0.1,0.5,0.7,1))),lwd=2,lty=c(rep(1,6)),col=c(1:6))

plot(pvals,f_pip(10,pvals),lwd=2,lty=1,type='l',col=1,ylim=c(0.5,0.7),ylab="PIP",xlab="p-value",xlim=c(0,1),main="PIP vs p-value")
lines(pvals,f_pip(20,pvals),lwd=2,col=2)
lines(pvals,f_pip(30,pvals),lwd=2,col=3)
lines(pvals,f_pip(40,pvals),lwd=2,col=4)
lines(pvals,f_pip(50,pvals),lwd=2,col=5)
lines(pvals,f_pip(60,pvals),lwd=2,col=6)
abline(v=0.05,lwd=2,lty=2)
legend("topright",c(paste0("n= ",seq(10,60,10)),"5% cut-off"),lwd=2,lty=c(rep(1,6),2),col=c(1:6,1))






# 
# plot(log(pvals),-10*log(1+4*qnorm(f_pip(10,pvals))^2),lwd=2,lty=1,type='l',col=1,ylab="PIP",xlab="P-values",main="PIP vs p")
# lines(log(pvals),-20*log(1+4*qnorm(f_pip(20,pvals))^2),lwd=2,col=2)
# lines(log(pvals),-30*log(1+4*qnorm(f_pip(30,pvals))^2),lwd=2,col=3)
# lines(log(pvals),-40*log(1+4*qnorm(f_pip(40,pvals))^2),lwd=2,col=4)
# lines(log(pvals),-50*log(1+4*qnorm(f_pip(50,pvals))^2),lwd=2,col=5)
# lines(log(pvals),-60*log(1+4*qnorm(f_pip(60,pvals))^2),lwd=2,col=6)
# abline(a=0,b=1)
# #abline(v=0.05,lwd=2,lty=2)
# legend("topright",c(paste0("n= ",seq(10,60,10)),"5% cut-off"),lwd=2,lty=c(rep(1,6),2),col=c(1:6,1))
# 
# par(mfrow=c(1,1))
# 
# exp(-1000*log(1+1/16))

