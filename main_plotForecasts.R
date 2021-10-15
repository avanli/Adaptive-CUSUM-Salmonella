library('surveillance')

data(shadar)
#Number of weekly salmonella hadar cases in Germany 2001-2006
n=length(shadar$observed)
path="C:\\Users\\ardav\\Dropbox\\Publications\\MonitorSeasonalPoisson\\R\\case Study Salmonella\\"
source(paste(path,"fitINARMA_Seasonal.R",sep=""))


#### plot data and model - fit a seasonal INGARCH(1,1) on onservations 1 to 200 ####
# period is 52 weeks
P=52


getForecasts=function(){
  out_garch11=fitINARMA_Seasonal(y,P,4)
  
  # predictions
  b0=as.numeric(out_garch11$model$coefficients[1])
  b_cos=as.numeric(out_garch11$model$coefficients[2])
  b_sin=as.numeric(out_garch11$model$coefficients[3])
  phi1=as.numeric(out_garch11$phi)
  theta1=as.numeric(out_garch11$theta)
  dt=as.numeric(out_garch11$dt)
  # model summary
  t=seq(1,N,1)
  # harmonic mean
  b1=phi1-theta1
  b2=theta1
  f_harmonic=b0/(1-b1-b2)+b_cos*cos(2*pi*t/P)+b_sin*sin(2*pi*t/P)+dt*t
  #f_harmonic=b0+b_cos*cos(2*pi*t/P)+b_sin*sin(2*pi*t/P)+dt*t
  eps_hat=0
  for (i in 2:N){
    eps_hat[i]=y[i]-phi1*y[i-1]+theta1*eps_hat[i-1]-
      b0-b_cos*cos(2*pi*(i-1)/P)-b_sin*sin(2*pi*(i-1)/P)-dt*(i-1)
  }
  eps_hat_l=eps_hat[1:N-1] 
  yl=y[1:N-1]
  # one-step ahead forecasts for time periods 2,3,...,N
  y_hat=b0+b_cos*cos(2*pi*t[2:N]/P)+b_sin*sin(2*pi*t[2:N]/P)+phi1*yl+theta1*eps_hat_l+dt*t[2:N]

  return(list(y_hat=y_hat))
  }
N1=240
y=shadar$observed[1:N1]
res1=getForecasts()
err1=y[2:N1]-res1$y_hat[2:N1]
N2=200
y=shadar$observed[1:N2]
res2=getForecasts()
err2=y[2:N2]-res2$y_hat[2:N2]
N3=170
y=shadar$observed[1:N3]
res3=getForecasts()
err3=y[2:N3]-res3$y_hat[2:N3]
png("Plot3.png", width = 7, height = 5, units = 'in', res = 600)
#windows(7,5)
par(mfrow=c(1,1),mar=c(3,3,1.5,0)+0.1,mgp=c(1.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(shadar$observed,pch=NA,xlab="week",ylab="Case Counts",xlim=c(1,295),ylim=c(-1,25))
lines(shadar$observed,type='h',col='grey')
lines(t[2:N1],res1$y_hat[2:N1],col='darkorange3',lwd=2)
lines(t[2:N2],res2$y_hat[2:N2],col='tan2',lwd=2,lty=1)
#lines(t[2:N3],res3$y_hat[2:N3],col='tan2')
dev.off()

windows(7,5)
par(mfrow=c(1,3),mar=c(3,3,1.5,0)+0.1,mgp=c(1.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
hist(err1,xlim=c(-20,30),breaks=seq(-10,25,5))
hist(err2,xlim=c(-20,30),breaks=seq(-10,25,5))
hist(err3,xlim=c(-20,30),breaks=seq(-10,25,5))

