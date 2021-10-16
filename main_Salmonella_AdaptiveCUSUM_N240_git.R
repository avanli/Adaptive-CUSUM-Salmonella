# Author: O. Arda Vanli
# Institution: FAMU-FSU College of Engineering
# Date: 10-14-2021
# baseline data is modeled with seasonal INGARCH(1,1) with linear secular trend intercept.
# surveillance method is tuned with seasonal INGARCH(1,1) with constant intercept.
# (1) to plot data and INGARCH models, go to chunk "plot data and model"
# (2) to directly generate CUSUM plots go to chunk "monitor observations from N+1 onwards"
# (3) to tune charts, go to chunk " find alarm thresholds"
library('surveillance')

data(shadar)
#Number of weekly salmonella hadar cases in Germany 2001-2006
n=length(shadar$observed)
flag_ChartTuning=FALSE        # set TRUE to run the chunk for chart tuning 
# used for chart ploting 
path="C:\\Users\\ardav\\Dropbox\\Publications\\MonitorSeasonalPoisson\\R\\case Study Salmonella\\"
# version 2 assumes intercept with secular trend - negative trend causes issues with long ARL0
source(paste(path,"get_ARL_CUSUM2.R",sep=""))                
source(paste(path,"get_ARL_CUSUMAdaptive2.R",sep=""))
source(paste(path,"get_ARL_CUSUMAdaptiveDES2.R",sep=""))
source(paste(path,"fitINARMA_Seasonal2.R",sep=""))
source(paste(path,"fitINARMA_Seasonal.R",sep=""))
source(paste(path,"plotAdaptiveCUSUM.R",sep=""))
source(paste(path,"plotCUSUM.R",sep=""))
#### plot data and model - fit a seasonal INGARCH(1,1) on onservations 1 to N ####
# period is 52 weeks
P=52
N=240
y=shadar$observed[1:N]
out_garch11=fitINARMA_Seasonal2(y,P,4)

# predictions
b0=as.numeric(out_garch11$model$coefficients[1])
b_cos=as.numeric(out_garch11$model$coefficients[2])
b_sin=as.numeric(out_garch11$model$coefficients[3])
phi1=as.numeric(out_garch11$phi)
theta1=as.numeric(out_garch11$theta)
dt=as.numeric(out_garch11$dt)

mdl.prmEst=c(as.numeric(out_garch11$model$coefficients[1]), phi1-theta1,theta1, 
             as.numeric(out_garch11$model$coefficients[2]), 
             as.numeric(out_garch11$model$coefficients[3]),P,as.numeric(out_garch11$dt))
# model summary
summary(out_garch11$model)
t=seq(1,N,1)
# harmonic mean
b1=phi1-theta1
b2=theta1
f_harmonic=b0/(1-b1-b2)+b_cos*cos(2*pi*t/P)+b_sin*sin(2*pi*t/P)+dt*t
eps_hat=0
for (i in 2:N){
  eps_hat[i]=y[i]-phi1*y[i-1]+theta1*eps_hat[i-1]-
    b0-b_cos*cos(2*pi*(i-1)/P)-b_sin*sin(2*pi*(i-1)/P)-dt*(i-1)
}
eps_hat_l=eps_hat[1:N-1] 
yl=y[1:N-1]
# one-step ahead forecasts for time periods 2,3,...,N
y_hat=b0+b_cos*cos(2*pi*t[2:N]/P)+b_sin*sin(2*pi*t[2:N]/P)+phi1*yl+theta1*eps_hat_l+dt*t[2:N]
windows(7,5)
par(mfrow=c(1,1),mar=c(3,3,1.5,0)+0.1,mgp=c(1.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(shadar$observed,pch=NA,xlab="week",ylab="Case Counts",xlim=c(1,295),ylim=c(-1,25))
lines(shadar$observed,type='h')
lines(t[2:N],y_hat,col='blue')
lines(t,f_harmonic,col='red',lwd=1)
legend("topright",legend=c('data','secular and seasonal mean','forecast'),
       col=c("black","red","blue"),lty=c(1,1,1),lwd=c(1,1,1))
windows(7,5)
resForecast=y[2:N]-y_hat
hist(resForecast,xlab = 'forecast error')
windows(7,5)
resHarmonic=y-f_harmonic
hist(resHarmonic,xlab = 'harmonic error')
# fit the no trend model
out_garch11_noTrend=fitINARMA_Seasonal(y,P,4)
summary(out_garch11_noTrend$model)
phi11=as.numeric(out_garch11_noTrend$phi)
theta11=as.numeric(out_garch11_noTrend$theta)
mdl.prmEst_NoTrend=c(as.numeric(out_garch11_noTrend$model$coefficients[1]), phi11-theta11,theta11, 
                     as.numeric(out_garch11_noTrend$model$coefficients[2]), 
                     as.numeric(out_garch11_noTrend$model$coefficients[3]),P)
#### monitor observations from N+1 onwards - compile the functions plotAdaptiveCUSUM and plotCUSUM ####
# hard coded model estimated from N =240 data. model with time trend: 
# lambda0[i]=b0+b1*y[i-1]+b2*lambda0[i-1]+b_cos*cos(2*pi/P*i)+b_sin*sin(2*pi/P*i)+dt*i
mdl.prmEst=c(4.03926280,  0.21294125, -0.05447749, -1.21130202, -0.47453502, 52.00000000, -0.01048445)
b0=mdl.prmEst[1];b1=mdl.prmEst[2];b2=mdl.prmEst[3];b_cos=mdl.prmEst[4];b_sin=mdl.prmEst[5];dt=mdl.prmEst[7]
# hard coded model estimated from N =240 data. model without time trend: 
# lambda0[i]=b00+b11*y[i-1]+b22*lambda0[i-1]+b_cos2*cos(2*pi/P*i)+b_sin2*sin(2*pi/P*i)
mdl.prmEst_NoTrend=c(0.96763979,  0.29234841,  0.41568894, -0.45458747, -0.06964397, 52.00000000)
b00=mdl.prmEst_NoTrend[1];b11=mdl.prmEst_NoTrend[2];b22=mdl.prmEst_NoTrend[3];
b_cos2=mdl.prmEst_NoTrend[4];b_sin2=mdl.prmEst_NoTrend[5];

# plot from plotTimeBegin, time_inj is when monitoring starts
time_inj=240 # monitor starts from 240 
plotTimeBegin=150
# monitoring with DEWMA
thresh.resid=4;flag_DEWMA=1
sm_prm=0.1;sm_prm2=0.05; h_sD=5.4 
resDEWMA=plotAdaptiveCUSUM(time_inj,b00,b11,b22,b_cos2,b_sin2,P,dt,shadar$observed,h_sD,flag_DEWMA,sm_prm,sm_prm2,thresh.resid)
# monitoring with EWMA
thresh.resid=Inf;flag_DEWMA=0
sm_prm=0.1; h_sE=0.93
resEWMA=plotAdaptiveCUSUM(time_inj,b00,b11,b22,b_cos2,b_sin2,P,dt,shadar$observed,h_sE,flag_DEWMA,sm_prm,sm_prm2,thresh.resid)
# monitoring with step cusum
delta_c=b00+1  ;k_inarch=delta_c/b00-1   ; beta_c=((1+k_inarch)*b00-b00);h_sC1=4.11
resCUSUMSmall=plotCUSUM(time_inj,b00,b11,b22,b_cos2,b_sin2,P,dt,beta_c,shadar$observed,h_sC1)
delta_c=b00+4  ;k_inarch=delta_c/b00-1   ; beta_c=((1+k_inarch)*b00-b00);h_sC2=4.07
resCUSUMLarge=plotCUSUM(time_inj,b00,b11,b22,b_cos2,b_sin2,P,dt,beta_c,shadar$observed,h_sC2)
delta_c=b00+0.5  ;k_inarch=delta_c/b00-1   ; beta_c=((1+k_inarch)*b00-b00);h_sC3=3.415
resCUSUMSmall2=plotCUSUM(time_inj,b00,b11,b22,b_cos2,b_sin2,P,dt,beta_c,shadar$observed,h_sC3)
windows(7,10)
shiftScale=2
par(mfrow=c(5,1),mar=c(3,3,1.5,0.5)+0.1,mgp=c(1.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
# panels 1 and 2 - DEWMA
plot(shadar$observed,pch=NA,xlab="week",ylab="Cases and DEWMA",xlim=c(plotTimeBegin,295),ylim=c(-4,22))
lines(shadar$observed,type='h')
points(time_inj, -1.5, pch=2,col='blue', xpd = TRUE,cex=1.2)
lines(seq(1,resDEWMA$RL,1),resDEWMA$yhat,col='blue',lwd=1)
lines(seq(1,resDEWMA$RL,1),resDEWMA$shift_hat*shiftScale,col='orange',lwd=1)
points(resDEWMA$RL, -1.5, pch=2,col='green', xpd = TRUE,cex=1.2,lwd=1)
plot(seq(1,resDEWMA$RL,1), resDEWMA$S[1:resDEWMA$RL]/h_sD,
     main="",xlab='time',ylab='DESCUSUM',xlim=c(plotTimeBegin,295),type='l')
points(time_inj, -0.2, pch=2,col='blue', xpd = TRUE,cex=1.2)
points(resDEWMA$RL, -0.2, pch=2,col='green', xpd = TRUE,cex=1.2,lwd=1)
text(resDEWMA$RL, -0.5,paste("alarm @ ",resDEWMA$RL,sep=""), xpd = TRUE)
abline(h=1)
# panels 3 and 4 - EWMA
shiftScale=1
plot(shadar$observed,pch=NA,xlab="week",ylab="Cases and EWMA",xlim=c(plotTimeBegin,295),ylim=c(-4,22))
lines(shadar$observed,type='h')
points(time_inj, -1.5, pch=2,col='blue', xpd = TRUE,cex=1.2)
lines(seq(1,resEWMA$RL,1),resEWMA$yhat,col='blue',lwd=1)
lines(seq(1,resEWMA$RL,1),resEWMA$shift_hat*shiftScale,col='orange',lwd=1)
points(resEWMA$RL, -1.5, pch=2,col='green', xpd = TRUE,cex=1.2,lwd=1)
plot(seq(1,resEWMA$RL,1), resEWMA$S[1:resEWMA$RL]/h_sE,
     main="",xlab='time',ylab='SESCUSUM',xlim=c(plotTimeBegin,295),type='l')
points(time_inj, -0.2, pch=2,col='blue', xpd = TRUE,cex=1.2)
points(resEWMA$RL, -0.2, pch=2,col='green', xpd = TRUE,cex=1.2,lwd=1)
text((resEWMA$RL-.75), -0.8,paste("alarm @ ",resEWMA$RL,sep=""), xpd = TRUE)
abline(h=1)
# panel 5 - step CUSUM
plot(seq(1,resCUSUMSmall$RL,1), resCUSUMSmall$S[1:resCUSUMSmall$RL]/h_sC1,type='l',
     xlab='time',ylab='Step shift CUSUM',xlim=c(plotTimeBegin,295))
abline(h=1)
lines(seq(1,resCUSUMLarge$RL,1), resCUSUMLarge$S[1:resCUSUMLarge$RL]/h_sC2,lty=2)
lines(seq(1,resCUSUMSmall2$RL,1), resCUSUMSmall2$S[1:resCUSUMSmall2$RL]/h_sC3,lty=3)
points(time_inj, -.3, pch=2,col='blue', xpd = TRUE,cex=1.2)
points(resCUSUMSmall2$RL, -.3, pch=2,col='green', xpd = TRUE,cex=1.2)
points(resCUSUMSmall$RL, -.3, pch=2,col='green', xpd = TRUE,cex=1.2)
points(resCUSUMLarge$RL, -.3, pch=2,lwd=1,col='green', xpd = TRUE,cex=1.2)
text((resCUSUMSmall$RL-.75), -0.5,paste("alarms @ ",resCUSUMLarge$RL,"&",resCUSUMSmall$RL,"&",resCUSUMSmall2$RL,sep=""), xpd = TRUE)

legend(150,0.7,legend=c('kappa*=0.5','kappa*=1.0','kappa*=4.0'),
       col=c("black","black","black"),
       lty=c(3,1,2),lwd=c(1,1,1),cex=0.75) #ncol=2)
#### find alarm thresholds using parameters estimated from N data: if flag_ChartTuning=TRUE####
if(flag_ChartTuning){
n_MC=1000
time_inj=1  # zero state ARL
#mdl.prm=c(b0,b1,b2,b_cos,b_sin,P,dt)
flag_trendShift=FALSE 
beta=0  # in control ARL
#EWMA
sm_prm=0.1  #EWMA smoothing parameter
thresh.resid=Inf
# h_s=0.93  gives ARL0=408.929 with  sm_prm=0.1,thresh.resid=Inf
h_s=.93
set.seed(98765)
res=get_ARL_CUSUMAdaptive2(beta,mdl.prmEst_NoTrend, h_s, time_inj, 
                          n_MC,flag_trendShift,sm_prm,thresh.resid)
paste("h_s=",h_s, ", sm_prm=", sm_prm,", thresh.resid=",thresh.resid, ", ARL=",res$ARL, ", seRL=",res$seRL,sep="")
# DEWMA
sm_prm=0.1  #EWMA smoothing parameter
sm_prm2=0.05
thresh.resid=4
# h_s=5.4 gives ARL0=407.561 with  sm_prm=0.1,sm_prm2=0.05,thresh.resid=4
h_s=5.4
set.seed(98765)
res=get_ARL_CUSUMAdaptiveDES2(beta,mdl.prmEst_NoTrend, h_s, time_inj, 
                             n_MC,flag_trendShift,sm_prm,thresh.resid,sm_prm2)
paste("h_s=",h_s, ", sm_prm=", sm_prm,", sm_prm2=", sm_prm2,", thresh.resid=",thresh.resid, ", ARL=",res$ARL, ", seRL=",res$seRL,sep="")
# step shift CUSUM
flag_trendCUSUM=FALSE
b00=mdl.prmEst_NoTrend[1]
delta_c=b00+.5  #b00+c means it is sensitive to c cases increase
k_inarch=delta_c/b00-1     
beta_c=((1+k_inarch)*b00-b00)
# h_s=3.415 gives ARL0=399.13400      with delta_c=b00+.5
# h_s=4.11 gives ARL0= 401.18200    with delta_c=b00+1
# h_s=4.07 gives ARL0= 399.62800    with delta_c=b00+4
# h_s=0.27 gives ARL0= 399.62800     with delta_c=b00+10
h_s=3.415
# for reproducibility
set.seed(98765)
res=get_ARL_CUSUM2(beta,beta_c,mdl.prmEst_NoTrend, h_s, time_inj, 
                  n_MC,flag_trendCUSUM,flag_trendShift)
print(c(res$ARL, res$seRL,h_s))
}
