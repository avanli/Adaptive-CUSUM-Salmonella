# Adaptive CUSUM for step shift detection
# version 2: it incorporates linear trend in intercept in addition to harmonic components. Used in Salmonella case study.
get_ARL_CUSUMAdaptive2<-function(beta,mdl.prm,h_s, time_inj, n_MC,flag_trendShift,sm_prm,thres.resid){
# if flag_trendShift theta=beta else beta is the ooc intercept value for step shift
# sm_prm: EWMA forecasting smoothing parameter
  
delta0=mdl.prm[1]
b1=mdl.prm[2]
b2=mdl.prm[3]
eta=mdl.prm[4]
psi=mdl.prm[5]
T=mdl.prm[6]
dt=mdl.prm[7]
RL=NA
tau_est=NA

j=1
while (j <=n_MC){
  S=0
  i=1
  lambda0=delta0/(1-b1-b2)
  y=lambda0
  lambda1=lambda0
  lambda1_c=lambda0
  tau=0
  yhat=lambda0  # initialize EWMA forecast
  muhat=0  # initialize EWMA smoother
  delta_hat=delta0
  shift_hat=0
  muhatVec=0
  yhatINGARCH=lambda0
  while (S[i]<h_s| i<=time_inj){
    i=i+1
    if (i>time_inj){
      if(flag_trendShift){
        theta=beta
        delta1=delta0+theta*(i-time_inj)  # trend shift
      }else{
        delta1=delta0+beta  # step shift
      }
    }else{
      delta1=delta0
      #lambda1[i]=lambda0[i]
    }
    # baseline rate
    lambda0[i]=delta0+b1*y[i-1]+b2*lambda0[i-1]+eta*cos(2*pi/T*i)+psi*sin(2*pi/T*i)+dt*i
    if(lambda0[i]<0){lambda0[i]=0}
    # actual rate
    lambda1[i]=delta1+b1*y[i-1]+b2*lambda1[i-1]+eta*cos(2*pi/T*i)+psi*sin(2*pi/T*i)+dt*i
    if (lambda1[i]>0){
      y[i]=rpois(1,lambda=lambda1[i])  # Poisson Data
    }else{
      y[i]=0
    }
    # INGARCH forecast of level
    yhatINGARCH[i]=lambda0[i]
    resid=y[i]-yhatINGARCH[i]-muhat*(1-b1-b2)    #V1
    #resid=y[i]-yhatINGARCH[i]-muhat             #V2
    
    # Huber's score function
    if(resid< (-thres.resid)){modified.resid=resid+(1-sm_prm)*thres.resid}
    else if(abs(resid)<thres.resid){modified.resid=sm_prm*resid}
    else{modified.resid=resid-(1-sm_prm)*thres.resid}
    muhat=muhat+modified.resid

    yhat[i]=muhat*(1-b1-b2)+yhatINGARCH[i]      #V1
    #yhat[i]=muhat+yhatINGARCH[i]               #V2
    
    shift_hat[i]=muhat*(1-b1-b2)
    
    muhatVec[i]=muhat
    k=i-1
    if (lambda1[i]>0 & lambda0[i]>0 & yhat[i]>0){
      weight=shift_hat[i] 
      increment=weight*(y[i]*log(yhat[i]/lambda0[i])-(yhat[i]-lambda0[i]))
      S[i]=max(0,S[i-1]+increment)
    }else{
      S[i]=S[i-1]
    }
    # reset tau (change point estimator for each subregion) to k if cusum = 0
    if (S[i]==0){tau=k}
    #print(i)
  }
  if (i>time_inj){
    #print(paste("MC run ",j,sep=""))
    RL[j]=i-time_inj
    tau_est[j]=tau
    j=j+1
  }
}

ARL=mean(RL)
seRL=sqrt(var(RL)/n_MC)

# bias and mse of change point estimates
bias=sum(tau_est-time_inj)/n_MC
mse=1/n_MC*(t(tau_est-time_inj)%*%(tau_est-time_inj))


# return results
result=list(ARL=ARL,seRL=seRL,tau_est=tau_est,bias=bias,
            mse=mse,SLastRun=S,YLastRun=y,lam0Out=lambda0,yhat=yhat,delta_hat=delta_hat,shift_hat=shift_hat,muhatVec=muhatVec)
return(result)

}