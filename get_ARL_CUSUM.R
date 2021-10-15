# average  run length designed to detect trends or step shifts under trend shifts
get_ARL_CUSUM<-function(beta,beta_c, mdl.prm,h_s, time_inj, n_MC,flag_trendCUSUM,flag_trendShift){

# if flag_trendCUSUM theta_c=beta_c else beta_c is tuned ooc intercept value for step CUSUM
# if flag_trendShift theta=beta else beta is the ooc intercept value for step shift
  
delta0=mdl.prm[1]
b1=mdl.prm[2]
b2=mdl.prm[3]
eta=mdl.prm[4]
psi=mdl.prm[5]
T=mdl.prm[6]
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
    }
    # baseline rate
    lambda0[i]=delta0+b1*y[i-1]+b2*lambda0[i-1]+eta*cos(2*pi/T*i)+psi*sin(2*pi/T*i)
    # actual rate
    lambda1[i]=delta1+b1*y[i-1]+b2*lambda1[i-1]+eta*cos(2*pi/T*i)+psi*sin(2*pi/T*i)
    if (lambda1[i]>0){
      y[i]=rpois(1,lambda=lambda1[i])  
    }else{
      y[i]=0
    }
    k=i-1
    if(flag_trendCUSUM){
      # trend shift cusum
      theta_c=beta_c
      delta1_c=delta0+theta_c*(k-tau+1)            # use change point estimate
      #delta1_c=delta0+theta_c*(k-time_inj+1)        # use true change point 
    }else{
      # sustained shift cusum
      delta1_c=delta0+beta_c
    }
    lambda1_c[i]=delta1_c+b1*y[i-1]+b2*lambda1_c[i-1]+eta*cos(2*pi/T*i)+psi*sin(2*pi/T*i)
    if (lambda1[i]>0 & lambda0[i]>0 & lambda1_c[i]>0){
      S[i]=max(0,S[i-1]+y[i]*log(lambda1_c[i]/lambda0[i])-(lambda1_c[i]-lambda0[i]))
    }else{
      S[i]=S[i-1]
    }
    # reset tau (change point estimator for each subregion) to k if cusum = 0
    if (S[i]==0){tau=k}
    #print(c(i,S[i]))
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
result=list(ARL=ARL,seRL=seRL,tau_est=tau_est,bias=bias,mse=mse,SLastRun=S,YLastRun=y,lamOut=lambda1)
return(result)

}