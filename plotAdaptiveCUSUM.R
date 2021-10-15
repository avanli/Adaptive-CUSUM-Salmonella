plotAdaptiveCUSUM<-function(time_inj,b0,b1,b2,b_cos,b_sin,P,dt,y,h_s,flag_DEWMA,sm_prm,sm_prm2,thresh.resid){
  
  
  
  lambda0=b0/(1-b1-b2)
  i=1
  
  yhatINGARCH=lambda0
  yhat=lambda0
  shift_hat=0
  S=0
  muhat=0
  SlopePrev=0
  LevelPrev=lambda0
  
  while (S[i]<h_s|i<=time_inj){
    i=i+1
    k=i-1
    # baseline rate
    #lambda0[i]=b0+b1*y[i-1]+b2*lambda0[i-1]+b_cos*cos(2*pi/P*i)+b_sin*sin(2*pi/P*i)+dt*i
    lambda0[i]=b0+b1*y[i-1]+b2*lambda0[i-1]+b_cos*cos(2*pi/P*i)+b_sin*sin(2*pi/P*i)
    # INGARCH forecast of leve
    yhatINGARCH[i]=lambda0[i]
    resid=y[i]-yhatINGARCH[i]-muhat*(1-b1-b2)   
    
    # use Huber's score function
    if(resid< (-thresh.resid)){modified.resid=resid+(1-sm_prm)*thresh.resid}
    if(abs(resid)<thresh.resid){modified.resid=sm_prm*resid}
    if(resid> thresh.resid){modified.resid=resid-(1-sm_prm)*thresh.resid}
    LevelNew=muhat+modified.resid
    
    SlopeNew=sm_prm2*(LevelNew-LevelPrev)+(1-sm_prm2)*SlopePrev
    
    if (flag_DEWMA){
      muhat=LevelNew+SlopeNew
    }else{
      muhat=muhat+modified.resid
    }
    
    yhat[i]=muhat*(1-b1-b2)+yhatINGARCH[i]     
    shift_hat[i]=muhat*(1-b1-b2)
    
    if (lambda0[i]>0 & yhat[i]>0&i>time_inj){
      weight=muhat*(1-b1-b2)  
      increment=weight*(y[i]*log(yhat[i]/lambda0[i])-(yhat[i]-lambda0[i]))
      #DES CUSUM
      S[i]=max(0,S[i-1]+increment)
    }else{
      S[i]=S[i-1]
    }
    
    LevelPrev=LevelNew
    SlopePrev=SlopeNew
    
  }
  # signal time
  RL=i
  
  return(list(S=S,RL=RL,shift_hat=shift_hat,yhat=yhat))
}
