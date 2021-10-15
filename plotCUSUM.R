plotCUSUM<-function(time_inj,b0,b1,b2,b_cos,b_sin,P,dt,beta_c,y,h_s){
  
  i=1
  S=0
  lambda0=b0/(1-b1-b2)
  lambda1_c=lambda0
  while (S[i]<h_s|i<=time_inj){
    i=i+1
    k=i-1
    
    # step shift cusum
    delta1_c=b0+beta_c
    # baseline rate
    lambda0[i]=b0+b1*y[i-1]+b2*lambda0[i-1]+b_cos*cos(2*pi/P*i)+b_sin*sin(2*pi/P*i)
    lambda1_c[i]=delta1_c+b1*y[i-1]+b2*lambda1_c[i-1]+b_cos*cos(2*pi/P*i)+b_sin*sin(2*pi/P*i)
    
    if (lambda0[i]>0 &i>time_inj){
      S[i]=max(0,S[i-1]+y[i]*log(lambda1_c[i]/lambda0[i])-(lambda1_c[i]-lambda0[i]))
    }else{
      S[i]=S[i-1]
    }
    
  }
  
  
  # signal time
  RL=i
  return(list(S=S,RL=RL))
}
