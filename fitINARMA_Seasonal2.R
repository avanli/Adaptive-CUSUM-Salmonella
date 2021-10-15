fitINARMA_Seasonal2<-function(y,P,iter_max){
  #P: period
  N=length(y)
  yr=y[2:N]
  yl=y[1:N-1]
  y_hat=matrix(NA,N-1,1)
  
  t=seq(1,N-1,1)
  
  m1 <- glm(yr ~cos(2*pi*t/P)+sin(2*pi*t/P)+yl+t , family=poisson(link=identity),start=c(1,0,0,0,0))
  
  d=m1$coefficients[1]
  b_cos=m1$coefficients[2]
  b_sin=m1$coefficients[3]
  phi=m1$coefficients[4]
  dt=m1$coefficients[5]
  theta=0
  
  #d=d0
  #phi=phi0
  #theta=theta0
  
  
  iter=2
  while (iter<=iter_max)
  {
    eps_hat=0
    #INARMA(1,1)
    for (i in 2:N){
      eps_hat[i]=y[i]-phi[iter-1]*y[i-1]+theta[iter-1]*eps_hat[i-1]-d[iter-1]-dt[iter-1]*(i-1)-
        b_cos[iter-1]*cos(2*pi*(i-1)/P)-b_sin[iter-1]*sin(2*pi*(i-1)/P)
    }
    eps_hat_l=eps_hat[1:N-1] 
    m2<- glm(yr ~cos(2*pi*t/P)+sin(2*pi*t/P)+yl+eps_hat_l+t , family=poisson(link=identity),start=c(1,0,0,0,0,0))
    
    
    d[iter]=m2$coefficients[1]
    b_cos[iter]=m2$coefficients[2]
    b_sin[iter]=m2$coefficients[3]
    phi[iter]=m2$coefficients[4]
    theta[iter]=-m2$coefficients[5]
    dt[iter]=m2$coefficients[6]
    
    iter=iter+1
  }
  
  # one-step ahead forecasts for time periods 2,3,...,N
  y_hat=dt[iter_max]*t+d[iter_max]+b_cos[iter_max]*cos(2*pi*t/P)+
    b_sin[iter_max]*sin(2*pi*t/P)+phi[iter_max]*yl+theta[iter_max]*eps_hat_l
  # Mean square error
  MSE=(y[2:N])%*%(y_hat)
  
  # return fitted model and coefficients
  result=list(model=m2,d=d[iter_max],phi=phi[iter_max],theta=theta[iter_max], dt=dt[iter_max],fitted=y_hat,MSE=MSE)
  return(result)
  
}
