W.dens <- function(x, mu, beta){
  u = mu
  b = beta
  A = abs(sqrt(u-u*b)-sqrt(b-u*b))
  B = abs(sqrt(u-u*b)+sqrt(b-u*b))
  fx = numeric(length(x))
  for(i in 1:length(x)){
    if( (x[i]>0) & (x[i]<1)){
      fx[i] = suppressWarnings(sqrt((x[i]-A)*(x[i]+A)*(B-x[i])*(B+x[i]))/(pi*b*x[i]*(1-x[i])*(1+x[i])))    
    }else if(x[i]==0){
      fx[i] = max(0,1-u/b) 
    }else if(x[i]==1){
      fx[i] = max(0,1-(1-u)/b)
    }
  }
  fx[is.nan(fx)] = 0
  return(fx)
}

K <- function(a,b,x){
  return(1/(pi*b)*(x*sqrt(a-x^2*U(a,b,x)^2)-a*acos(x*U(a,b,x)/sqrt(a))+b*acos(x*(U(a,b,x)-1)/sqrt(b))))
}

U <- function(a,b,x){
  return(0.5+0.5*(a-b)/x^2)
}

W.dist <- function(x,mu,beta){
  u = mu
  b = beta
  tmp = suppressWarnings((1-u)*K(u-u*b,b-u*b,x)+u*(1-K(1-b-u+u*b,u*b,sqrt(1-x^2))))
  if(all(is.nan(tmp))){
    y = seq(0,1,1e-4)
    tmp2 = suppressWarnings((1-u)*K(u-u*b,b-u*b,y)+u*(1-K(1-b-u+u*b,u*b,sqrt(1-y^2))))
    tmp2[is.nan(tmp2)] = 0  
    idx = which(diff(tmp2)< -1e-3)
    if(length(idx)>0) tmp2[min(idx+1):length(y)] = tmp2[min(idx)]
    tmp2[length(y)] = 1
    tmp = numeric(length(x))
    for(i in 1:length(x)){
      ida = max(which(y < x[i]))
      idb = min(which(y >= x[i]))
      a = y[ida]
      b = y[idb]
      w = (b-x[i])/(b-a)
      tmp[i] = w*tmp2[ida]+(1-w)*tmp2[idb]  
    }
  }else if(any(is.nan(tmp))){
    tmp[is.nan(tmp)] = 0  
  }
  
  idx = which(diff(tmp)< -1e-3)
  if(length(idx)>0) tmp[min(idx+1):length(x)] = tmp[min(idx)]
  
  # if(max(x==1)) tmp[length(x)] = 1
  return(tmp)
}



W.quan <- function(q,mu,beta, res=1e-4){
  
  res=1e-4
  x = seq(0,1,res)  
  tmp = numeric(length(q))
  for(i in 1:length(q)){
    Wtmp = W.dist(x,mu,beta)
    if(any(Wtmp<q[i])){
      ida = max(which(Wtmp < q[i]))
    }else{
      ida = 1
    }
    if(any(Wtmp>=q[i])){
      idb = min(which(Wtmp >= q[i]))  
    } else{
      idb = min(which(Wtmp==max(Wtmp)))
    }
    a = Wtmp[ida]
    b = Wtmp[idb]
  
    if(b>a){
      w = (b-q[i])/(b-a)
      tmp[i] = w*x[ida]+(1-w)*x[idb]  
    }else{
      tmp[i] = x[idb]
    }
    
  }
  return(tmp)
}
