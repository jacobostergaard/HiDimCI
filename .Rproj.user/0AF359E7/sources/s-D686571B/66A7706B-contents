# Tools for simulation

# Sample a multi-dim process starting at 0.
sampleX <- function(B,N,sig){
  p  = ncol(B)
  X  = matrix(0,nr=p,nc=N)
  # mu = rep(1,p)
  # mu = -runif(p,min = .25, max=1.25)
  mu = 0
  dt = 1
  for(i in 2:N){
    dW   = dt*rnorm(p, sd=sig)
    dX   = (B%*%X[,i-1]-mu)*dt+dW
    X[,i] = X[,i-1]+dX
  }
  return(X)
}
