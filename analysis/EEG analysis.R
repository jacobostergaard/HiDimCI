library(HiDimCI)
library(Matrix)
library(igraph)
library(misc)
misc::clean_up()

filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")

library(matlab)
library(R.matlab)


res = R.matlab::readMat("~/Downloads/eeg.set")

M = res$EEG[[16]]

plot_trial <- function(i=NULL){
  if(is.null(i)){
    Y = t(apply(M,c(1,2),mean))
  }else{
    Y = t(M[,,i])
  }
  Y = apply(Y,2,diff)
  matplot(Y, lty=1, type='l', col=1)
  return(Y)
}
Z = plot_trial()

# Hilbert transformation using Fourier transforms
hilbert <- function(inVec,n=length(inVec)){
  #if(!exists("n")) n = length(inVec)
  tmp = fft(inVec[1:n])
  hVec = rep(0,n)
  hVec[c(1,n/2+1)] = c(1,1)
  hVec[2:(n/2)] = 2
  outVec = fft(tmp*hVec,inverse=T)/n

  return(outVec[1:n])
}

# Unwrap phases function
unwrapPhase <- function(X,Y){

  p = ncol(X) # no. of oscillators
  N = nrow(X) # no. of observations

  phi = matrix(0,nr=N,nc=p)
  for(i in 1:p){

    phi[,i] = signal::unwrap((atan2(Y[,i],X[,i])+2*pi)%%(2*pi))
  }
  phi = data.frame(phi)
  names(phi) = paste0("phi",1:p)
  invisible(phi)
}

X = Y = Z
for(i in 1:ncol(Z)){
  Y[,i] = Im(hilbert(Z[,i]))
}
Y
P = unwrapPhase(X,Y)
P = P-P$Cz
channels = unlist(res$EEG[[23]][seq(1,640,10)])
names(P) = channels
P = P[,-which(names(P)=="Cz")]
matplot(P, lty=1, type='l', col=1)

# str(res)
# res$EEG
# length(res$EEG[[23]])
vecm = johansen(t(P), r=40, Psi = TRUE)

Pi = vecm$alpha%*%t(vecm$beta)
Pi.sym = .5*(Pi+t(Pi))
plot_adj_matrix(Pi.sym)

cls = cluster(Pi.sym)$grps

cls$channels = channels[cls$idx]
plot_adj_matrix(Pi.sym[cls$idx,cls$idx])
names(P)
cls

plot()
Y[,1]

unwrapPhase(hilbert(Y[,1]))
hilbert(Y[,1])
plot(Z[,1])
