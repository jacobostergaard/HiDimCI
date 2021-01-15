# Wrapper functions for simulating data and performing bootstrap test and inference (restrictions on matrices)
# Output here is more user friendly and useable.


johansen <-function(X,r=1,A=NULL,H=NULL, Psi = FALSE, dt=1){
  # Wrapper function for the johansen cpp function.

  X = as.matrix(X)
  N = ncol(X)-1
  p = nrow(X)

  if(is.null(A)){
    A = as.matrix(diag(p))
  }
  if(is.null(H)){
    H = as.matrix(diag(p))
  }
  df = r*(p-ncol(A))+r*(p-ncol(H))

  out = HiDimCI::johansenCpp(X,r,as.matrix(A),as.matrix(H), Psi, dt)

  if(r > 0){
    a.hat = matrix(out[,1:r],nr=p)         # alpha estimate
    b.hat = matrix(out[,(r+1):(2*r)],nr=p) # beta estimate
    rownames(a.hat) = paste0("x",1:p)
    colnames(a.hat) = paste0("r",1:r)
    rownames(b.hat) = paste0("x",1:p)
    colnames(b.hat) = paste0("r",1:r)
  }
  P.hat = out[,2*r+1]                                   # deterministic trend estimate
  O.hat = out[,(2*r+2):(2*r+1+p)]                       # omega estimate
  test  = out[,(2*r+2+p)]                               # test statistics
  eigs  = out[,(2*r+2+p+1)]                             # eigenvalues
  res   = out[,(2*r+2+p+2):(2*r+2+p+N+1)]               # residuals

  names(P.hat)    = paste0("x",1:p)
  rownames(O.hat) = paste0("x",1:p)
  colnames(O.hat) = paste0("x",1:p)
  names(test)     = paste0("r=",0:(p-1))
  names(eigs)     = paste0("l",1:p)
  rownames(res)   = paste0("x",1:p)

  if(r==0){
    a.hat = b.hat = NULL
  }
  return(list(N=N, p=p, r=r ,alpha = a.hat, beta = b.hat, Psi=P.hat, Omega = O.hat, test=test, lambda=eigs, A=A, H=H, df=df, dt=dt, res=res, data=X))
}



bootstrap <-function(X,B=1000, dt =1, verbose=TRUE, ncores=1){
  # Wrapper function for the boostrap cpp function. Default is 1000 bootstrap samples.

  X = as.matrix(X)
  N = ncol(X)-1
  p = nrow(X)
  r = 1
  tmp = johansenCpp(X,r,diag(p),diag(p), FALSE, dt)
  out = list()
  out$test = tmp[,(2*r+2+p)]

  if(ncores>1){

    # Setup list of batches (individual B-samples to send to workers)
    b0 = ncores # desired number of batches

    b1 = B %% b0
    b2 = floor(B/b0)
    b  = c(rep(b2,b0)) # this is the list of smaller bootstrap samples to send to workers

    if(b1>0){ # not all B samples have been distributed in each batch, iterate until all have been placed
      btmp = b1
      i=1
      while(btmp>0){
        b[i] = b[i]+1 # place 1 sample in a batch and move on to the next one, until all are placed
        i = i+1
        btmp = btmp-1
      }
    }

    # Parallel part:
    cl= makeCluster(ncores, setup_timeout=.5)
    # parallel::clusterExport(cl, c("bootstrapCpp","X","dt"))
    parallel::clusterExport(cl, c("bootstrapCpp"))
    tmp = parallel::parLapply(cl, b, function(i){
                                          return(bootstrapCpp(X,B=i,dt, verbose=FALSE))
                                        })
    stopCluster(cl)


    # Extract results
    Mtmp = numeric(0)
    for(i in 1:length(tmp)){
      Mtmp = cbind(Mtmp,tmp[[i]])
    }

    out$boot = Mtmp
  } else{
    out$boot = bootstrapCpp(X,B=B,dt, verbose)
  }

  return(out)
}


