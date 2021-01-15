library(HiDimCI)
library(misc)
library(Matrix)
library(doParallel)
filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")
clean_up(filelib)
set.seed(1234)

make_angles <- function(X, P, verbose=FALSE){
  p    = nrow(X)
  out  = numeric(p)
  for(i in 1:p){
    if(verbose){
      cat("\rr=",i, sep="")
    }
    M      = HiDimCI::johansen(X,i)
    Phat   = M$alpha%*%t(M$beta)
    out[i] = acos(HiDimCI::matrix_angle(Phat,P))
  }
  return(out)
}

Ns  = c(500,1000,2000,5000)
Ss  = c(.5,1,2)
for(i in 1:3){
  for(j in 1:3){
    nrep = 50
    b    = make_batch(nrep,ncores=8)
    N    = Ns[i]
    p    = 100
    d    = 8
    sig  = Ss[j]

    f <- function(N,p,d,sig){
      B    = as.matrix(HiDimCI::Kuramoto.matrix(p,d))
      X    = HiDimCI::sampleX(B,N,sig)
      Pols = S_moments(X)$OLS
      res  = make_angles(X, Pols, FALSE)
      return(res)
    }

    cat("\nN:",N,"sig:",sig,"\n")

    # Parallel part:
    cl= makeCluster(length(b), setup_timeout=.5)
    parallel::clusterExport(cl, c("f","N","p","d","sig","S_moments", "make_angles"))
    tmp = parallel::parLapply(cl, b, function(k){
      out = replicate(k, f(N,p,d, sig))
      return(out)
    })
    stopCluster(cl)

    Mtmp = numeric(0)
    for(k in 1:length(tmp)){
      Mtmp = cbind(Mtmp,tmp[[k]])
    }
    angleOLS = rbind(pi/2,Mtmp)

    filename  = paste0(filelib,"data/angles2_N",N,"_p",p,"_d",d,"_sig",sig,".Rda")
    filename = sub("\\.","",filename)
    save(angleOLS, file=filename)

    msg = paste0("N=",N," sig=",sig," done!")
    notify(msg)
  }
}

notify("All angles done!")



