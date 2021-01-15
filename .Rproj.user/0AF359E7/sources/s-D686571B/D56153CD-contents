library(HiDimCI)
library(misc)
library(Matrix)
library(doParallel)
filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")

# setwd(filelib)
clean_up(filelib)
set.seed(1234)

run_and_save <- function(N,p,d){
  sig     = 1
  mc_iter = 300
  dt      = 1
  # p-floor(p/d)*d    # This gives: p-floor(p/d)*d independent processes

  B      = as.matrix(HiDimCI::Kuramoto.matrix(p,d))
  Brank  = rankMatrix(B)[1]
  X      = HiDimCI::sampleX(B,N,sig)

  scramble    = sample(1:p)
  unscramble  = sort(scramble, index.return=TRUE)$ix
  X       = X[scramble,]
  filename  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,".Rda")

  # if(!file.exists(filename)){

    sim = HiDimCI::bootstrap(X, mc_iter,1, verbose=TRUE, ncores = detectCores())
    # sim = list()
    sim$scramble = scramble
    sim$unscramble = unscramble
    sim$X = X
    sim$B = B

    save(sim, file = filename)
  # }else{
    # load(filename)
  # }
}

d=8
# system.time({
#   run_and_save(2000,100,8)
# })


for(N in c(2000,1000,500)){
  p=100
  run_and_save(N,p,d)
  p=50
  run_and_save(N,p,d)
}


notify("Bootstrap test runs done!")


