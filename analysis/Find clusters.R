library(HiDimCI)
library(misc)
library(Matrix)
library(doParallel)
library(igraph)
filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")

# setwd(filelib)
clean_up(filelib)
set.seed(1234)

run_and_save <- function(N,p,d, scramble=FALSE){
  # N       = 2000 # 500,1000,2000 or 5000
  # p       = 50
  # d       = 8
  sig     = 1
  mc_iter = 300
  dt      = 1
  # p-floor(p/d)*d    # This gives: p-floor(p/d)*d independent processes

  B      = as.matrix(HiDimCI::Kuramoto.matrix(p,d))
  Brank  = rankMatrix(B)[1]
  X      = HiDimCI::sampleX(B,N,sig)

  if(scramble){
    scramble_idx    = sample(1:p)
    unscramble_idx  = sort(scramble_idx, index.return=TRUE)$ix
    X       = X[scramble_idx,]
    # X       = X[unscramble_idx,]
    filename  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,"_scrambled.Rda")
  } else{
    scramble_idx    = 1:p
    unscramble_idx    = 1:p
    filename  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,".Rda")
  }



  # if(!file.exists(filename)){
  b.test = HiDimCI::bootstrap_parallel(X, mc_iter,1, verbose=TRUE, cls = detectCores())
  b.test$scramble_idx = scramble_idx
  b.test$X = X
  b.test$B = B
  save(b.test, file = filename)
  # }else{
  # load(filename)
  # }
}
load_or_run <- function(N,d,p,scrambled){
  if(scrambled){
    fn  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,"_scrambled.Rda")
  } else {
    fn  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,".Rda")
  }

  if(file.exists(fn)){
    load(fn,envir = .GlobalEnv)
  } else{
    run_and_save(N,p,d, scrambled)
    load(fn, evir = .GlobalEnv)
  }

}

load_or_run(N=2000, d=8, p=100, scrambled=TRUE)
# load_or_run(N=1000, d=8, p=100, scrambled=TRUE)

# for(N in c(2000,1000,500)){
#   run_and_save(N,p,d, TRUE)
#   run_and_save(N,p,d, FALSE)
# }

# notify("Bootstrap test runs done!")
b.test$rank = rankTest(b.test)
b.test$rank$r

p = nrow(b.test$X)
B = b.test$B
Brank = rankMatrix(B)
X = b.test$X
r = Brank[1] #b.test$rank$r
r = b.test$rank$r
VECM_est = johansen(X,r)

par(bty='n')
boxplot(x=t(b.test$boot), col=add.alpha('dodgerblue',.75), border=add.alpha('black',.35), pch=16, cex=0.5, xaxt='n', outline=FALSE, range=0.95, lwd=2)
lines(1:p,b.test$test, lwd=2, col=add.alpha('red',.75), type='b', pch=16, cex=0.75)
abline(v=Brank+1, lwd=2, lty=3, col=add.alpha('black',.75))
axis(1, at=pretty(c(0,p))+1,labels = pretty(c(0,p)))


# Unrestricted estimator
P.hat = VECM_est$alpha%*%t(VECM_est$beta)
# Project estimated P.hat onto a symmetric subspace
P.prj = low_rank(sym_project(P.hat),r)
# Symmetric estimator of Pi
P.sym =  P.sym.hat(X,r)

test_matrix(P.hat,B)
test_matrix(P.prj,B)
test_matrix(P.sym,B)
rankMatrix(P.sym)
# logL(X,P.hat)-logL(X,M)


# plot_adj_matrix(B)
# plot_adj_matrix(P.hat)
# plot_adj_matrix(P.prj)
# plot_adj_matrix(P.sym)

idx = sort(b.test$scramble_idx, index.return=TRUE)$ix
idx = order(b.test$scramble_idx)
# plot_adj_matrix(P.hat[idx,idx])
# plot_adj_matrix(P.prj[idx,idx])
# plot_adj_matrix(P.sym[idx,idx])


cls = cluster(P.sym)$grps

# How many are put in each group:
tbl = table(cls$cluster)

Z = matrix(NA, nr=max(cls$cluster), nc=max(tbl))
for(i in 1:nrow(Z)){
  tmp = sort(cls$idx[cls$cluster==i])
  Z[i,1:length(tmp)] = tmp
}
Z = Z[order(Z[,1]),] # order by first column

plot_adj_matrix(P.sym[cls$idx,cls$idx])

Y = matrix(idx[1:96],nr=12,nc=8, byrow = TRUE)
tmp = matrix(NA, nr=4, nc=ncol(Y))
tmp[,1] = idx[97:100]

Y = t(apply(Y,1,sort)) # sort rows
Y = rbind(Y,tmp)

Y = Y[order(Y[,1]),] # order by first column
Z
xtable::xtable(cbind(Y,Z))
xtable::xtable(Z)

# max(cluster(P.sym)$mod)


fn = paste0(filelib,"plots/P_recovered.pdf")
# pdf(file=fn, height=10, width = 10)
  plot_adj_matrix(P.sym[cls$idx,cls$idx])
# dev.off()
fn = paste0(filelib,"plots/P_unscrambled.pdf")
# pdf(file=fn, height=10, width = 10)
  plot_adj_matrix(P.sym[idx,idx])
# dev.off()


  M = P.sym
  tmp   = diag(M) # extract diagonal, unconnected processes should have approx. zero
  p     = length(tmp)
  r.idx = which(tmp > mean(tmp)+2*sd(tmp) | tmp < mean(tmp)-2*sd(tmp) ) # find outliers, i.e. entries approx zero

  G       = graph.adjacency(abs(M[-r.idx,-r.idx]), mode="undirected", weighted=TRUE)
  cl = igraph::fastgreedy.community(G)

  plot(cl$modularity  )
