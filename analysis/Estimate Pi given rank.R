library(HiDimCI)
library(misc)
library(Matrix)
filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")
set_normal_plot_bg()
clean_up(filelib)
set.seed(1234)

N = 2000
p = 100
d = 8

filename  = paste0(filelib,"data/Bootstrap_N",N,"_p",p,"_d",d,".Rda")
load(filename)

idx   = sim$unscramble
r.hat = rankTest(sim)$r
# r.hat = rankMatrix(sim$B)
r.hat   = rankMatrix(sim$B)


M.hat = johansen(sim$X,r.hat)

# Unrestricted OLS estimator
  P.ols = S_moments(sim$X)$OLS
# Unrestricted Johansen estimator
  M.hat = johansen(sim$X,r.hat-10)
  P.hat72 = M.hat$alpha%*%t(M.hat$beta)
  M.hat = johansen(sim$X,r.hat+10)
  P.hat92 = M.hat$alpha%*%t(M.hat$beta)
  M.hat = johansen(sim$X,r.hat)
  P.hat = M.hat$alpha%*%t(M.hat$beta)
# Project estimated P.hat onto a symmetric subspace
  P.hat1 = project_lift(target=P.hat, rank=r.hat)
  P.hat172 = project_lift(target=P.hat, rank=r.hat-10)
  P.hat192 = project_lift(target=P.hat, rank=r.hat+10)
# Symmetric estimator of Pi
  P.ols1 = project_lift(target=P.ols, rank=r.hat)
  P.ols172 = project_lift(target=P.ols, rank=r.hat-10)
  P.ols192 = project_lift(target=P.ols, rank=r.hat+10)

  blank <- function(){
    plot(0,0, type='n', bty='n', ann=FALSE, axes=FALSE)
  }

  fn = paste0(filelib,"plots/estimates.pdf")

  pdf(file=fn, height=11, width=8)
  layout(matrix(c( 1, 2, 2, 3, 3, 4,
                   5, 5, 6, 6, 7, 7,
                   8, 8, 9, 9,10,10,
                   11,11,12,12,13,13,
                   14,15,16,16,17,18), nr=5, nc=6, byrow=TRUE), heights = c(3,3,3,3,1))

    par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    blank()
    plot_matrix(sim$B)
    plot_matrix(P.ols[idx,idx])
    blank()

    plot_matrix(P.hat72[idx,idx])
    plot_matrix(P.hat[idx,idx])
    plot_matrix(P.hat92[idx,idx])

    plot_matrix(P.hat172[idx,idx])
    plot_matrix(P.hat1[idx,idx])
    plot_matrix(P.hat192[idx,idx])

    plot_matrix(P.ols172[idx,idx])
    plot_matrix(P.ols1[idx,idx])
    plot_matrix(P.ols192[idx,idx])
    blank()
    blank()
    plot_matrix(P.ols[idx,idx], colbar = "horiz")
    blank()
    blank()
  dev.off()








  test_matrix(P.ols,sim$B)
  test_matrix(P.hat,sim$B)
  test_matrix(P.hat1,sim$B)
  test_matrix(P.ols1,sim$B)
