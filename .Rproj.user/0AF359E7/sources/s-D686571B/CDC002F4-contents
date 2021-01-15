library(HiDimCI)
library(misc)
library(Matrix)
library(doParallel)
filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")
clean_up(filelib)
set.seed(1234)

p  = 100
d  = 8
Ns = c(500,1000,2000,5000)
Ss = c(.5,1,2)
angles = list()
for(i in 1:3){
  tmp = list()
  for(j in 1:3){
    N   = Ns[i]
    sig = Ss[j]
    fn  = paste0(filelib,"data/angles2_N",N,"_p",p,"_d",d,"_sig",sig,".Rda")
    fn  = sub("\\.","",fn)
    load(fn)
    tmp[[j]] = angleOLS
  }
  names(tmp) = paste0("s",Ss)
  angles[[i]] = tmp
}
names(angles) = paste0("N",Ns)

layout(1)
cols = c('orange','steelblue','tomato2')
plot(0,0, type='n', xlim=c(0,100), ylim=c(0,pi/2), ann=FALSE, bty='n')
for(i in 1:3){
  for(j in 1:3){
    tmp = angles[[i]][[j]]
    tmp = tmp[1:p,]
    tmp = apply(tmp,1, function(x) quantile(x, probs=c(.025,.5,.975)))
    polygon(c(1:100,100:1), c(tmp[1,],rev(tmp[3,])), border=NA, col=add.alpha(cols[i],.25))
    # lines(tmp[2,], col=cols[i],lwd=2)
  }
}
abline(v=84, lty=3, lwd=2)
legend("bottomleft",c(paste0("N=",c(Ns)), expression("True"~B),"True rank"), bty='n', col=c(cols,'black','black'), lty=c(1,1,1,1,3), lwd=2)
mtext("Estimation rank")
mtext("Angle between")

B = as.matrix(HiDimCI::Kuramoto.matrix(p,d))
x = y = 1:p
for(i in 1:length(x)){
  ri = x[i]
  y[i] = acos(matrix_angle(low_rank(B,ri),B))
}
lines(x,y, lwd=2, col=add.alpha('black',.5))
