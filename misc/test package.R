library(HiDimCI)
library(misc)
library(Matrix)

clean_up(filelib)

set.seed(1234)

N       = 2000 # 500,1000,2000 or 5000
p       = 100
d       = 8
sig     = 1
mc_iter = 300
p-floor(p/d)*d # This gives: p-floor(p/d)*d independent processes

dt= 1
B = as.matrix(HiDimCI::Kuramoto.matrix(p,d))
Brank = rankMatrix(B)[1]
e = dt*matrix(rnorm(N*p),nr=p,nc=N)

X   = HiDimCI::sampleVAR1(N,rep(0,p),B+diag(p),rep(0,p),sig*diag(p),e,dt)


b100 = Kuramoto.beta(p,rep(d,floor(p/d)),rep(d,floor(p/d)))
b100 = diag(p)
# johansen(X,r = Brank,H = b100)
M = johansen(X,r = Brank)

plot(79:89,M$test[80:90], type='l', bty='n')

dX = diff(t(X))
X1 = t(X[,-N])
R0 = t(dX-apply(dX,2,sum)/N)
R1 = t(X1-apply(X1,2,sum)/N)

S00 = Sij(R0,R0)
S01 = Sij(R0,R1)
S10 = Sij(R1,R0)
S11 = Sij(R1,R1)
S00_1 = solve(S00)

S10S00S01   = S10%*%S00_1%*%S01
S11.b       = t(b100)%*%S11%*%b100
S10S00S01.b = t(b100)%*%S10S00S01%*%b100
e1    = eigen(S11.b)$values
W     = eigen(S11.b)$vectors
S11.5 = W%*%diag(1/sqrt(e1))%*%t(W)
eigen(S11.5%*%S10S00S01.b%*%S11.5)
l     = sort(eigen(S11.5%*%S10S00S01.b%*%S11.5)$values,decreasing=TRUE)
U     = eigen(S11.5%*%S10S00S01.b%*%S11.5)$vectors
V     = S11.5%*%U

