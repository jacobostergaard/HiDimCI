library(misc)
misc::clean_up()

filelib = icloud_lib("GitHub/Source/R/packages/HiDimCI/")

sp500 <- read.csv(paste0(filelib,"data/all_stocks_5yr.csv"))
comps <- read.csv(paste0(filelib,"data/constituents.csv"))

names(comps) = tolower(names(comps))
names(sp500)[7] = "ticker"


sp500 = dplyr::left_join(sp500, comps, by=c("ticker"="symbol"))


head(sp500)
dim(sp500)
str(sp500)
# summary(sp500)

tmp = unique(sp500[,8:9])
table(tmp$sector)

# tmp = subset(sp500, sector%in%c("Financials","Industrials"))
tmp = subset(sp500, sector%in%c("Real Estate","Health Care"))

# tmp = sp500


nms = unique(tmp$ticker)
df = data.frame(unique(sp500$date))
names(df) = "date"
for(i in 1:length(nms)){
  df = dplyr::left_join(df, sp500[sp500$ticker==nms[i],c(1,5)], by="date" )
}
names(df) = c("date",nms)

dim(df)

df2 = cbind(df$date[-1],apply(df[,-1],2, function(x) round(diff(x),3) ))
df2 = as.data.frame(df2)
names(df2)[1] = "date"
dim(df2)


df = df[,!apply(apply(df,2,is.na),2,any)]
X  = t(log(as.matrix(df[,-1])))

system.time({
  B.test = HiDimCI::bootstrap(X=X,B=50, dt =1, verbose=TRUE, ncores=6)
})
rankTest(B.test)

M = B.test
p = length(M$test)
pVals = numeric(0)
for(i in 1:p){
  test.cdf = ecdf(M$boot[i,])
  pVals[i] = 1-test.cdf(M$test[i])
}
names(pVals) = paste0("r=",(0:(p-1)))
if(any(pVals<0.05)){
  r.est = min(which(pVals<0.05))-2
} else{
  r.est = p
}
names(r.est) = ""
out = list(r=r.est,pVal = pVals)


r = r.est
SP500.est = johansen(X,r)

M = SP500.est$alpha%*%t(SP500.est$beta)
# M = .5*(M+t(M))
plot_adj_matrix(M)




cls = cluster(M)$grps

cls$ticker = rownames(X)[cls$idx]
tmp = unique(sp500[,7:9])
cls = dplyr::left_join(cls, tmp, by="ticker" )

length(unique(cls$cluster))

unique(cls$sector)
aggregate(idx ~ cluster+sector, data = cls, length)
# tmp = subset(cls, sector=="Materials")
# length(unique(tmp$cluster))

G = graph.adjacency(abs(M), mode="undirected", weighted=TRUE)

E(G)$width <- E(G)$weight + min(E(G)$weight) + 1 # offset=1
plot(G)

# qgraph(M,edge.labels=FALSE)  #if you want the weights on the edges as well

cls$name


unique(cls[,c(2,5)])

cls[cls$cluster %in% c(2,3,4,5,6,7,8),]



#
# y  = log(as.matrix(df[,-1]))
# dy = apply(y,2, diff)
#
# # y = apply(df[,-1],2, diff)
# # dy = apply(y, 2, diff)
#
#
#
#
# S00 = S01 = S10 = S11 = 0
# N = nrow(dy)
# i=2
# for(i in 2:N){
#   S00 = S00+dy[i,]%*%t(dy[i,])
#   S01 = S01+dy[i,]%*%t(y[i-1,])
#   S10 = S10+y[i-1,]%*%t(dy[i,])
#   S11 = S11+y[i-1,]%*%t(y[i-1,])
# }
#
# S00 = S00/N
# S01 = S01/N
# S10 = S10/N
# S11 = S11/N
#
#
# # tmp = eigen(S11)
# # tmp$values
# Matrix::rankMatrix(S11)
# dim(S11)
# # plot(tmp$values)
#
# P.sym = .5*(S01%*%solve(S11)+solve(S11)%*%S10)
#
# nms = rownames(P.sym)
# tmp = unique(sp500[,c(7,9)])
# nms = as.data.frame(nms)
# names(nms) = "ticker"
#
# nms = dplyr::left_join(nms,tmp)
# idx = order(nms$sector)
# P.sym = P.sym[idx,idx]
#
# nms = nms[idx,]
#
# unique(nms$sector)
# length(idx)
#
# dim(P.sym)
# r = 450
# M = svd(P.sym, nu = r, nv=r)
# P.sym2 = M$u%*%diag(M$d[1:r])%*%t(M$v)
#
# idx = which(nms$sector %in% c("Energy","Real Estate"))
# length(idx)
# # image(P.sym)
# image(P.sym2[idx,idx])
#
