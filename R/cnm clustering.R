
cluster <- function(M){
  tmp   = diag(M) # extract diagonal, unconnected processes should have approx. zero
  p     = length(tmp)

  # find outliers, i.e. entries approx zero using "leave-one-out" estimation of the mean, analogous to the dfbetas measure
  m = numeric(length(tmp))
  for(i in 1:length(tmp)){
    m[i] = mean(tmp)-mean(tmp[-i])
  }
  r.idx = which(m > 3*sd(m) |  m < -3*sd(m))
  # r.idx = which(tmp > mean(tmp)+2*sd(tmp) | tmp < mean(tmp)-2*sd(tmp) )

  G       = graph.adjacency(abs(M[-r.idx,-r.idx]), mode="undirected", weighted=TRUE)
  cluster = igraph::fastgreedy.community(G)

  grps = data.frame( idx = c( (1:p)[-r.idx], r.idx ), cluster= c(cluster$membership, max(cluster$membership)+1:length(r.idx)) )

  grps = grps[order(grps$cluster, grps$idx),]

  return(list(grps=grps, mod=cluster$modularity))
}

