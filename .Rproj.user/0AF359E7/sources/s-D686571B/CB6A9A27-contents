Z

make_adj_matrix <- function(Z){

  tmp     = as.numeric(t(Z))
  ord     = tmp[!is.na(tmp)]
  inv_ord = sort(ord, index.return=TRUE)$ix
  p = max(Z, na.rm = TRUE)
  A = matrix(0, nr=p, nc=p)
  for(i in 1:p){
    idx = which(Z==i, arr.ind = TRUE)
    A[i,Z[idx[1],]] = 1
  }
  A  = A[ord,ord]
  ps = apply(Z, 1, function(x) sum(!is.na(x)))
  rs = ps-1

  diag(A) = rep(-rs,ps)
  A = A[inv_ord,inv_ord]

  out = list(A=A, p=ps, ord=ord, inv_ord = inv_ord)
  return(out)
}
tmp = make_adj_matrix(Z)

image(t(A[tmp$ord,tmp$ord]))

O = (1:16)/16
O[which(tmp$p == 1)] = 0
tmp2 = rep(O, tmp$p)
B = A%*%diag(tmp2[tmp$inv_ord])

image(t(B))
image(t(B[tmp$ord,tmp$ord]))

cls_tmp = rep(1:length(tmp$p),tmp$p)
M_tmp = P.sym[tmp$ord, tmp$ord]
nd = length(tmp$p)
O  = numeric(nd)
for(i in 1:nd){
  i_tmp  = which(cls_tmp==i)
  M_tmp2 = M_tmp[i_tmp, i_tmp]
  O[i] = (sum(M_tmp2)-sum(diag(M_tmp2)))/(length(i_tmp)^2-length(i_tmp))

}
O[is.infinite(O)] = 0
O

f <- function(O){
  O[which(tmp$p == 1)] = 0
  tmp2 = rep(O, tmp$p)
  B = A%*%diag(tmp2[tmp$inv_ord])
  return(-logL(X,B)  )
}

f(O)
O2 = O
O2[O2>0] = 0.02040816
f(O2)
O3 = O2
O3[3] = O2[3]*.92
f(O3)

system.time({
  res = nlm(f,O)
})

res$par
res$value
