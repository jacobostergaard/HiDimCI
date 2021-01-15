# Different variations of the 'B' matrix

Kuramoto.matrix <- function(p,d, sparse=FALSE, r.sparse=NULL, scramble=FALSE){
  # A Kuramoto type block matrix where 'p' is the desired dimension of the matrix
  # and 'd' the desired size of each block.
  # Within each block, the diagonal value corresponds to the sum of the off-diagonal row entries,
  # and with opposite signs, such that row sums are always zero.
  # When p mod d is non zero, the matrix is padded with zeros in the lower right to be p-times-p.
  # The input "scramble" will randomly interchange rows.

  if(sum(d)>p){
    p = sum(d)
  }

  out = matrix(0,nr=p,nc=p)

  if(length(d)>1){
    out = matrix(NA,nr=0,nc=0)
    for(j in 1:length(d)){
      out = Matrix::bdiag(out,Kuramoto.matrix(d[j],d[j]))
    }
  } else{
    if(d>1){
      a = 1/(d-1)
      b = a/(d-1)
      tmp = numeric(d)
      Tmp = matrix(nr=d,nc=d)
      for(i in 1:d){
        tmp[i] = -a
        tmp[-i] = b
        Tmp[i,] = tmp
      }
    } else{
      Tmp = as.matrix(0)
    }
    out = Tmp
    if(d<p/2){
      for(i in 1:(floor(p/d)-1))
        out = Matrix::bdiag(out,Tmp)
    }
  }

  if(ncol(out)!=p){
    for(i in 1:(p-ncol(out))){
      out = Matrix::bdiag(out,0)
    }
  }

  if(sparse){
    tmp.sparse = Matrix::rsparsematrix(p, p, nnz = r.sparse, symmetric=TRUE)
    out[tmp.sparse==0] = 0
  }

  if(scramble){
    out = out[sample(nrow(out)),]
  }

  return(out)
}


bimodal.matrix <- function(p,r.sparse,mu1=5, mu2=-5,sig=1, sparse=TRUE){
  # This function returns a sparse matrix where the entries are bimodal normal distributed
  rbimodal <- function(n,mu1=5,mu2=-5, sig=1){
    tmp = rnorm(n)
    for(i in 1:n){
      if(sign(tmp[i])>0){
        tmp[i] = rnorm(1,mean = mu1, sd=sig)
      } else if(sign(tmp[i])<0){
        tmp[i] = rnorm(1,mean = mu2, sd=sig)
      } else{
        tmp[i] = 0
      }
    }
    return(tmp)
  }

  if(sparse)
    out = Matrix::rsparsematrix(p, p, nnz = r.sparse,  rand.x = rbimodal, symmetric=TRUE)
  return(out)
}

