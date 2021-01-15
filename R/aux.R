make_batch <- function(nrep, ncores){
  # Setup list of batches (individual B-samples to send to workers)
  b0 = ncores # desired number of batches
  b1 = nrep %% b0
  b2 = floor(nrep/b0)
  b  = c(rep(b2,b0)) # this is the list of smaller bootstrap samples to send to workers

  if(b1>0){ # not all B samples have been distributed in each batch, iterate until all have been placed
    btmp = b1
    i=1
    while(btmp>0){
      b[i] = b[i]+1 # place 1 sample in a batch and move on to the next one, until all are placed
      i = i+1
      btmp = btmp-1
    }
  }

  b = b[b>0]

  return(b)
}
