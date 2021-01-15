// Bootstrapping for critical values for Johansens rank test for a VECM
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional...

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"
#include "misc.h"
#include "johansen.h"

using namespace Rcpp;
using namespace arma;
using namespace std;




// The the function returns ...
// [[Rcpp::export]]
arma::mat sampleVAR1(int N, arma::vec X0, arma::mat A, arma::vec mu, arma::mat D, arma::mat e, double dt){
  /* Inputs are
  *    N    Number of timesteps (of size dt)
  *    X0   Initial values (vector)
  *    A    Coefficient matrix for lagged values
  *    mu   Vector of constant trenda
  *    D    Cholesky decomp of covariance matrix, such that DD' is the covariance matrix
  *    e    Errors
  *    dt   Time resolution.
  */

  int p = X0.n_rows;        // Dimension of the system
  mat X = zeros<mat>(p,N);  // Output

  X.col(0) = X0;

  for(int n=1;n<N;n++){
    X.col(n) = (A*X.col(n-1)+mu)*dt + D*e.col(n);
  }

  return X;
}
