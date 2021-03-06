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

// Declare functions
arma::mat sampleVAR1(int N, arma::vec X0, arma::mat A, arma::vec mu, arma::mat D, arma::mat e, double dt);


// Bootstrapping loop, using the Johansen procedure for estimation of statistic-distribution
arma::mat bootLoop(arma::mat X, int r, int B, double dt, bool verbose){
  /* Inputs are
  *    X     Time series observations: dimensions should be p-by-N => each row is a marginal time-series
  *    r     The rank to use for estimation of parameters
  *    B     Number of bootstrap samples
  *    dt    Time resolution.
  */

  bool estPsi = FALSE;

  // Overall properties of the system
  int N = X.n_cols-1;   // Number of observations
  int p = X.n_rows;     // Dimension of the system

  double pct = (double) r/p;

  mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot, A, e, Mboot;
  mat out = zeros<mat>(1,B);
  mat Ip = eye<mat>(p,p);
  vec X0 = zeros<vec>(p);

  // Get residuals from model estimation
  if(r == 0){
    Minit   = var(X, estPsi, dt);
    Pi      = zeros<mat>(p,p);
  } else{

    Minit   = vecm(X,r,Ip,Ip, estPsi, dt);
    alpha   = Minit.cols( 0         , r-1       );
    beta    = Minit.cols( r         , 2*r-1     );
    Pi      = alpha*beta.t();
  }
  res     = Minit.cols( 2*r+p+2+1 , 2*r+p+2+N );
  Psi     = Minit.cols( 2*r       , 2*r       );
  S2      = Minit.cols( 2*r+1     , 2*r+p   );


  // Center residuals
  res = (res-repmat(mean(res,1),1,res.n_cols));
  // Run bootstrap samples, each sample is called Xboot, each model is Mboot
  for(int b=0;b<B;b++){
    e = res%randn(p,N);   // Wild bootstrap innovations
    A = Pi+Ip;            // VAR matrix

    Xboot = sampleVAR1(N, X0, A, Psi, S2, e, dt);     // Sample a VAR(1) with parameter estimates
    Mboot = vecm(Xboot,1,Ip,Ip, estPsi, dt);                // Find test statistics
    // Mboot = johansenCpp(Xboot,r,Ip,Ip, estPsi, dt);                // Find test statistics

    out(0,b) =  Mboot(r,2*1+p+1);
    if(verbose){
      display_progress(pct + (double) (b+1)/(B*p) );
    }
  }

  return out;
}



// The the function returns ...
// [[Rcpp::export]]
arma::mat bootstrapCpp(arma::mat X, int B, double dt, bool verbose){

  // Include support for restricted alpha/beta matrices!
  int p = X.n_rows;
  mat est,boot,out;

  if(B > 0){
    // Run bootstrap method to determine test-statistic distributions!
    out = zeros<mat>(p,B);
    for(int r=0;r<p;r++){
      boot = bootLoop(X,r,B,dt, verbose);
      out.row(r) = boot;
          // display_progress((double) (r+1)/p);
    }
  } else{
    out = zeros<mat>(p,1);
  }
  return out;
}
