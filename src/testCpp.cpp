// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "misc.h"
#include "tools.h"
#include "johaTools.h"
#include "johansen.h"

using namespace std;
using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat test_vecm(arma::mat X, int r, arma::mat A, arma::mat B, double dt){

  /* Inputs:
  *    X     Time series observations: dimensions should be p-by-N => each row is a marginal time-series
  *    r     The rank to use for estimation of parameters
  *    A,B   Restrictions for alpha/beta matrices. For no restrictions, input identity matrices.
  *    dt    Time resolution.
  */

  // Overall properties of the system
  int N = X.n_cols-1;   // Number of observations (-1 since we look at the differenced process)
  int p = X.n_rows;     // Dimension of the system
  int s = B.n_cols;     // Restrictions for beta matrix

  // Find Z0, Z1 and Z2 for input data matrix X
  mat Z0 = zeros<mat>(p,N);
  mat Z1 = zeros<mat>(p,N);
  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n);  // The differenced process dX[t] (calculated as X[t+1]-X[t])
    Z1.col(n) = X.col(n);             // The lagged process X[t-1]
  }
  mat Z2 = ones<mat>(1,N);


  // Find moment matrices M_{ij}
  mat M22,M02,M12,M22_1;
  M22   = M(Z2,Z2);
  M02   = M(Z0,Z2);
  M12   = M(Z1,Z2);
  M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0,R1;
  R0 = Z0-M02*M22_1*Z2;
  R1 = Z1-M12*M22_1*Z2;

  // Find covariance matrices S_{ij}
  mat S00,S11,S01,S10,S00_1, S10S00S01, S11_B;
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);


  // For restricted A: correct residuals
  mat Ap = M_perp(A);
  mat Ab = M_bar(A);
  mat Rt0,Rt1,St00,St11,St01,St10,St00_1;


  // Find the matrix to solve for the eigenvalue problem
  mat solveMat;

  if(Ap.size() != 1){ // A.perp is not degenerate => correct for alpha restrictions
    Rt0       = R0-S00*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;
    Rt1       = R1-S10*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;
    St00      = S(Rt0,Rt0);
    St11      = S(Rt1,Rt1);
    St01      = S(Rt0,Rt1);
    St10      = S(Rt1,Rt0);
    St00_1    = inv(Ab.t()*S00*Ab);
    S10S00S01 = B.t()*St10*Ab*St00_1*Ab.t()*St01*B;

  } else{ // Solve for the eigenvalues (and vectors) as usual (including beta restrictions)
    St00      = S(R0,R0);
    St11      = S(R1,R1);
    St01      = S(R0,R1);
    St10      = S(R1,R0);
    St00_1    = inv(S00);
    S10S00S01 = B.t()*S10*S00_1*S01*B;
  }

  S11_B       = B.t()*St11*B;

  // solve for eigenvalues
  cx_vec eigval_cx;
  cx_mat eigvec_cx;

  // first for S11...
  eig_gen(eigval_cx, eigvec_cx, S11_B);
  // C++ function returns complex vectors/matrices, so extract real parts
  vec p_val = real(eigval_cx);
  mat W     = real(eigvec_cx);

  // find decomposition of S11 from diagonalization
  mat S11_5 = W*diagmat(1/sqrt(p_val))*W.t();

  // then for S11_5...S11_5...
  eig_gen(eigval_cx, eigvec_cx, S11_5*S10S00S01*S11_5);
  // C++ function returns complex vectors/matrices, so extract real parts
  vec l_val = real(eigval_cx);
  mat U     = real(eigvec_cx);

  // convert to the correct eigenvectors
  mat V = S11_5*U;

  // set the values into old placeholders
  vec eigval = l_val;
  mat eigvec = V;

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  l_val = eigval(sort_index(eigval,"descend"));

  solveMat    = inv(B.t()*St11*B)*S10S00S01;
  eig_gen(eigval_cx, eigvec_cx, solveMat);

  eigval = real(eigval_cx);
  eigvec = real(eigvec_cx);

  eigval = eigval(sort_index(eigval,"descend"));

  return eigval-l_val;
}
