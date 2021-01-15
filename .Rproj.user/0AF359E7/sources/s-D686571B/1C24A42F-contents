// Perform Johansen procedure for alpha/beta-restricted VECM models.
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional... for now.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "johaTools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::mat vecm(arma::mat X, int r, arma::mat A, arma::mat B, bool Psi, double dt){

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
    if(Psi){
      R0 = Z0-M02*M22_1*Z2;
      R1 = Z1-M12*M22_1*Z2;
    } else{
      R0 = Z0;
      R1 = Z1;
    }


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

  // solveMat    = inv(B.t()*St11*B)*S10S00S01;
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
    eigval = eigval(sort_index(eigval,"descend"));




  // Normalize eigenvectors
    for(int i=0;i<s;i++){

      if(Ap.size() != 1){ // A.perp is not degenerate!
        double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*St11*B)*eigvec.col(i)));
        eigvec.col(i) = eigvec.col(i)/sqrt(nf);
      } else{
        double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*S11*B)*eigvec.col(i)));
        eigvec.col(i) = eigvec.col(i)/sqrt(nf);
      }
  }

  // To use cumsum for the teststats, the eigenvalues must be sorted "in reverse", this is ascending...
    vec testStat = -N*cumsum(log(1-eigval(sort_index(eigval,"ascend"))));
    testStat = sort(testStat,"descend");

  // Normalize beta with the identity matrix of size 'r' in the upper 'r' rows using c {p x r} = (I{r x r},0{()p-r) x r})
    mat b_hat = B*eigvec.cols(0,r-1);
    mat Ir = eye<mat>(r,r);
    mat c = join_cols(Ir,zeros<mat>(p-r,r));
    if(det(c.t()*b_hat)!=0){ // Normalize b_hat if possible...
      b_hat = b_hat*inv(c.t()*b_hat);
    }

// mat b_print = b_hat.submat(0,0,6,6);
// cout << round(100*b_print)/100 << endl;



  // Find OLS estimates of alpha (loadings), psi (deterministic trends) and omega (covariance)
    mat BS11B_1 = inv(b_hat.t()*S11*b_hat);
    mat a_hat = A*Ab.t()*S01*b_hat*BS11B_1;
    mat Psi_hat = (M02*M22_1-a_hat*b_hat.t()*M12*M22_1);
    if(!Psi){
      // Set estimate at 0 if not to be estimated...
      Psi_hat = 0*Psi_hat;
    }
    mat Omega_hat = S00-S01*b_hat*BS11B_1*b_hat.t()*S10;

  // Calculate residuals
    mat Pi_hat = a_hat*b_hat.t();
    mat res = zeros<mat>(p,N);
    for(int n=0;n<N;n++){
      res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n))-Psi_hat;
    }

  // Set the output with parameter estimates (alpha, beta, psi, omega), test statistics, eigenvalues and residuals
    int outRows = p;
    int outCols = a_hat.n_cols+b_hat.n_cols+1+p+2+N;
    mat out = zeros<mat>(outRows,outCols);

  // Insert alpha estimate in output
    for(int i=0; i<a_hat.n_cols;i++){
      out.col(i) = a_hat.col(i)/dt;
    }

  // Insert beta estimate in output
    for(int i=0;i<b_hat.n_cols;i++){
      int j = a_hat.n_cols;
      out.col(j+i) = b_hat.col(i);
    }

  // Insert Psi estimate in output
    int j = a_hat.n_cols+b_hat.n_cols;
    out.col(j) = Psi_hat/dt;

  // Insert Omega estimate in output
    for(int i=0;i<p;i++){
      j = a_hat.n_cols+b_hat.n_cols+1;
      out.col(j+i) = Omega_hat.col(i)/dt;
    }

  // Insert test statistic in output, append zeros to match dimensions...
    testStat = join_cols(testStat,zeros<mat>(p-s,1));
    eigval = join_cols(eigval,zeros<mat>(p-s,1));

    j = a_hat.n_cols+b_hat.n_cols+1+p;
    out.col(j) = testStat;
    j = a_hat.n_cols+b_hat.n_cols+1+p+1;
    out.col(j) = eigval;

  // Insert residuals in output
    for(int i=0;i<N;i++){
      j = a_hat.n_cols+b_hat.n_cols+1+p+2;
      out.col(j+i) = res.col(i);
    }
  return out;
}


// Estimation for a VAR model when r=0...
arma::mat var(arma::mat X, bool Psi, double dt){

  /* Inputs:
   *    X     Time series observations: dimensions should be p-by-N => each row is a marginal time-series
   *    dt    Time resolution.
   */

  // Overall properties of the system
    int N = X.n_cols-1;
    int p = X.n_rows;

  // Setup placeholder for difference process and deterministic trends
    mat Z0 = zeros<mat>(p,N);
    mat Psi_hat = zeros<mat>(p,1);

  // Estimation using Least Squares, formula is: Psi_hat = T^-1*sum_{t=1}^T dX_t (this is px1 dimensional)
    if(Psi){
      // Estimate if given as TRUE
      for(int n=0;n<N;n++){
        Z0.col(n) = X.col(n+1)-X.col(n);
        Psi_hat += Z0.col(n);
      }
      Psi_hat = Psi_hat/N;
    } else{
      for(int n=0;n<N;n++){
        Z0.col(n) = X.col(n+1)-X.col(n);
      }
    }



  // Calculate residuals and covariance estimator (see Lütkepohls book, p. 75)
    mat res = zeros<mat>(p,N);
    mat Omega = zeros<mat>(p,p);
    for(int n=0;n<N;n++){
      res.col(n) = Z0.col(n)-Psi_hat;
      Omega += res.col(n)*res.col(n).t();
    }
    Omega = (Omega/N)*(N/(N-1));

    // cout << "res" << endl;
    // cout << res << endl;

  // Calculate r-hypotheses statistics
    int r_tmp = 1;
    mat tmp = vecm(X,r_tmp, eye<mat>(p,p), eye<mat>(p,p), Psi, dt);
    mat test  = tmp.cols(2*r_tmp+p+1,2*r_tmp+p+1);
    mat eigs  = tmp.cols(2*r_tmp+p+2,2*r_tmp+p+2);
    mat joha = join_rows(test,eigs);

  // Model estimates
    mat est = join_rows(Psi_hat/dt,Omega/dt);

  // Add test statistics for r- hypotheses
    mat est_test = join_rows(est,joha);

  // Output estimates and residuals
    mat out  = join_rows(est_test,res);

  // Estimation using MLE is not implemented yet.

  return out;
}




// [[Rcpp::export]]
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, bool Psi, double dt=1){
// Johansen estimation procedure, returns parameter estimate, test statistics, eigenvalues and residuals

  int N = X.n_cols-1;   // Number of observations
  int p = X.n_rows;     // Dimension of the system

  mat out = zeros<mat>(p,N);
  if(r > 0){
    out = vecm(X,r,A,B,Psi,dt);
  } else{
    out = var(X,Psi,dt);
  }
  return out;
}
