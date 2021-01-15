#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "johaTools.h"


// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

// Declare functions

// [[Rcpp::export]]
arma::mat S_eigen(arma::mat X){
    // Find Z0, Z1 and Z2 for input data matrix X
    int p = X.n_rows;     // Dimension of the system
    int N = X.n_cols-1;   // Number of observations

    mat Z0 = zeros<mat>(p,N);
    mat Z1 = zeros<mat>(p,N);
    mat Z2 = ones<mat>(1,N);
    // mat Z2 = mu;

    for(int n=0;n<N;n++){
      Z0.col(n) = X.col(n+1)-X.col(n);  // This is dX
      Z1.col(n) = X.col(n);             // This is X[t-1]
    }

    // Find M_{ij} matrices
    mat M00,M11,M22,M01,M02,M12,M22_1;
    // M00 = M(Z0,Z0);
    // M11 = M(Z1,Z1);
    M22 = M(Z2,Z2);
    // M01 = M(Z0,Z1);
    M02 = M(Z0,Z2);
    M12 = M(Z1,Z2);

    M22_1 = inv(M22);

    // Find residuals R0,R1
    mat R0,R1;
    R0 = Z0-M02*M22_1*Z2;
    R1 = Z1-M12*M22_1*Z2;

    // R0 = Z0-sum(Z0.t()).t()/N*ones<vec>(p);
    // R1 = Z1-sum(Z1.t()).t()/N*ones<vec>(p);

    // Find matrices S_{ij}
    mat S00,S11,S01,S10,S00_1;
    S00 = S(R0,R0);
    S11 = S(R1,R1);
    S01 = S(R0,R1);
    S10 = S(R1,R0);
    S00_1 = inv(S00);

    mat solveMat;
    // mat cholS;

    solveMat = inv(S11)*S10*S00_1*S01;
    // cholS = chol(S11,"lower");

    cx_vec eigval_cx;
    cx_mat eigvec_cx;

    eig_gen(eigval_cx, eigvec_cx, solveMat);

    // C++ function returns complex vectors/matrices, so extract real parts
    vec eigval = real(eigval_cx);
    mat eigvec = real(eigvec_cx);

    // Sort by eigenvalues, descending (sort vectors first!)
    eigvec = eigvec.cols(sort_index(eigval,"descend"));
    eigval = eigval(sort_index(eigval,"descend"));

    // Normalize eigenvectors
    // for(int i=0;i<s;i++){
    //     double nf = as_scalar(sqrt(eigvec.col(i).t()*S11*eigvec.col(i)));
    //     eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    // }

    return eigval;
}

