// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef johaTools_H
#define johaTools_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::vec f(arma::vec phi,arma::vec gam, arma::mat alpha, arma::mat beta, arma::vec omega=zeros<vec>(1));
double logLik(mat Z0, mat Z1, mat Z2, mat Psi, mat alpha, mat beta, mat Omega);
mat M(mat X,mat Y);
mat S(mat X, mat Y);
mat M_perp(mat M);
mat M_bar(mat M);
// This is the end of the header guard
#endif
