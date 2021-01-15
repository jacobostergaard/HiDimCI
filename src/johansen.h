// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef johansen_H
#define johansen_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat vecm(arma::mat X, int r, arma::mat A, arma::mat B, bool Psi, double dt);
arma::mat var(arma::mat X, bool Psi, double dt);
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, bool Psi, double dt);

// This is the end of the header guard
#endif
