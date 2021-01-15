// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef TOOLS_H
#define TOOLS_H


// This is the content of the .h file, which is where the declarations go
arma::mat phase(arma::vec z);
arma::mat amplitude(arma::vec z);
arma::mat polarToXY(arma::vec gam, arma::vec phi);

// This is the end of the header guard
#endif
