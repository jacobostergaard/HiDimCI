#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;

// Function returns the p phases (phi) in [0,2pi) for a 2p-dim input vector z of (x,y) coordinates.
arma::mat phase(arma::vec z) {
  int p = z.size()/2; // Number of bi-variate oscillators
  arma::vec out(p, fill::zeros);
  for(int i=0; i<p; i++){
    out(i) = atan2(z(2*i+1), z(2*i));        // 2-input atan function, return angles in [-pi,pi)
out(i) = fmod(2*M_PI+out(i), 2*M_PI); // Shift output to [0,2pi)
  }
  return out;
}


// Function returns the p amplitudes (gamma)for a 2p-dim input vector z of (x,y) coordinates.
arma::mat amplitude(arma::vec z) {
  int p = z.size()/2; // Number of bi-variate oscillators
  arma::vec out(p, fill::zeros);
  for(int i=0; i<p; i++){
    out(i) = sqrt(pow(z(2*i), 2)+pow(z(2*i+1), 2));
  }
  return out;
}


// Return Cartesian (x,y) coordinates, from vectors of p-dim polar coordinate input: amplitude and phase.
arma::mat polarToXY(arma::vec gam, arma::vec phi){
  int p = phi.size();
  arma::vec out(2*p, fill::zeros);
  for(int i=0; i<p; i++){
    out(2*i) = gam(i)*cos(phi(i));
    out(2*i+1) = gam(i)*sin(phi(i));
  }
  return out;
}

// Return the sign of input
int sgn(int val) {
  return (val > 0) - (val < 0);
}
