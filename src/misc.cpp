// Misc tools

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;

// print progress in percent
void display_progress(double pct){

  // change to percent
  pct = round(pct*100);

  if(pct==100){
    cout << "\rProgress: " << pct << "% completed";
  } else if(pct >= 10){
    cout << "\rProgress:  " << pct << "% completed";
  } else{
    cout << "\rProgress:   " << pct << "% completed";
  }

}
