#define R_NO_REMAP

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

//#include <math.h>           // math routines
#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>


#include <R.h>              // R header
#include "RcppArmadillo.h"  // RcppArmadillo library header


#define Malloc(type,n) (type *)malloc((n)*sizeof(type))







using namespace std;
using namespace Rcpp;



// [[Rcpp::export]]
arma::vec SVDS(arma::sp_mat  &X) {
  arma::mat U;
  arma::vec S;
  arma::mat V;

  arma::svds(U,S,V,X, 1,0.000001);
  return S;
}
