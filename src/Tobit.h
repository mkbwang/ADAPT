

#include <RcppArmadillo.h>
#include <nloptrAPI.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

#ifndef TOBIT_H
#define TOBIT_H

//TODO: constructor
struct tobitinput{
  vec Y; // log count ratio
  Col<int> Delta; // censorship
  mat X; // covariate matrix
  tobitinput() = default;
  tobitinput(const vec &dependent_vec,const Col<int> &censor_vec,const mat &covar_mat):
    Y(dependent_vec), Delta(censor_vec), X(covar_mat){}
};


//TODO: constructor
struct tobitoutput{
  vec params; // estimated parameter
  double llk; // log likelihood
  nlopt_result status; // optimization status
  tobitoutput(const vec &estimate, double maxllk, nlopt_result outcome):
    params(estimate), llk(maxllk), status(outcome){}
};

// tobit loglikelihood function
double tobitllk(unsigned ndim, const double* params, double* grad, void* input);

// tobit estimation
tobitoutput estimation(void *input, bool isnull=false);

#endif //TOBIT_H
