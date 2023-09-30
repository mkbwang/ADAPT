

#include <RcppArmadillo.h>
#include <nloptrAPI.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

#ifndef TOBIT_H
#define TOBIT_H

// input for fitting tobit models
struct tobitinput{
  vec Y; // log count ratio
  vec Delta; // censorship
  mat X; // covariate matrix
  size_t n_sample;
  double stepsize; // gradient ascent step size
  tobitinput() = default;
  tobitinput(const vec &dependent_vec, const vec &censor_vec, const mat &covar_mat):
    Y(dependent_vec), Delta(censor_vec), X(covar_mat), n_sample(covar_mat.n_rows), stepsize(0.1/covar_mat.n_rows){}
};


// output of tobit models
struct tobitoutput{
  vec params; // estimated parameter
  double llk; // log likelihood
  tobitoutput(vec &estimate, double maxllk):
    params(estimate), llk(maxllk){}
};

// tobit loglikelihood function
double tobitllk_vanilla(unsigned ndim, const double* params, double* grad, void* input); // vanilla llk, not used but a good reference
double tobitllk_firth(unsigned ndim, const double* params, double* grad, void* input); // penalized llk with Firth penalty


// tobit estimation
tobitoutput estimation(void *input, bool null=false);

#endif //TOBIT_H
