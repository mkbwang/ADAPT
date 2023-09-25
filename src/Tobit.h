

#include <RcppArmadillo.h>
#include <nloptrAPI.h>



#ifndef TOBIT_H
#define TOBIT_H

//TODO: constructor
struct tobitinput{
  arma::vec Y; // log count ratio
  arma::Col<int> Delta; // censorship
  arma::mat X; // covariate arma::matrix
  tobitinput() = default;
  tobitinput(const arma::vec &dependent_vec,const arma::Col<int> &censor_vec,const arma::mat &covar_mat):
    Y(dependent_vec), Delta(censor_vec), X(covar_mat){}
};


//TODO: constructor
struct tobitoutput{
  arma::vec params; // estiarma::mated parameter
  double llk; // log likelihood
  nlopt_result status; // optimization status
  tobitoutput(const arma::vec &estimate, double maxllk, nlopt_result outcome):
    params(estimate), llk(maxllk), status(outcome){}
};

// tobit loglikelihood function
double tobitllk(unsigned ndim, const double* params, double* grad, void* input);

// tobit estiarma::mation
tobitoutput estimation(void *input, bool null=false);

#endif //TOBIT_H
