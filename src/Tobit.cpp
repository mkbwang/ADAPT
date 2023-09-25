#include "Tobit.h"


double tobitllk(unsigned ndim, const double* params, double* grad, void* input){

  auto inputdata= (const tobitinput *) input;
  arma::vec rho(params, ndim-1); // beta/sigma
  double omega = exp(params[ndim-1]); // 1/sigma

  arma::vec z = omega*inputdata->Y - inputdata->X * rho;
  arma::vec cumnorm_z = normcdf(z);

  // calculate loglikelihood
  arma::vec llk = inputdata->Delta % (-0.5*square(z) -0.5*log(2*arma::datum::pi) + log(omega))  + \
    (1-inputdata->Delta) % log(cumnorm_z); // log likelihood contributed from each observation

  // calculate gradient
  arma::vec deriv_rho(grad, ndim-1, false, true);
  arma::vec deriv_z = -inputdata->Delta % z + \
    1/sqrt(2*arma::datum::pi) * (1-inputdata->Delta) / cumnorm_z % exp(-0.5*square(z));

  deriv_rho =  - inputdata->X.t() * deriv_z; // gradient for the elements of rho
  grad[ndim-1] = (accu(inputdata->Delta / omega) + accu(deriv_z % inputdata->Y)) * omega;

  // return log likelihood
  return accu(llk);
}



tobitoutput estimation(void *input, bool null){

  auto inputdata= (const tobitinput *) input;
  const unsigned int n_dim = inputdata->X.n_cols+1;

  // optimization object
  nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, n_dim);
  // nlopt_set_lower_bound(opt, n_dim-1, 0); // the inverse scale parameter is positive
  if (null) { // fitting null model
    nlopt_set_lower_bound(opt, 1, 0);
    nlopt_set_upper_bound(opt, 1, 0);
  }
  nlopt_set_max_objective(opt, tobitllk, input);
  nlopt_set_xtol_rel(opt, 1e-4);
  // set up the parameter arma::vector to estimate
  arma::vec param_estimate(n_dim, arma::fill::zeros);
  param_estimate(0) = mean(inputdata->Y)/stddev(inputdata->Y); // initialize the intercept
  param_estimate(n_dim-1) = -log(stddev(inputdata->Y)); // initialize the inverse of standard deviation
  double *param_pt = param_estimate.memptr();
  double llk; //loglikelihood

  nlopt_result status = nlopt_optimize(opt, param_pt, &llk);
  nlopt_destroy(opt);

  // transform the estimates of rho and omega into beta and sigma
  param_estimate.subvec(0, n_dim-2) = param_estimate.subvec(0, n_dim-2)/exp(param_estimate(n_dim-1));
  param_estimate(n_dim-1) = 1/exp(param_estimate(n_dim-1));

  tobitoutput output(param_estimate, llk, status);

  return output;
}

