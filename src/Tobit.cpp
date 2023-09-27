#include "Tobit.h"


double tobitllk(unsigned ndim, const double* params, double* grad, void* input){

  auto inputdata= (const tobitinput *) input;
  vec rho(params, ndim-1); // beta/sigma
  double omega = exp(params[ndim-1]); // 1/sigma

  vec z = omega*inputdata->Y - inputdata->X * rho;
  vec cumnorm_z = normcdf(z);

  // calculate loglikelihood
  vec llk = inputdata->Delta % (-0.5*square(z) -0.5*log(2*datum::pi) + log(omega))  + \
    (1-inputdata->Delta) % log(cumnorm_z); // log likelihood contributed from each observation

  // calculate gradient
  vec deriv_rho(grad, ndim-1, false, true);
  vec deriv_z = -inputdata->Delta % z + \
    1/sqrt(2*datum::pi) * (1-inputdata->Delta) / cumnorm_z % exp(-0.5*square(z));

  deriv_rho =  - inputdata->X.t() * deriv_z; // gradient for the elements of rho
  grad[ndim-1] = (accu(inputdata->Delta / omega) + accu(deriv_z % inputdata->Y)) * omega;

  // return log likelihood
  return accu(llk);
}



tobitoutput estimation(void *input, bool isnull){

  auto inputdata= (const tobitinput *) input;
  const unsigned int n_dim = inputdata->X.n_cols+1;

  // optimization object
  nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, n_dim);
  std::vector<double> lower_bounds(n_dim, -HUGE_VAL);
  std::vector<double> upper_bounds(n_dim, +HUGE_VAL);
  if (isnull) { // fitting isnull model
    lower_bounds[1] = 0;
    upper_bounds[1] = 0;
  }
  nlopt_set_lower_bounds(opt, &lower_bounds[0]);
  nlopt_set_upper_bounds(opt, &upper_bounds[0]);

  nlopt_set_max_objective(opt, tobitllk, input);
  nlopt_set_ftol_rel(opt, 5e-4);
  nlopt_set_ftol_abs(opt, 5e-3);
  nlopt_set_maxeval(opt, 40);
  // set up the parameter vector to estimate
  vec param_estimate(n_dim, fill::zeros);
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

