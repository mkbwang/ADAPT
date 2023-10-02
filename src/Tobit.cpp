#include "Tobit.h"

// vanilla Tobit likelihood
double tobitllk_vanilla(unsigned ndim, const double* params, double* grad, void* input){

  auto inputdata= (const tobitinput *) input;
  vec rho(params, ndim-1); // beta/sigma
  double omega = exp(params[ndim-1]); // 1/sigma=omega=exp(phi), precision
  vec z = omega*inputdata->Y - inputdata->X * rho;
  vec cumnorm_z = normcdf(z);
  vec exp_2z2 = exp(-0.5*square(z));

  // calculate vanilla loglikelihood
  vec llk = inputdata->Delta % (-0.5*square(z) -0.5*log(2*datum::pi) + log(omega))  + \
    (1-inputdata->Delta) % log(cumnorm_z); // log likelihood contributed from each observation

  // calculate gradient from vanilla loglikelihood
  vec deriv_rho(grad, ndim-1, false, true);
  vec deriv_z = -inputdata->Delta % z + \
    1/sqrt(2*datum::pi) * (1-inputdata->Delta) / cumnorm_z % exp_2z2;


  deriv_rho =  - inputdata->X.t() * deriv_z * inputdata->stepsize; // gradient for the elements of rho
  grad[ndim-1] = accu(inputdata->Delta) * inputdata->stepsize+ accu(deriv_z % inputdata->Y)*omega * inputdata->stepsize; // gradient for phi

  // return log likelihood
  return accu(llk) ;
}


double tobitllk_firth(unsigned ndim, const double* params, double* grad, void* input){

  auto inputdata= (const tobitinput *) input;
  vec rho(params, ndim-1); // beta/sigma
  double omega = exp(params[ndim-1]); // 1/sigma=omega=exp(phi), sqrt(precision)
  vec z = omega*inputdata->Y - inputdata->X * rho;
  vec cumnorm_z = normcdf(z);
  vec exp_2z2 = exp(-0.5*square(z));

  // calculate vanilla loglikelihood
  vec llk = inputdata->Delta % (-0.5*square(z) -0.5*log(2*datum::pi) + log(omega))  + \
    (1-inputdata->Delta) % log(cumnorm_z); // log likelihood contributed from each observation

  // calculate gradient from vanilla loglikelihood
  vec deriv_rho(grad, ndim-1, false, true);
  vec deriv_z = -inputdata->Delta % z + \
    1/sqrt(2*datum::pi) * (1-inputdata->Delta) / cumnorm_z % exp_2z2;

  deriv_rho =  - inputdata->X.t() * deriv_z * inputdata->stepsize; // gradient for the elements of rho
  grad[ndim-1] = accu(inputdata->Delta)* inputdata->stepsize+ accu(deriv_z % inputdata->Y)*omega* inputdata->stepsize; // gradient for phi

  // add gradients from the Firth penalty
  mat information(ndim, ndim, fill::randu); // information matrix
  // negative 2nd derivative of log likelihood over z
  vec neg_deriv_z2 =  inputdata->Delta + (1-inputdata->Delta) / (2*datum::pi) % square(exp_2z2) / square(cumnorm_z) + \
    (1 - inputdata->Delta) / sqrt(2*datum::pi) % z % exp_2z2 / cumnorm_z;

  information(span(0, ndim-2), span(0, ndim-2)) = \
    inputdata->X.t() * diagmat(neg_deriv_z2) * inputdata->X; // 2nd degree derivative over rho
  information(ndim-1, ndim-1) =  accu(neg_deriv_z2 % square(inputdata->Y)) + accu(inputdata->Delta) / pow(omega, 2); // 2nd degree over phi
  information(span(0, ndim-2), ndim-1) = - inputdata->X.t() * (neg_deriv_z2 % inputdata->Y);
  information(ndim-1, span(0, ndim-2)) = information(span(0, ndim-2), ndim-1).as_row();

  if (!information.is_sympd()){ // keep an eye on possible numerical issues
    throw std::overflow_error("Ill formed information matrix");
  }

  mat inv_information = inv_sympd(information);

  // negative 3rd derivative of log likelihood over z
  vec neg_deriv_z3 = (1-inputdata->Delta)/sqrt(2*datum::pi) / cumnorm_z % exp_2z2 %
    (-1 / datum::pi / square(cumnorm_z) % square(exp_2z2) - 3 / sqrt(2*datum::pi) / cumnorm_z % z % exp_2z2 + 1 - square(z));

  mat information_deriv(ndim, ndim, fill::randu); // derivative of information matrices
  for (size_t k=0; k<ndim-1; k++){
    vec xvec = inputdata->X.col(k);
    information_deriv(span(0, ndim-2), span(0, ndim-2)) = - inputdata->X.t() * diagmat(neg_deriv_z3 % xvec) * inputdata->X;
    information_deriv(ndim-1, ndim-1) = - accu(neg_deriv_z3 % xvec % square(inputdata->Y));
    information_deriv(span(0, ndim-2), ndim-1) = inputdata->X.t() * (neg_deriv_z3 % xvec % inputdata->Y);
    information_deriv(ndim-1, span(0, ndim-2)) = information_deriv(span(0, ndim-2), ndim-1).as_row();
    grad[k] += 0.5 * trace(inv_information* information_deriv)* inputdata->stepsize; // gradient of the kth covariate coefficient
  }
  information_deriv(span(0, ndim-2), span(0, ndim-2)) = inputdata->X.t() * diagmat(neg_deriv_z3 % inputdata->Y) * inputdata->X;
  information_deriv(ndim-1, ndim-1) = accu(neg_deriv_z3 % inputdata->Y % square(inputdata->Y) ) - 2*accu(inputdata->Delta)/pow(omega, 3);
  information_deriv(span(0, ndim-2), ndim-1) = - inputdata->X.t() * (neg_deriv_z3 % square(inputdata->Y));
  information_deriv(ndim-1, span(0, ndim-2)) = information_deriv(span(0, ndim-2), ndim-1).as_row();
  grad[ndim-1] += 0.5 * trace(inv_information* information_deriv)*omega * inputdata->stepsize; // gradient of the log precision

  double firth_penalty = 0.5 * log_det_sympd(information);

  // return log likelihood
  return accu(llk) + firth_penalty;
}


tobitoutput estimation(void *input, bool null){

  auto inputdata= (const tobitinput *) input;
  const unsigned int n_dim = inputdata->X.n_cols+1;

  // optimization object
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, n_dim);
  // nlopt_set_lower_bound(opt, n_dim-1, 0); // the inverse scale parameter is positive
  std::vector<double> lower_bounds(n_dim, -HUGE_VAL);
  std::vector<double> upper_bounds(n_dim, +HUGE_VAL);
  if (null) { // fitting null model
    lower_bounds[1] = 0;
    upper_bounds[1] = 0;
  }
  nlopt_set_lower_bounds(opt, &lower_bounds[0]);
  nlopt_set_upper_bounds(opt, &upper_bounds[0]);
  nlopt_set_maxeval(opt, 40);
  nlopt_set_max_objective(opt, tobitllk_firth, input);
  nlopt_set_ftol_rel(opt, 1e-4);
  nlopt_set_ftol_abs(opt, 1e-5);
  nlopt_set_vector_storage(opt, 5); // specific for L-BFGS, number of past gradients
  // set up the parameter vector to estimate
  vec param_estimate(n_dim, fill::zeros);
  param_estimate(0) = mean(inputdata->Y)/stddev(inputdata->Y); // initialize the intercept
  param_estimate(n_dim-1) = -log(stddev(inputdata->Y)); // initialize the inverse of standard deviation
  double *param_pt = param_estimate.memptr();
  double llk; //loglikelihood

  nlopt_optimize(opt, param_pt, &llk);
  nlopt_destroy(opt);

  // transform the estimates of rho and omega into beta and sigma
  param_estimate.subvec(0, n_dim-2) = param_estimate.subvec(0, n_dim-2)/ exp(param_estimate(n_dim-1));
  param_estimate(n_dim-1) = 1/exp(param_estimate(n_dim-1));

  tobitoutput output(param_estimate, llk);

  return output;
}

