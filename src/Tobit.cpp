#include "Tobit.h"

// First define the functions in the vanilla model
tobit_vanilla::tobit_vanilla(const vec& Nmt_input, const vec& Denom_input, const vec& delta_input, const mat&X_input,
                             double tolerance, size_t maxiter):
  model(zeros(size(delta_input)), delta_input, X_input, tolerance, maxiter),
  Nmt(Nmt_input), Denom(Denom_input), Delta_orig(delta_input)
{
  Y_orig = -log(Denom);
  uvec existence_indices = find(Delta_orig > 0);
  Y_orig.elem(existence_indices) += log(Nmt.elem(existence_indices));
  reorder(false);
  reset();
}

void tobit_vanilla::reset(bool reduced, uvec null_indices){

  params = vec(P+1, fill::zeros);
  params(0) = mean(Y)/stddev(Y);
  params(P) = 1/stddev(Y);
  deriv_z = vec(N, fill::zeros);
  deriv_2z = vec(N, fill::zeros);
  uvec fixed = uvec(P+1, fill::zeros);
  isreduced = reduced;
  if (isreduced){
    fixed.elem(null_indices) = uvec(null_indices.n_elem, fill::ones);
  }
  subindices = find(fixed == 0);
  step_working_params = vec(subindices.n_elem, fill::zeros);
  score = vec(P+1, fill::zeros);
  working_score = score(subindices);
  hessian = mat(P+1, P+1, fill::zeros);
  working_hessian = hessian(subindices, subindices);
  update_utils();
  update_llk();
  iter_counter = 0;
  convergence_code = SUCCESS;
} // reset the parameters

void tobit_vanilla::reorder(bool bootstrap){

  if(bootstrap){
    uvec boot_indices = randi<uvec>(N, distr_param(0, N-1));
    Delta = Delta_orig(boot_indices);
    if (all(Delta == 0)){ // do not allow all data to be censored
      throw std::runtime_error("Absence in all samples");
    }
    uvec existence_indices = find(Delta > 0); // samples that are chosen to have observed counts
    vec boot_Y = Y_orig.elem(boot_indices.elem(existence_indices)); // bootstrapped ratios for nonzero counts
    vec boot_denom = Denom.elem(existence_indices); // the denominators of the bootstrapped samples with nonzero counts
    vec boot_counts = round(boot_denom % exp(boot_Y)); // the nonzero counts of the bootstrapped samples
    if (all(boot_counts == 0)){ // do not allow all data to be censored
      throw std::runtime_error("Absence in all samples");
    }
    // filter out the extra sampling zeros
    existence_indices = existence_indices.elem(find(boot_counts > 0));
    boot_counts = boot_counts.elem(find(boot_counts > 0));
    // modify the existence vector
    Delta = zeros(N);
    Delta.elem(existence_indices) = ones(existence_indices.n_elem);
    Y = -log(Denom);
    Y.elem(existence_indices) += log(boot_counts);
  } else{
    Y = Y_orig;
    Delta = Delta_orig;
  }

}//bootstrap the responses

// update the vectors related to z
void tobit_vanilla::update_utils(){
  Z = params(P)*Y -  X*params.head(P);
  cumnorm_z = normcdf(Z);
  exp_2z2 = exp(-0.5*square(Z));
}


double tobit_vanilla::tobit_vanilla_llk(){
  vec llks = Delta % (-0.5*square(Z) - 0.5*log(2*datum::pi) + log(params(P))) + \
    (1 - Delta) % log(cumnorm_z);
  return accu(llks);
}

void tobit_vanilla::update_llk(){
  llk = tobit_vanilla_llk();
}

vec tobit_vanilla::tobit_vanilla_score() {
  deriv_z = -Delta % Z + 1/sqrt(2*datum::pi) * (1-Delta) / cumnorm_z % exp_2z2;
  vec new_score(P+1, fill::zeros);
  new_score.head(P) = - X.t() * deriv_z;
  new_score(P) = accu(Delta)/params(P) + accu(deriv_z % Y) ;
  return new_score;
}


void tobit_vanilla::update_score(){
  score = tobit_vanilla_score();
  if(isreduced){ // reduced model
    working_score = score(subindices);
  } else{
    working_score = score;
  }
}

mat tobit_vanilla::tobit_vanilla_hessian(){

  mat new_hessian(P+1, P+1, fill::zeros);
  deriv_2z = -Delta - (1 - Delta) / (2*datum::pi) % square(exp_2z2) / square(cumnorm_z) - \
    (1 - Delta) / sqrt(2*datum::pi) % Z % exp_2z2 / cumnorm_z;
  new_hessian(span(0, P-1), span(0, P-1)) =  X.t() * diagmat(deriv_2z) * X;
  new_hessian(P, P) = - 1/pow(params(P), 2) * accu( Delta) + accu(deriv_2z % square(Y));
  new_hessian(span(0, P-1), P) = - X.t() * (deriv_2z % Y);
  new_hessian(P, span(0, P-1)) = new_hessian(span(0, P-1), P).as_row();

  return new_hessian;
}


int tobit_vanilla::update_hessian(){
  hessian = tobit_vanilla_hessian();
  if(isreduced){
    working_hessian = hessian(subindices, subindices);
  } else{
    working_hessian = hessian;
  }

  mat information = -working_hessian;
  if (information.is_sympd()){
    return SUCCESS;
  } else{
    return FAIL;
  }
}


void tobit_vanilla::update_param(){
  step_working_params = solve(-working_hessian, working_score, arma::solve_opts::likely_sympd);

  // prevent the inverse scale parameter to become smaller than zero
  double inv_scale_step = step_working_params(step_working_params.n_elem-1);
  if (params(P) + inv_scale_step < 0){
    step_working_params = step_working_params * fabs(params(P)/inv_scale_step) / 2;
  }

  if (isreduced) {// reduced model
    params(subindices) = params(subindices) + step_working_params;
  } else{// full model
    params= params + step_working_params;
  }
}

int tobit_vanilla::fit(){

  while(iter_counter < maxiter){
    update_score();
    int hessian_check = update_hessian();
    if (hessian_check > 0){
      convergence_code = FAIL;
      break;
    }
    update_param();
    iter_counter++;
    update_utils();
    double old_llk = llk;
    update_llk();
    if(fabs((llk - old_llk)/old_llk) < tolerance){
      break;
    }
  }

  if(iter_counter == maxiter){
    convergence_code = STOPEARLY;
  }
  return convergence_code;
}

int tobit_vanilla::return_iterations() {
  return iter_counter;
}

double tobit_vanilla::return_prevalences() {
  return accu(Delta)/N;
}

// extra functions/overrides for tobit_firth
tobit_firth::tobit_firth(const vec& Nmt_input, const vec& Denom_input, const arma::vec &delta_input, const arma::mat &X_input,
                         double tolerance, size_t maxiter):
  tobit_vanilla(Nmt_input, Denom_input, delta_input, X_input, tolerance, maxiter){
  deriv_3z = vec(N, fill::zeros);
  step_working_score = vec(subindices.n_elem, fill::zeros);
  hessian = tobit_vanilla_hessian();
}

double tobit_firth::tobit_firth_llk() {
  hessian = tobit_vanilla_hessian();
  return 0.5*log_det_sympd(-hessian);
}

void tobit_firth::update_llk() {
  llk = tobit_vanilla_llk() + tobit_firth_llk();
}

vec tobit_firth::tobit_firth_score() {
  inv_information = inv_sympd(-hessian);
  deriv_3z = (1-Delta)/sqrt(2*datum::pi) / cumnorm_z % exp_2z2 %
    (1 / datum::pi / square(cumnorm_z) % square(exp_2z2) + 3 / sqrt(2*datum::pi) / cumnorm_z % Z % exp_2z2 - 1 + square(Z));

  vec firth_score(P+1, fill::zeros);
  mat information_deriv(size(inv_information), fill::zeros);
  for (size_t k=0; k<P; k++){ // scores for all the effect sizes
    vec xvec = X.col(k);
    information_deriv(span(0, P-1), span(0, P-1)) = X.t() * diagmat(deriv_3z % xvec) * X;
    information_deriv(P, P) = accu(deriv_3z % xvec % square(Y));
    information_deriv(span(0, P-1), P) =  - X.t() * (deriv_3z % xvec % Y);
    information_deriv(P, span(0, P-1)) = information_deriv(span(0, P-1), P).as_row();
    firth_score(k) = 0.5 * trace(inv_information* information_deriv); // gradient of the kth effect size
  }
  // score for the inverse scale
  information_deriv(span(0, P-1), span(0, P-1)) = -X.t() * diagmat(deriv_3z % Y) * X;
  information_deriv(P, P) = - accu(deriv_3z % Y % square(Y) ) - 2*accu(Delta)/pow(params(P), 3);
  information_deriv(span(0, P-1), P) = X.t() * (deriv_3z % square(Y));
  information_deriv(P, span(0, P-1)) = information_deriv(span(0, P-1), P).as_row();
  firth_score(P) = 0.5 * trace(inv_information* information_deriv);

  return firth_score;
}

void tobit_firth::update_score() {
  score = tobit_vanilla_score() + tobit_firth_score();
  if(isreduced){ // reduced model
    step_working_score = score(subindices) - working_score;
    working_score = score(subindices);
  } else{
    step_working_score = score - working_score;
    working_score = score;
  }
}


int tobit_firth::update_hessian(){

  if (iter_counter == 0){ // initial starting point
    // use the hessian from the vanilla tobit model as approximation
    if (isreduced){
      working_hessian = hessian(subindices, subindices);
    } else{
      working_hessian = hessian;
    }
  } else{ // BFGS update
    vec approx_score_step = working_hessian * step_working_params;
    working_hessian += step_working_score * step_working_score.t() / dot(step_working_score, step_working_params) - \
      approx_score_step * approx_score_step.t() / dot(approx_score_step, step_working_params);
  }

  mat working_information = -working_hessian;
  if (working_information.is_sympd()){
    return SUCCESS;
  } else{
    return FAIL;
  }

}

