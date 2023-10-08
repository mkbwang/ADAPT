

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

#ifndef TOBIT_H
#define TOBIT_H

#define SUCCESS 0;
#define FAIL 1;
#define STOPEARLY 2;


class model{
protected:
  vec Y; // response variable, dimension N
  vec Delta; // indicator of nonzero counts, dimension N
  mat X; // covariate matrix (including intercept), dimension N*P
  double tolerance; // relative tolerance of llk to stop the NR/BFGS iterations
  size_t maxiter; // max iterations
  size_t N; // number of individuals N
  size_t P; // number of effect sizes to estimate P
  vec params; // parameters to estimate, including rho(effect sizes) and phi(-log scale), dimension P+1
  double llk; // log likelihood, scalar

public:
  model(const vec& Y_input, const vec& delta_input, const mat&X_input,
        double tolerance = 1e-4, size_t maxiter=50):
  Y(Y_input), Delta(delta_input), X(X_input),tolerance(tolerance), maxiter(maxiter),
  N(X_input.n_rows), P(X_input.n_cols),  params(vec(X_input.n_cols+1, fill::zeros)),
  llk(0){};

  // update some utility variables for calculation of likelihood, score and hessians
  virtual void update_utils() = 0;
  //calculate log likelihood
  virtual void update_llk() = 0;
  //update derivatives
  virtual void update_score() = 0;
  //update hessian
  virtual int update_hessian() = 0; // check if negative hessian (information) is positive definite
  //update parameter
  virtual void update_param() = 0;
  //master function for fitting the models
  virtual int fit() = 0;
  //return estimated values
  vec return_param(){
    return params;
  };
  double return_llk(){
    return llk;
  };
  virtual ~model() = default; //destructor

};

class tobit_vanilla: public model{
protected:
  bool isreduced; // whether you are fitting a reduced model
  vec Y_orig; // original response vector (in contrast with potentially bootstrapped Y)
  vec Delta_orig; // original incidence vector (in contrast with potentially boostrapped Delta)
  vec Z; // exp(phi)*Y - X^t rho, dimension N
  vec cumnorm_z; // cumulative distribution of standard normal for each individual z, dimension N
  vec exp_2z2; // exp(-0.5*z^2), dimension N
  vec deriv_z; // first derivative of llk over each individual z, dimension N
  vec deriv_2z; // second derivative of llk over each individual z, dimension N
  uvec subindices; // indices of parameters to be estimated (for null models)
  vec step_working_params; // most recent step of the params
  vec score; // derivative of llk over all the parameters parameter (score equation), dimension P+1
  vec working_score; // shorter than score when fitting a reduced model
  mat hessian; // second derivative (hessian) of llk over each parameter, dimension (P+1)*(P+1)
  mat working_hessian; // smaller than hessian when fitting a reduced model
  size_t iter_counter;
  int convergence_code; // an integer representing the convergence codes
public:
  tobit_vanilla(const vec& Y_input, const vec& delta_input, const mat&X_input,
                double tolerance = 1e-4, size_t maxiter=50);
  void reset(bool reduced=false, uvec null_indices = {1});
  void reorder(bool bootstrap=false);
  void update_utils() override;
  void update_llk() override;
  double tobit_vanilla_llk() ;
  void update_score() override;
  vec tobit_vanilla_score() ;
  int update_hessian() override;
  mat tobit_vanilla_hessian();
  void update_param() override;
  int fit()override;
  int return_iterations();
  double return_prevalences();
};

class tobit_firth: public tobit_vanilla{
protected:
  vec deriv_3z; // third derivative of llk over each individual z, dimension N
  vec step_working_score; // most recent step of the score vector
  mat inv_information; // inverse of information matrix
public:
  tobit_firth(const vec& Y_input, const vec& delta_input, const mat&X_input,
              double tolerance = 1e-4, size_t maxiter=50);

  void update_llk() override;
  double tobit_firth_llk();
  void update_score() override;
  vec tobit_firth_score();
  int update_hessian() override;

};

#endif //TOBIT_H
