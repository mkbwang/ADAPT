

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "Tobit.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]


// fit Tobit model for all the count ratios
struct CRestimate: public RcppParallel::Worker
{
  // inputs
  const mat& counts; // count matrix
  const vec& refcounts; // sum of counts in the reference set
  const mat& Delta; // Matrix of indicator whether
  const mat& X; // Covariate matrix

  //outputs
  RcppParallel::RMatrix<double> results; //  effect size estimate and test statistic

  CRestimate (const mat& count_mat, const vec& refcounts, const mat& Delta, const mat& X, NumericMatrix& results):
      counts(count_mat), refcounts(refcounts), Delta(Delta), X(X), results(results){}

  void operator()(size_t begin, size_t end){

    for(size_t j=begin; j<end;j++){
      try{
        vec count_nmt= counts.col(j);
        vec delta_vec = Delta.col(j);

        // estimation for the full model
        tobit_firth  tbmodel{count_nmt, refcounts, delta_vec, X, 1e-5, 50};
        int convergence = tbmodel.fit();
        if(convergence != 0){
          throw std::runtime_error("Full model didn't converge");
        }
        vec estimates = tbmodel.return_param();
        size_t length=estimates.n_elem;
        vec beta = estimates.head(length-1) / estimates(length-1);
        results(j, 0) = beta(1);
        double full_llk = tbmodel.return_llk();

        // estimation for the reduced model
        tbmodel.reset(true, {1});
        convergence = tbmodel.fit();
        if(convergence != 0){
          throw std::runtime_error("Reduced model didn't converge");
        }
        double reduced_llk = tbmodel.return_llk();
        double real_teststat = 2*(full_llk - reduced_llk);
        results(j, 1) = real_teststat;
      }catch(std::runtime_error& err){
        results(j, 2) = 1;
      }
    }
  }
};




// [[Rcpp::export]]
NumericMatrix cr_estim(arma::mat& count_mat, arma::vec& refcounts, arma::mat& Delta, arma::mat& X) {

  size_t n_taxa = count_mat.n_cols;
  NumericMatrix statinference(n_taxa, 3);
  CRestimate cr_obj(count_mat, refcounts, Delta, X, statinference);
  parallelFor(0, n_taxa, cr_obj, 20);
  return statinference;

}

