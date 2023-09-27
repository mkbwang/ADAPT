

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "Tobit.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]


struct CRratios: public RcppParallel::Worker
{
  // inputs
  const mat& Y; // matrix of log count ratios between individual taxon and reference
  const mat& Delta; // Matrix of indicator whether
  const mat& X; // Covariate matrix
  const size_t n_gene; // number of genes
  const size_t n_sample; // number of samples
  const size_t n_boot; // number of bootstrap times

  //outputs
  RcppParallel::RMatrix<double> results; //  effect size estimate and test statistic

  CRratios(const mat& Y, const mat& Delta, const mat& X,
           const size_t n_gene, const size_t n_sample, const size_t n_boot, NumericMatrix& results):
      Y(Y), Delta(Delta), X(X), n_gene(n_gene), n_sample(n_sample), n_boot(n_boot), results(results){}

  void operator()(size_t begin, size_t end){

    for(size_t j=begin; j<end;j++){
      vec yvec= Y.col(j); // assume that rows are genes and cols are samples
      vec delta_vec = Delta.col(j);

      // first cary out estimation using the real data
      tobitinput  singleinput(yvec, delta_vec, X);
      tobitoutput full_estimates = estimation(&singleinput, false);
      tobitoutput null_estimates = estimation(&singleinput, true);
      results(j, 0) = full_estimates.params(1);
      results(j, 1) = full_estimates.status;
      double real_teststat = 2*(full_estimates.llk - null_estimates.llk);
      results(j, 2) = real_teststat;
      // bootstrap
      //double avg_boot_teststat=0;
      vec teststat(n_boot, fill::zeros);
      tobitinput boot_input(yvec, delta_vec, X);
      for (unsigned int k=0; k<n_boot; k++){
        Col<uword> selected_indices = randi<uvec>(n_sample, distr_param(0, n_sample-1));
        boot_input.Y = singleinput.Y.elem(selected_indices);
        boot_input.Delta = singleinput.Delta.elem(selected_indices);
        tobitoutput boot_full_estimates = estimation(&boot_input, false);
        tobitoutput boot_null_estimates = estimation(&boot_input, true);
        teststat(k) = 2*(boot_full_estimates.llk - boot_null_estimates.llk);
      }
      results(j, 3)= mean(teststat); // Bartlett corrected test statistic
    }
  }
};

//' @export
// [[Rcpp::export]]
NumericMatrix cr_lrt(arma::mat& Y, arma::mat& Delta, arma::mat& X, size_t n_gene, size_t n_sample, size_t n_boot=1000) {

  NumericMatrix statinference(n_gene, 4);
  CRratios cr_obj(Y, Delta, X, n_gene, n_sample, n_boot, statinference);
  parallelFor(0, n_gene, cr_obj, 20);
  return statinference;

}


