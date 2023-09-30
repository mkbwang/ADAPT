

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
  const mat& Y; // matrix of log count ratios between individual taxon and reference
  const mat& Delta; // Matrix of indicator whether
  const mat& X; // Covariate matrix

  //outputs
  RcppParallel::RMatrix<double> results; //  effect size estimate and test statistic

  CRestimate (const mat& Y, const mat& Delta, const mat& X, NumericMatrix& results):
      Y(Y), Delta(Delta), X(X), results(results){}

  void operator()(size_t begin, size_t end){

    for(size_t j=begin; j<end;j++){
      vec yvec= Y.col(j); // assume that rows are genes and cols are samples
      vec delta_vec = Delta.col(j);

      // first cary out estimation using the real data
      tobitinput  singleinput(yvec, delta_vec, X);
      tobitoutput full_estimates = estimation(&singleinput, false);
      tobitoutput null_estimates = estimation(&singleinput, true);
      results(j, 0) = full_estimates.params(1);
      double real_teststat = 2*(full_estimates.llk - null_estimates.llk);
      results(j, 1) = real_teststat;
    }
  }
};


// estimate scaling coefficient for the test statistics
struct scaleestimate: public RcppParallel::Worker{
  // inputs
  const mat& Y; // matrix of log count ratios between individual taxon and reference
  const mat& Delta; // Matrix of indicator whether
  const mat& X; // Covariate matrix
  const size_t n_sample; // number of samples
  const size_t n_boot; // number of bootstrap times

  //outputs
  RcppParallel::RMatrix<double> scale_results; //  effect size estimate and test statistic

  scaleestimate (const mat& Y, const mat& Delta, const mat& X, const size_t n_boot, NumericMatrix& scale_results):
    Y(Y), Delta(Delta), X(X), n_sample(X.n_rows), n_boot(n_boot), scale_results(scale_results){}


  void operator()(size_t begin, size_t end){

    for(size_t j=begin; j<end;j++){
      vec yvec= Y.col(j); // assume that rows are genes and cols are samples
      vec delta_vec = Delta.col(j);

      vec teststat(n_boot, fill::zeros);
      tobitinput boot_input(yvec, delta_vec, X);
      for (unsigned int k=0; k<n_boot; k++){
        Col<uword> selected_indices = randi<uvec>(n_sample, distr_param(0, n_sample-1));
        boot_input.Y = yvec.elem(selected_indices);
        boot_input.Delta = delta_vec.elem(selected_indices);
        tobitoutput boot_full_estimates = estimation(&boot_input, false);
        tobitoutput boot_null_estimates = estimation(&boot_input, true);
        teststat(k) = 2*(boot_full_estimates.llk - boot_null_estimates.llk);
      }

      scale_results(j, 1) = mean(teststat);
    }
  }

};


//' @export
// [[Rcpp::export]]
NumericMatrix cr_estim(arma::mat& Y, arma::mat& Delta, arma::mat& X) {

  size_t n_gene = Y.n_cols;
  NumericMatrix statinference(n_gene, 2);
  CRestimate cr_obj(Y, Delta, X, statinference);
  parallelFor(0, n_gene, cr_obj, 20);
  return statinference;

}


//' @export
// [[Rcpp::export]]
NumericMatrix boot_estim(arma::mat& Y, arma::mat& Delta, arma::mat& X, size_t boot_replicate=1000, size_t n_boot_gene=100){

  size_t n_gene = Y.n_cols;
  size_t n_selected_gene = (n_gene > n_boot_gene)? n_boot_gene : n_gene; // number of subset of genes
  NumericMatrix scale_chisq(n_selected_gene, 2); // store the results. First column is the gene indices. Second column is the estimates
  uvec indices = (n_gene > n_boot_gene)? randperm(n_gene, n_selected_gene) : linspace<uvec>(0, n_gene-1, n_gene);
  scale_chisq(_, 0) = as<NumericVector>(wrap(indices));
  mat subset_Y = Y.cols(indices);
  mat subset_Delta = Delta.cols(indices);

  scaleestimate chisq_estimate(subset_Y, subset_Delta, X, boot_replicate, scale_chisq);
  parallelFor(0, n_selected_gene, chisq_estimate, 20);

  return scale_chisq;

}


