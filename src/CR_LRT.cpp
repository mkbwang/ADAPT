

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


// estimate scaling coefficient for the test statistics
struct scaleestimate: public RcppParallel::Worker{
  // inputs
  const mat& counts; // count matrix
  const vec& refcounts; // sum of counts in the reference set
  const mat& Delta; // Matrix of indicator whether
  const mat& X; // Covariate matrix
  const size_t n_sample; // number of samples
  const size_t n_boot; // number of bootstrap times

  //outputs
  RcppParallel::RMatrix<double> scale_results; //  effect size estimate and test statistic

  scaleestimate (const mat& count_mat, const vec& refcounts, const mat& Delta, const mat& X, const size_t n_boot, NumericMatrix& scale_results):
    counts(count_mat), refcounts(refcounts), Delta(Delta), X(X), 
    n_sample(X.n_rows), n_boot(n_boot), scale_results(scale_results){}


  void operator()(size_t begin, size_t end){

    for(size_t j=begin; j<end;j++){
        vec count_nmt= counts.col(j); 
        vec delta_vec = Delta.col(j);

        vec teststat(n_boot, fill::zeros);
        vec failure(n_boot, fill::zeros); // numerical failure indicator
        tobit_firth tbmodel_boot{count_nmt, refcounts, delta_vec, X, 1e-5, 50};
        for (unsigned int k=0; k<n_boot; k++){
          try{
            tbmodel_boot.reorder(true);
            tbmodel_boot.reset(false, {1}); // full model with bootstrapped responses
            int convergence = tbmodel_boot.fit();
            if(convergence != 0){
              throw std::runtime_error("Full model didn't converge");
            }
            double full_llk = tbmodel_boot.return_llk();

            tbmodel_boot.reset(true, {1}); // reduced model with bootstrapped responses
            convergence = tbmodel_boot.fit();
            if (convergence != 0){
              throw std::runtime_error("Reduced model didn't converge");
            }

            double reduced_llk = tbmodel_boot.return_llk();
            teststat(k) = 2*(full_llk - reduced_llk);
          } catch(std::runtime_error& err){
            failure(k) = 1;
          }
        }

        scale_results(j, 1) = mean(teststat(find(failure == 0)));
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



// [[Rcpp::export]]
NumericMatrix boot_estim(arma::mat& count_mat, arma::vec& refcounts, arma::mat& Delta, arma::mat& X, size_t boot_replicate=2000, size_t n_boot_taxa=500){

  size_t n_taxa = count_mat.n_cols;
  size_t n_selected_taxa = (n_taxa > n_boot_taxa)? n_boot_taxa : n_taxa; // number of subset of genes
  // store the results. First column is the gene indices. Second column is the estimates.
  NumericMatrix scale_chisq(n_selected_taxa, 2);
  uvec indices = (n_taxa > n_boot_taxa)? randperm(n_taxa, n_selected_taxa) : linspace<uvec>(0, n_taxa-1, n_taxa);
  scale_chisq(_, 0) = as<NumericVector>(wrap(indices));
  mat subset_counts = count_mat.cols(indices);
  mat subset_Delta = Delta.cols(indices);

  scaleestimate chisq_estimate(subset_counts, refcounts, subset_Delta, X, boot_replicate, scale_chisq);
  parallelFor(0, n_selected_taxa, chisq_estimate, 20);

  // add indices by one for exporting to R
  scale_chisq(_, 0) = scale_chisq(_, 0) + 1;

  return scale_chisq;
}


