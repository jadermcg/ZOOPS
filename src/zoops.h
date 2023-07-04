#pragma once
#include <RcppArmadillo.h>
#include <strings.h>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]

Rcpp::List logzoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
Rcpp::List zoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
