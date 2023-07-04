#include "zoops.h"
#include "utils.h"
#include "prob_utils.h"

//' Runs Expectation Maximization ZOOPS and reestimates model parameters.
//'@name zoops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability to each sequence has a motif.
//'@return Updated PWM model.
//[[Rcpp::export]]
Rcpp::List zoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  /**
   * Posteriori
   */
  arma::vec z(m);
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  double new_w = 0.0;
  
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    
    new_alpha.fill(1e-100);
    new_w = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      /**
       * E-STEP
       */
      double marginal_prob = 0.0;
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        z[j] = w * probSeqGivenPos(seq, alpha, beta, j);
        marginal_prob += z[j];
      }
      
      // zoops
      double q = m * probSeqGivenBeta(seq, beta) * (1-w);
      marginal_prob += q;
      /**
       * M-STEP
       */
      z = z / marginal_prob;
      new_w += arma::accu(z);
      soft_update(new_alpha, seq, z.t());
    }
    double sumcol = arma::accu(new_alpha.col(0));
    alpha = new_alpha / sumcol;
    w = new_w / n;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}

//' Runs Expectation Maximization ZOOPS and reestimates model parameters.
//'@name logzoops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability to each sequence has a motif.
//'@return Updated PWM model.
//[[Rcpp::export]]
Rcpp::List logzoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  
  /**
   * Posteriori
   */
  arma::vec z(m);
  
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  double new_w = 0.0;
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    
    new_alpha.fill(1e-100);
    new_w = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      double marginal_prob = 0.0;
      
      /**
       * E-STEP
       */
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        z[j] = std::log(w) + probSeqGivenPosLog(seq, alpha, beta, j);
        marginal_prob += std::exp(z[j]);
      }
      
      // zoops
      double q = std::log(m) + probSeqGivenBetaLog(seq, beta) + std::log(1-w);
      marginal_prob += std::exp(q);
      /**
       * M-STEP
       */
      z = arma::exp(z - std::log(marginal_prob));
      new_w += arma::accu(z);
      soft_update(new_alpha, seq, z.t());
    }
    double sumcol = arma::accu(new_alpha.col(0));
    alpha = new_alpha / sumcol;
    w = new_w / n;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}
