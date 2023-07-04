// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Q
double Q(const std::vector<std::string>& fasta, const arma::mat& alpha, const arma::mat& beta);
RcppExport SEXP _ZOOPS_Q(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(Q(fasta, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// LL
double LL(const std::vector<std::string>& fasta, const arma::mat& alpha, const arma::mat& beta);
RcppExport SEXP _ZOOPS_LL(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(LL(fasta, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// createMarkovChain
arma::mat createMarkovChain(const std::vector<std::string>& fasta, const int tau);
RcppExport SEXP _ZOOPS_createMarkovChain(SEXP fastaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(createMarkovChain(fasta, tau));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenAlpha
double probSeqGivenAlpha(const std::string& kmer, const arma::mat& alpha);
RcppExport SEXP _ZOOPS_probSeqGivenAlpha(SEXP kmerSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenAlpha(kmer, alpha));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenAlphaLog
double probSeqGivenAlphaLog(const std::string& kmer, const arma::mat& alpha);
RcppExport SEXP _ZOOPS_probSeqGivenAlphaLog(SEXP kmerSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenAlphaLog(kmer, alpha));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenBeta
double probSeqGivenBeta(const std::string& seq, const arma::mat& beta);
RcppExport SEXP _ZOOPS_probSeqGivenBeta(SEXP seqSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenBeta(seq, beta));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenBetaLog
double probSeqGivenBetaLog(const std::string& seq, const arma::mat& beta);
RcppExport SEXP _ZOOPS_probSeqGivenBetaLog(SEXP seqSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenBetaLog(seq, beta));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenPos
double probSeqGivenPos(const std::string& seq, const arma::mat& alpha, const arma::mat& beta, const int pos);
RcppExport SEXP _ZOOPS_probSeqGivenPos(SEXP seqSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenPos(seq, alpha, beta, pos));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenPosLog
double probSeqGivenPosLog(const std::string& seq, const arma::mat& alpha, const arma::mat& beta, const int pos);
RcppExport SEXP _ZOOPS_probSeqGivenPosLog(SEXP seqSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenPosLog(seq, alpha, beta, pos));
    return rcpp_result_gen;
END_RCPP
}
// computeIC
double computeIC(const arma::mat alpha, const arma::mat beta);
RcppExport SEXP _ZOOPS_computeIC(SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeIC(alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// computeICU
double computeICU(const arma::mat alpha);
RcppExport SEXP _ZOOPS_computeICU(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeICU(alpha));
    return rcpp_result_gen;
END_RCPP
}
// soft_update
void soft_update(arma::mat& alpha, const std::string& seq, const arma::rowvec& posteriori);
RcppExport SEXP _ZOOPS_soft_update(SEXP alphaSEXP, SEXP seqSEXP, SEXP posterioriSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type posteriori(posterioriSEXP);
    soft_update(alpha, seq, posteriori);
    return R_NilValue;
END_RCPP
}
// hard_update
void hard_update(arma::mat& alpha, const std::string& seq, const arma::rowvec& posteriori);
RcppExport SEXP _ZOOPS_hard_update(SEXP alphaSEXP, SEXP seqSEXP, SEXP posterioriSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type posteriori(posterioriSEXP);
    hard_update(alpha, seq, posteriori);
    return R_NilValue;
END_RCPP
}
// kmers2alpha
arma::mat kmers2alpha(const std::vector<std::string>& kmers);
RcppExport SEXP _ZOOPS_kmers2alpha(SEXP kmersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type kmers(kmersSEXP);
    rcpp_result_gen = Rcpp::wrap(kmers2alpha(kmers));
    return rcpp_result_gen;
END_RCPP
}
// alpha2kmers
std::vector<std::string> alpha2kmers(const arma::mat& alpha, const std::vector<std::string>& fasta);
RcppExport SEXP _ZOOPS_alpha2kmers(SEXP alphaSEXP, SEXP fastaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha2kmers(alpha, fasta));
    return rcpp_result_gen;
END_RCPP
}
// fasta2alpha
arma::mat fasta2alpha(const std::vector<std::string>& fasta, const int k);
RcppExport SEXP _ZOOPS_fasta2alpha(SEXP fastaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(fasta2alpha(fasta, k));
    return rcpp_result_gen;
END_RCPP
}
// corr
std::string corr(const std::string& a, const std::string b);
RcppExport SEXP _ZOOPS_corr(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(corr(a, b));
    return rcpp_result_gen;
END_RCPP
}
// corr_freq
double corr_freq(const std::string& correlation_str);
RcppExport SEXP _ZOOPS_corr_freq(SEXP correlation_strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type correlation_str(correlation_strSEXP);
    rcpp_result_gen = Rcpp::wrap(corr_freq(correlation_str));
    return rcpp_result_gen;
END_RCPP
}
// fast_corr_freq
double fast_corr_freq(const std::string& a, const std::string b);
RcppExport SEXP _ZOOPS_fast_corr_freq(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_corr_freq(a, b));
    return rcpp_result_gen;
END_RCPP
}
// computeDKL
double computeDKL(const arma::mat& alpha, const arma::mat& beta, const std::string& kmer);
RcppExport SEXP _ZOOPS_computeDKL(SEXP alphaSEXP, SEXP betaSEXP, SEXP kmerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDKL(alpha, beta, kmer));
    return rcpp_result_gen;
END_RCPP
}
// computeDKLU
double computeDKLU(const arma::mat& alpha, const std::string& kmer);
RcppExport SEXP _ZOOPS_computeDKLU(SEXP alphaSEXP, SEXP kmerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDKLU(alpha, kmer));
    return rcpp_result_gen;
END_RCPP
}
// computeSCORE
int computeSCORE(const arma::mat& alpha);
RcppExport SEXP _ZOOPS_computeSCORE(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeSCORE(alpha));
    return rcpp_result_gen;
END_RCPP
}
// zoops
Rcpp::List zoops(const std::vector<std::string>& fasta, arma::mat alpha, const arma::mat& beta, const double cutoff, int niter, double w);
RcppExport SEXP _ZOOPS_zoops(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cutoffSEXP, SEXP niterSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(zoops(fasta, alpha, beta, cutoff, niter, w));
    return rcpp_result_gen;
END_RCPP
}
// logzoops
Rcpp::List logzoops(const std::vector<std::string>& fasta, arma::mat alpha, const arma::mat& beta, const double cutoff, int niter, double w);
RcppExport SEXP _ZOOPS_logzoops(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cutoffSEXP, SEXP niterSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(logzoops(fasta, alpha, beta, cutoff, niter, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ZOOPS_Q", (DL_FUNC) &_ZOOPS_Q, 3},
    {"_ZOOPS_LL", (DL_FUNC) &_ZOOPS_LL, 3},
    {"_ZOOPS_createMarkovChain", (DL_FUNC) &_ZOOPS_createMarkovChain, 2},
    {"_ZOOPS_probSeqGivenAlpha", (DL_FUNC) &_ZOOPS_probSeqGivenAlpha, 2},
    {"_ZOOPS_probSeqGivenAlphaLog", (DL_FUNC) &_ZOOPS_probSeqGivenAlphaLog, 2},
    {"_ZOOPS_probSeqGivenBeta", (DL_FUNC) &_ZOOPS_probSeqGivenBeta, 2},
    {"_ZOOPS_probSeqGivenBetaLog", (DL_FUNC) &_ZOOPS_probSeqGivenBetaLog, 2},
    {"_ZOOPS_probSeqGivenPos", (DL_FUNC) &_ZOOPS_probSeqGivenPos, 4},
    {"_ZOOPS_probSeqGivenPosLog", (DL_FUNC) &_ZOOPS_probSeqGivenPosLog, 4},
    {"_ZOOPS_computeIC", (DL_FUNC) &_ZOOPS_computeIC, 2},
    {"_ZOOPS_computeICU", (DL_FUNC) &_ZOOPS_computeICU, 1},
    {"_ZOOPS_soft_update", (DL_FUNC) &_ZOOPS_soft_update, 3},
    {"_ZOOPS_hard_update", (DL_FUNC) &_ZOOPS_hard_update, 3},
    {"_ZOOPS_kmers2alpha", (DL_FUNC) &_ZOOPS_kmers2alpha, 1},
    {"_ZOOPS_alpha2kmers", (DL_FUNC) &_ZOOPS_alpha2kmers, 2},
    {"_ZOOPS_fasta2alpha", (DL_FUNC) &_ZOOPS_fasta2alpha, 2},
    {"_ZOOPS_corr", (DL_FUNC) &_ZOOPS_corr, 2},
    {"_ZOOPS_corr_freq", (DL_FUNC) &_ZOOPS_corr_freq, 1},
    {"_ZOOPS_fast_corr_freq", (DL_FUNC) &_ZOOPS_fast_corr_freq, 2},
    {"_ZOOPS_computeDKL", (DL_FUNC) &_ZOOPS_computeDKL, 3},
    {"_ZOOPS_computeDKLU", (DL_FUNC) &_ZOOPS_computeDKLU, 2},
    {"_ZOOPS_computeSCORE", (DL_FUNC) &_ZOOPS_computeSCORE, 1},
    {"_ZOOPS_zoops", (DL_FUNC) &_ZOOPS_zoops, 6},
    {"_ZOOPS_logzoops", (DL_FUNC) &_ZOOPS_logzoops, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ZOOPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}