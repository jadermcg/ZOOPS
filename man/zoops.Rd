% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{zoops}
\alias{zoops}
\title{Runs Expectation Maximization ZOOPS and reestimates model parameters.}
\usage{
zoops(fasta, alpha, beta, cutoff, niter)
}
\arguments{
\item{fasta}{Dataset of sequences.}

\item{alpha}{PWM model.}

\item{beta}{0-order Markov Chain.}

\item{cutoff}{Cutoff for EM convergence.}

\item{niter}{Maximum number of iterations.}
}
\value{
Updated PWM model.
}
\description{
Runs Expectation Maximization ZOOPS and reestimates model parameters.
}
