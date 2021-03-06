% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmlate.R
\name{kmlate}
\alias{kmlate}
\title{Kaplan-Meier Local Average Treatment Effect}
\usage{
kmlate(out, delta, treat, z, xpscore, b = 1000, ci = c(0.9, 0.95, 0.99),
  trunc = NULL, cores = 1)
}
\arguments{
\item{out}{vector containing the outcome of interest}

\item{delta}{vector containing the censoring indicator (1 if observed, 0 if censored)}

\item{treat}{vector containing the treatment indicator (1 if treated, 0 if control)}

\item{z}{vector containing the binary instrument}

\item{xpscore}{matrix (or data frame) containing the covariates (and their
transformations) to be included in the instrument propensity score estimation.
Instrument Propensity score estimation is based on Logit.}

\item{b}{The number of bootstrap replicates to be performed. Default is 1,000.}

\item{ci}{A scalar or vector with values in (0,1) containing the confidence level(s)
of the required interval(s). Default is a vector with
0,90, 0.95 and 0.99}

\item{trunc}{scalar that defined the truncation parameter. Default is NULL, which does not perform any kind of
truncation in the computation of the ATE. When trunc is different than NULL, all outcomes which values greater
than trunc are truncated.}

\item{cores}{number of processesors to be used during the bootstrap (default is 1).
If cores>1, the bootstrap is conducted using snow}
}
\value{
a list containing the local average treatment effect estimate, late,
       and the bootstrapped \emph{ci} confidence
       confidence interval, late.lb (lower bound), and late.ub (upper bound).
}
\description{
\emph{kmlate} computes the Local Average Treatment Effect for possibly right-censored outcomes.
The estimator relies on the availability of an Instrumental variable Z, and on a monotonicity assumption.
To implement the estimator, we make use of an instrumental propensity score approach.
For details of the estimation procedure, see Sant'Anna (2016a), 'Program Evaluation with
Right-Censored Data'.
}

