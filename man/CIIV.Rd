% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIIV_Functions.R
\name{CIIV}
\alias{CIIV}
\title{The Confidence Interval Method for Selecting Valid Instrumental Variables.}
\usage{
CIIV(
  Y,
  D,
  Z,
  X,
  alpha = 0.05,
  tuning = 0.1/log(length(Y)),
  intercept = TRUE,
  robust = TRUE
)
}
\arguments{
\item{Y}{A numeric vector of outcomes.}

\item{D}{A numeric vector of exposures/treatments.}

\item{Z}{A numeric matrix of instrumental variables, with each column
referring to one instrument.}

\item{X}{An optional numeric matrix of exogenous explanatory variables, with
each column referring to one variable.}

\item{alpha}{A numeric scalar between 0 and 1 specifying the significance
level for the confidence interval for the causal effect estimate (default =
0.05).}

\item{tuning}{A numeric scalar specifiying the threshold p-value for the
Saran/Hansen-J test (default = 0.1/log(n)).}

\item{intercept}{Logical. If intercept = TRUE, an intercept term is included
in the linear model (default = TRUE).}

\item{robust}{Logical. If robust = TRUE, the linear model is robust to
heteroskedasticity.}
}
\value{
Invalid instruments:
   Identities of the invalid instrumental variables
 selected by the algorithm.

Number of Invalid Instruments:
   The number of the selected invalid
 instrumental variables.

Coefficients:
   The matrix for the post-selection IV estimation results
 for the coefficients of the exposure/treatment variable and exogenous
 explanatory variables using the selected valid instruments as IV and
 controlling for the selected invalid instruments. The first two columns are
 2SLS estimates and their standard errors. If robust = TRUE, the second
 column is heteroskedasticity-robust stand errors. The two-step GMM estimates
 and their standard errors are also reported. If intercept = TRUE, the
 estimates for the intercept are reported in the last row.

Confidence Interval:
   The confidence interval for the 2SLS estimates for
 the coefficient of the exposure/treatment variable with significance level
 specified by alpha (default = 0.05). If robust = TRUE, the confidence
 interval for the two-step GMM estimate is also reported.

p-value of Sargan:
   The p-value for the Sargan overidentifying test for
 the selected valid instruments. If robust = TRUE, the p-value for the
 Hansen-J test is reported.
}
\description{
Implements the confidence interval (CI) downward testing procedure based on
the Sargan/Hansen-J test for overidentifying restrictions. The algorithm
selects the valid instrumental variables from a set of potential instruments
that may contain invalid ones. It also provides post-selection IV estimation
results using the selected valid instruments as IV and controlling for the
selected invalid instruments as explanatory variables.
}
\examples{
library(AER) library(sandwich)
# the MASS package is only needed to
run the working example
library(MASS)
#Generate data
n = 2000; L = 10; s = 3
pi = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
Z = matrix(rnorm(n*L),n,L)
D = 0.5 + Z \%*\% gamma + epsilon[,1]
Y = -0.5 + Z \%*\% pi + D * beta + epsilon[,2]
result = CIIV(Y,D,Z,robust=FALSE)
result
}