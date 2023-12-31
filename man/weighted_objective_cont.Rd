% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted_objective_cont.R
\name{weighted_objective_cont}
\alias{weighted_objective_cont}
\title{Weighted Objective Functions for Continuous Outcomes}
\usage{
weighted_objective_cont(
  sda = 1,
  sdb = 1,
  fa = 0.5,
  fb = 0.5,
  cme = 3.765,
  weight = 0.5,
  n_a = 10,
  n_b = 10,
  n = 20,
  zalpha = qnorm(0.05/2),
  zbeta = qnorm(0.2),
  nu = 1
)
}
\arguments{
\item{sda}{The standard deviation of the first group}

\item{sdb}{The standard deviation of the second group}

\item{fa}{Location-invariant value of the first group}

\item{fb}{Location invariant value of the second group}

\item{cme}{The set constraint to ensure the difference is at least as great as the clinically meaningful difference}

\item{weight}{The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation}

\item{n_a}{The number of individuals currently in the first group}

\item{n_b}{The number of individuals currently in the second group}

\item{n}{The total number of individuals currently in either group}

\item{zalpha}{The standardized set type 1 error rate}

\item{zbeta}{The standardized set type 2 error rate}

\item{nu}{The set clinically meaningful difference}
}
\value{
Returns optimal allocation, optimal sample size ratio, and allocation of the next patient for continuous normal distributions
}
\description{
Given the standard deviations, positive functions of theta, and sample size of two treatment arms along with the weighting and minimally clinically meaningful difference,
the optimal sample size, allocation ratio, and sample size ratio are found using this function. From these the minimum sample size of each group, the biased coin probability,
and the allocation of the next participant from tht probability are given.
}
\examples{
weighted_objective_cont()
}
