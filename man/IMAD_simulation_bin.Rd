% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMAD_simulation_bin.R
\name{IMAD_simulation_bin}
\alias{IMAD_simulation_bin}
\title{Title}
\usage{
IMAD_simulation_bin(
  n = c(20, 20),
  p = c(0.6, 0.5),
  p_clin = c(0.6, 0.5),
  weight = 0.5,
  alpha = 0.05,
  beta = 0.8
)
}
\arguments{
\item{n}{A vector allocating initial counts to each group}

\item{p}{A vector that takes the actual simulated difference between two groups}

\item{p_clin}{A vector that takes the minimum clinically significant success probabilities of each group}

\item{weight}{The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation}

\item{alpha}{The Type 1 error rate}

\item{beta}{The Type 2 error rate}
}
\value{
The number of individuals allocated to each group, total sample size, and minimum estimated sample size
}
\description{
Title
}
\examples{
IMAD_simulation_bin()
}
