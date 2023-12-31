% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMAD_next_val_cont.R
\name{IMAD_next_val_cont}
\alias{IMAD_next_val_cont}
\title{Allocation using IMAD for Continuous Variables}
\usage{
IMAD_next_val_cont(
  group_a,
  group_b,
  trtlvls = c("A", "B"),
  nu,
  pinought = 0.5,
  weight = 0.5,
  alpha = 0.05,
  beta = 0.2
)
}
\arguments{
\item{group_a}{A vector containing the values of the first group}

\item{group_b}{A vector containing the values of the second group}

\item{trtlvls}{A vector containing the labels of each respective group}

\item{nu}{The preset clinically meaningful difference to be looked for}

\item{pinought}{The target allocation ratio for the better performing group}

\item{weight}{The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation}

\item{alpha}{The Type 1 error rate}

\item{beta}{The Type 2 error rate}
}
\value{
The grouping of the next value and the expected sample sizes for each group given current means and standard deviations
}
\description{
Given two vectors of continuous normal distributions, this function takes the
information from the two, calculates the expected number of participants to be allocated in each group,
and gives the allocation of the next participant of the study.
}
\examples{
IMAD_next_val_cont(rnorm(20, 1, 2), rnorm(20, 0, 2), nu = 0.8, weight = 0.7)
}
