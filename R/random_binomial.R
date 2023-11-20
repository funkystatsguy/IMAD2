#' Get Random Output from Two Binomial Distributions
#'
#' @param n A vector containing the number of samples allocated to each group
#' @param p A vector containing the success probabilities of each group
#'
#' @return The binomial results from each group, the number of samples in each group, and the probabilities assigned to each group
#' @export
#'
#' @examples random_binomial()
#' @examples random_binomial(n = c(10, 10), p = c(0.9, 0.1))
random_binomial <- function(n = c(20, 20), p = c(0.5, 0.5)){

  #Set Probabilities of each group
  pa <- p[1]
  pb <- p[2]

  #Set initial burn in values
  nainit <- n[1]
  nbinit <- n[2]

  #Get total sample size at initial iteration
  ninit <- nainit + nbinit

  #Get initial random values
  a <- as.matrix(replicate(nainit, rbinom(1, 1, pa)))
  b <- as.matrix(replicate(nbinit, rbinom(1, 1, pb)))

  #Put together list of vectors
  t <- list(a, b, nainit, nbinit, ninit, p)
  names(t) <- c("a", "b", "nainit", "nbinit", "ninit", "p")

  return(t)

}
