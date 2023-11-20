#' Get Random Distribution from Two Normal Distributions
#'
#' @param n A vector containing the number of samples allocated to each group
#' @param mu A vector of means of the two normal distributions
#' @param sig A vector of standard deviations of the two normal distributions
#'
#' @return Results from the random distributions along with the preset values given
#' @export
#'
#' @examples random_normal()
#' @examples random_normal(n = c(10, 10), mu = c(-1, 0), sig = c(2, 2))
random_normal <- function(n = c(20, 20), mu = c(-10, 0), sig = c(1, 1)){

  #Randomize Variables
  #Set initial mean values
  mu_a <- mu[1]
  mu_b <- mu[2]

  #Set initial standard deviations
  sig_a <- sig[1]
  sig_b <- sig[2]

  #Set initial burn in values
  nainit <- n[1]
  nbinit <- n[2]

  #Get total sample size at initial iteration;
  ninit <- nainit + nbinit

  #Randomize set a and set b
  a <- as.matrix(rnorm(nainit, mu_a, sig_a))
  b <- as.matrix(rnorm(nbinit, mu_b, sig_b))

  #Put together list of vectors
  t <- list(a, b, nainit, nbinit, ninit, mu = mu, sig = sig)
  names(t) <- c("a", "b", "nainit", "nbinit", "ninit", "mu", "sig")

  #Return list of vectors;
  return(t)
}
