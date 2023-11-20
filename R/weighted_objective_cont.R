#' Weighted Objective Functions for Continuous Outcomes
#'
#' Given the standard deviations, positive functions of theta, and sample size of two treatment arms along with the weighting and minimally clinically meaningful difference,
#' the optimal sample size, allocation ratio, and sample size ratio are found using this function. From these the minimum sample size of each group, the biased coin probability,
#' and the allocation of the next participant from tht probability are given.
#'
#' @param sda The standard deviation of the first group
#' @param sdb The standard deviation of the second group
#' @param fa Location-invariant value of the first group
#' @param fb Location invariant value of the second group
#' @param cme The set constraint to ensure the difference is at least as great as the clinically meaningful difference
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param n_a The number of individuals currently in the first group
#' @param n_b The number of individuals currently in the second group
#' @param n The total number of individuals currently in either group
#' @param zalpha The standardized set type 1 error rate
#' @param zbeta The standardized set type 2 error rate
#' @param nu The set clinically meaningful difference
#'
#' @return Returns optimal allocation, optimal sample size ratio, and allocation of the next patient for continuous normal distributions
#' @export
#'
#' @examples weighted_objective_cont()
weighted_objective_cont <- function(sda = 1, sdb = 1, fa = .5, fb = .5, cme = 3.765,
                               weight = .5, n_a = 10, n_b = 10, n = 20,
                               zalpha = qnorm(0.05/2), zbeta = qnorm(0.2), nu = 1){
  #Set Weights;
  w1 = weight
  w2 = 1 - w1

  #Lambda Functions;
  lambda1 <- ((sda * sqrt(fa) + sdb * sqrt(fb)) / cme) ^ 2
  lambda2 <- ((sda + sdb) / cme) ^ 2

  #Objective function for RAR (phi1) and SSR(phi2)
  phi1 <- sqrt(lambda1) * sda * sqrt(fa) + sqrt(lambda1) * sdb * sqrt(fb)
  phi2 <- lambda2 * sda + lambda2 * sdb

  #Combine all elements to identify lambda for sample size of a and b
  lambda <- ((sda * sqrt((w1 * fa) / phi1 + w2 / phi2) + sdb * sqrt((w1 * fb) / phi1 + w2 / phi2)) / cme) ^ 2

  #Update optimal sample sizes for a and b
  naopt <- sqrt(lambda) * (sda / (sqrt((w1 * fa) / phi1 + w2 / phi2)))
  nbopt <- sqrt(lambda) * (sdb / (sqrt((w1 * fb) / phi1 + w2 / phi2)))

  #Updated allocation ratio;
  piopthat <- naopt / (naopt + nbopt)

  #Updated sample size ratio;
  rhoopthat <- nbopt / naopt

  #Re-Estimate sample size required for each group and total sample size;
  nahat <- ceiling(((zalpha + zbeta) ^ 2 * ((sdb ^ 2 / sda ^ 2) + rhoopthat) * (sda ^ 2)) / (rhoopthat * nu ^ 2) +
                     (((sdb ^ 2 / sda ^ 2) + rhoopthat ^ 3) * zalpha ^ 2) / (2 * rhoopthat * (sdb ^ 2 / sda ^ 2 + rhoopthat) ^ 2))
  nbhat <- ceiling(rhoopthat * nahat)
  nj <- nahat + nbhat

  #Get ERADE optimal Method;
  gamma <- 0.5
  piratio <- n_a / n
  piadjhat <- ifelse(piratio > piopthat, gamma * piopthat, ifelse(piopthat > piratio, 1 - gamma * (1 - piopthat), piopthat))

  #Get randomized next allocation;
  allocation <- ifelse(runif(1) < piadjhat, 1, 0)

  #Return necessary outcomes
  outcomes <- list(allocation, piadjhat, nj, lambda1, lambda2, phi1, phi2, lambda, nahat, nbhat, naopt, nbopt)
  names(outcomes) <- c("allocation", "pi", "estimated_n", "lambda1", "lambda2", "phi1", "phi2", "lambda", "estimated_na", "estimated_nb", "naopt", "nbopt")

  return(outcomes)
}
