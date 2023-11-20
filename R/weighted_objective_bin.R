#' Weighted Objective Functions for Binary Outcomes
#'
#' @param pa_clin The minimum clinically significant success probability of the first group
#' @param pb_clin The minimum clinically significant success probability of the second group
#' @param pa_hat The current expected probability of the first group
#' @param pb_hat The current expected probability of the second group
#' @param fa Location invariant value of the first group
#' @param fb Location invariant value of the second group
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param n_a The number of individuals currently in the first group
#' @param n_b The number of individuals currently in the second group
#' @param n The total number of individuals currently in either group
#' @param zalpha The standardized set type 1 error rate
#' @param zbeta The standardized set type 2 error rate
#'
#' @return Returns optimal allocation, optimal sample size ratio, and allocation of the next patient for binary outcomes
#' @export
#'
#' @examples weighted_objective_bin()
weighted_objective_bin <- function(pa_clin = 0.6, pb_clin = 0.5, pa_hat = 0.6, pb_hat = 0.5, fa = 0.5,
                                   fb = 0.5, weight = 0.5, n_a = 10, n_b = 10, n = 20,
                                   zalpha = qnorm(0.05/2), zbeta = qnorm(0.2)){

  #Set nu;
  nu = pa_clin - pb_clin

  #Set CME again;
  cme = ( nu / (zalpha + zbeta)) ^ 2

  #Weights
  w1 <- weight
  w2 <- 1 - w1

  #Set q values;
  qa_hat = 1 - pa_hat
  qb_hat = 1 - pb_hat
  qa_clin = 1 - pa_clin
  qb_clin = 1 - pb_clin

  #Lambda Functions
  lambda1 <- ((sqrt(pa_hat * qa_hat) * sqrt(fa) + sqrt(pb_hat * qb_hat) * sqrt(fb)) / cme ) ^ 2
  lambda2 <- ((sqrt(pa_hat * qa_hat) + sqrt(pb_hat * qb_hat)) / (cme)) ^ 2

  #Objective functions for RAR (phi1) and SSR(phi2)
  phi1 <- sqrt(lambda1) * sqrt(pa_hat * qa_hat) * sqrt(fa) + sqrt(lambda1) * (pb_hat * qb_hat) * sqrt(fb)
  phi2 <- sqrt(lambda2) * sqrt(pa_hat * qa_hat) + sqrt(lambda2) * sqrt(pb_hat * qb_hat)

  #Combine all elements to identify lambda for sample size of a and b
  lambda <- ((sqrt(pa_hat * qa_hat) * sqrt((w1 * fa) / phi1 + w2 / phi2) +
                sqrt(pb_hat * qb_hat) * sqrt((w1 * fa) / phi1 + w2 / phi2)) / cme) ^ 2

  #Update optimal sample sizes for a and b
  naopt <- sqrt(lambda) * sqrt(pa_hat * qa_hat) / sqrt((w1 * fa) / phi1 + w2 / phi2)
  nbopt <- sqrt(lambda) * sqrt(pb_hat * qb_hat) / sqrt((w1 * fb) / phi1 + w2 / phi2)

  #Update allocation ratio
  piopthat <- naopt / (naopt + nbopt)

  #Updated sample size ratio
  rhoopthat <- nbopt / naopt

  #Re-estimate sample size required for each group and total sample size;
  #WE CANNOT USE PA_HAT OR PB_HAT HERE!!!!
  nahat <- ceiling(((zalpha + zbeta) ^ 2 / nu ^ 2) * ((pb_clin * qb_clin) / rhoopthat + pa_clin * qa_clin))
  nbhat <- ceiling(rhoopthat * nahat)
  nj <- nahat + nbhat

  #Get ERADE optimal method
  gamma <- 0.5
  piratio <- n_a / n
  piadjhat <- ifelse(piratio > piopthat, gamma * piopthat, ifelse(piopthat > piratio, 1 - gamma * (1 - piopthat), piopthat))

  #Get Randomized next allocation
  allocation <- ifelse(runif(1) < piadjhat, 1, 0)

  #Return necessary outcomes
  outcomes <- list(allocation, piadjhat, nj, nahat, nbhat)
  names(outcomes) <- c("allocation", "pi", "estimated_n", "estimated_na", "estimated_nb")

  return(outcomes)
}
