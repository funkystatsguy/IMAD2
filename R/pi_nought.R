#' Title
#'
#' @param pa_clin The minimum clinically significant success probability of the first group
#' @param pb_clin The minimum clinically significant success probability of the second group
#' @param weight A weighting of Sample Size Re-estimation and Random Adaptive Re-Allocation
#' @param cme The set constraint to ensure the difference is at least as great as the clinically meaningful difference
#'
#' @return The optimal allocation ratio derived from the minimally significant difference in success probabilities between two groups
#' @export
#'
#' @examples pi_nought()
pi_nought <- function(pa_clin = 0.6, pb_clin = 0.5, weight = 0.5, cme = 0.001){

  #Weights
  w1 <- weight
  w2 <- 1 - w1

  #Set q values;
  qa_clin = 1 - pa_clin
  qb_clin = 1 - pb_clin

  #Lambda Functions
  lambda1 <- ((sqrt(pa_clin * qa_clin) * sqrt(fa) + sqrt(pb_clin * qb_clin) * sqrt(fb)) / cme ) ^ 2
  lambda2 <- ((sqrt(pa_clin * qa_clin) + sqrt(pb_clin * qb_clin)) / (cme)) ^ 2

  #Objective functions for RAR (phi1) and SSR(phi2)
  phi1 <- sqrt(lambda1) * sqrt(pa_clin * qa_clin) * sqrt(fa) + sqrt(lambda1) * (pb_clin * qb_clin) * sqrt(fb)
  phi2 <- sqrt(lambda2) * sqrt(pa_clin * qa_clin) + sqrt(lambda2) * sqrt(pb_clin * qb_clin)

  #Combine all elements to identify lambda for sample size of a and b
  lambda <- ((sqrt(pa_clin * qa_clin) * sqrt((w1 * fa) / phi1 + w2 / phi2) +
                sqrt(pb_clin * qb_clin) * sqrt((w1 * fa) / phi1 + w2 / phi2)) / cme) ^ 2

  #Get null optimal sample sizes for a and b
  nanull <- sqrt(lambda) * sqrt(pa_clin * qa_clin) / sqrt((w1 * fa) / phi1 + w2 / phi2)
  nbnull <- sqrt(lambda) * sqrt(pb_clin * qb_clin) / sqrt((w1 * fb) / phi1 + w2 / phi2)

  #Get null allocation ratio
  pinull <- nanull / (nanull + nbnull)

  #Export pinull
  return(pinull)

}
