#' Allocation using IMAD for Binary Variables
#'
#' @param group_a A vector containing the values of the first group
#' @param group_b A vector containing the values of the second group
#' @param trtlvls A vector containing the labels of each respective group
#' @param p_clin A vector that takes the minimum clinically significant success probabilities of each group
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param alpha The Type 1 error rate
#' @param beta The Type 2 error rate
#'
#' @return The grouping of the next value and the expected sample sizes for each group given current estimated probabilities
#' @export
#'
#' @examples IMAD_next_val_bin()
IMAD_next_val_bin <- function(group_a, group_b, trtlvls = c("A", "B"), p_clin, weight = 0.5, alpha = 0.05, beta = 0.2){


  #Get nu from clinician expected difference;
  nu = p_clin[1] - p_clin[2]

  #Get z values for alpha, beta, and combined
  zalpha <- qnorm(alpha/2)
  zbeta <- qnorm(beta)
  zc <- zalpha + zbeta

  #CME
  cme <- (nu / (zc)) ^ 2

  #Initial Sample sizes
  n_a <- length(group_a)
  n_b <- length(group_b)
  n <- n_a + n_b

  #Set initial odds ratio
  or_nought <- (p_clin[1] * (1 - p_clin[2])) / (p_clin[2] * (1 - p_clin[1]))

  #Set pi_nought;
  #pi_null <- pi_nought(pa_clin = p_clin[1], pb_clin = p_clin[2], weight = weight, cme = cme)
  pi_null = 0.5

  #Set probability of a
  pa_hat <- sum(group_a) / length(group_a)
  pa_hat <- ifelse(pa_hat == 1, .999, pa_hat)
  pa_hat <- ifelse(pa_hat == 0, .001, pa_hat)
  qa_hat <- 1- pa_hat

  #Set probability of b
  pb_hat <- sum(group_b) / length(group_b)
  pb_hat <- ifelse(pb_hat == 1, .999, pb_hat)
  pb_hat <- ifelse(pb_hat == 0, .001, pb_hat)
  qb_hat <- 1 - pb_hat

  #Set odds ratios of a and b
  or_a <- pa_hat * qb_hat / pb_hat * qa_hat
  or_b <- pb_hat * qa_hat / pa_hat * qb_hat

  #Set nj
  nj = (1 / log(1/or_nought ^ 2)) * 2 * log(sqrt((pa_hat * qa_hat) / (pb_hat * qb_hat)) * ((1 / pi_null) - 1))

  #Get function needed for a and b
  fa <- exp((nj) * log(1 / or_a))
  fb <- exp((nj) * log(1 / or_b))


  #Get weighted objective function
  wo <- weighted_objective_bin(pa_clin = p_clin[1], pb_clin = p_clin[2], pa_hat = pa_hat, pb_hat = pb_hat,
                               fa = fa, fb = fb, weight = weight, n_a = n_a, n_b = n_b,
                               n = n, zalpha = zalpha, zbeta = zbeta)

  allocation <- wo$allocation

  #Set correct allocation group
  grouping <- ifelse(allocation == 0, trtlvls[1], trtlvls[2])

  #Export grouping and estimated sample sizes
  outcomes <- list(grouping, wo$estimated_na, wo$estimated_nb, wo$estimated_n)
  names(outcomes) <- c("allocation_group", "estimated_groupa", "estimated_groupb", "estimated_total")

  #Export proper values;
  return(outcomes)
}
