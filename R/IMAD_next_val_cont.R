#' Allocation using IMAD for Continuous Variables
#'
#' Given two vectors of continuous normal distributions, this function takes the
#' information from the two, calculates the expected number of participants to be allocated in each group,
#' and gives the allocation of the next participant of the study.
#'
#'
#'
#' @param group_a A vector containing the values of the first group
#' @param group_b A vector containing the values of the second group
#' @param trtlvls A vector containing the labels of each respective group
#' @param nu The preset clinically meaningful difference to be looked for
#' @param pinought The target allocation ratio for the better performing group
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param alpha The Type 1 error rate
#' @param beta The Type 2 error rate
#'
#' @return The grouping of the next value and the expected sample sizes for each group given current means and standard deviations
#' @export
#'
#' @examples IMAD_next_val_cont(rnorm(20, 1, 2), rnorm(20, 0, 2), nu = 0.8, weight = 0.7)

IMAD_next_val_cont <- function(group_a, group_b, trtlvls = c("A", "B"), nu, pinought = 0.5, weight = 0.5, alpha = 0.05, beta = 0.2){
  #Standard Deviation Calculation
  #Alpha deviation;
  zalpha <- qnorm(alpha/2)
  #Beta deviation
  zbeta <- qnorm(beta)
  #Aggregate Deviation
  zc <- zalpha + zbeta
  #Standardize Clinically meaningful effect size
  #BE SURE TO RELABEL
  cme <- (nu / zc) ^ 2

  #Initial Sample Sizes
  n_a <- length(group_a)
  n_b <- length(group_b)
  n <- n_a + n_b

  #Mean of Intervention
  meana <- mean(group_a)
  #Mean of Control
  meanb <- mean(group_b)

  #SD of Intervention
  sda <- sd(group_a)
  #SD of Control
  sdb <- sd(group_b)

  #Set optimization function for each iteration;
  optfunc <- function(nu, sda, sdb, pinought, par){
    faz <- (-nu) / (par * sqrt(sda^2 + sdb^2))
    fbz <- (nu) / (par * sqrt(sda^2 + sdb^2))
    (pinought - (sda * sqrt(pnorm(fbz)))/ (sda * sqrt(pnorm(fbz)) + sdb * sqrt(pnorm(faz))))^2
  }

  #Get optimization for given pi_nought;
  optval <- optim(par = 1, optfunc, nu = nu, sda = sda,
                  sdb = sdb, pinought = pinought, method = "Brent",
                  lower = -10000, upper = 10000)

  #Get eta_j for each iteration;
  eta_j <- optval$par

  #Get estimated functions for a and b;
  #Identify z-scores for each variable;
  faz <- (meana - meanb) / (eta_j * sqrt(sda^2 + sdb^2))
  fab <- (meanb - meana) / (eta_j * sqrt(sda^2 + sdb^2))

  #Get normal CDF scores for a and b;
  fa <- pnorm(faz)
  fb <- pnorm(fab)

  #Get weighted objective function
  wo <- weighted_objective_cont(sda = sda, sdb = sdb, fa = fa, fb = fb, cme = cme,
                                n = n, n_a = n_a, n_b = n_b, zalpha = zalpha, zbeta = zbeta,
                                weight = weight, nu = nu)

  allocation <- wo$allocation
  grouping <- ifelse(allocation == 1, trtlvls[1], trtlvls[2])

  #Get Outcome Values
  outcome <- list(grouping, wo$pi, wo$estimated_na, wo$estimated_nb, wo$estimated_n)
  names(outcome) <- c("allocation_group", "allocation_ratio", "estimated_groupa", "estimated_groupb", "estimated_total")

  #Return value
  return(outcome)
}
