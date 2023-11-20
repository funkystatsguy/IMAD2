#' Simulation Code for IMAD Continuous Case
#'
#' Given a pair of preset means, standard deviations, and initial sample sizes along with weight and the minimally meaningful clinical difference,
#' a single simulation result of what a theoretical trial would look like along with
#'
#' @param mu A vector giving the expected means of two groups
#' @param sig  A vector giving the expected standard deviations of two groups
#' @param n A vector giving the initial counts given to each group
#' @param weight A weighting of Sample Size Re-estimation and Random Adaptive Re-Allocation
#' @param alpha The Type 1 error rate
#' @param beta The Type 2 error rate
#' @param nu  The preset clinically meaningful difference to be looked for
#'
#' @return Gives the total number of patients to be allocated in each group as well as total sample size for that iteration of the simulation
#' @export
#'
#' @examples IMAD_simulation_cont(mu = c(2, 0), sig = c(3, 3),  nu = 1.3, n = c(20, 20), weight = 0.3)
IMAD_simulation_cont <- function(mu, sig, n, nu, weight = 0.5, alpha = 0.05, beta = 0.2){
  #Bring in Random Number Generator
  random_norm <- random_normal(mu = mu, sig = sig, n = n)

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

  #Set values from random_norm
  #Intervention
  a <- random_norm$a
  #Control
  b <- random_norm$b

  #Initial Sample Sizes
  n <- random_norm$ninit
  n_a <- random_norm$nainit
  n_b <- random_norm$nbinit

  #Get random mu values
  mu = mu
  sig = sig


  repeat{

    #Get intervention values
    if(n_a == random_norm$nainit){
      a = random_norm$a
    } else {
      a = newval$a
    }

    #Get control values
    if(n_b == random_norm$nbinit){
      b = random_norm$b
    } else {
      b = newval$b
    }

    #Get n of a and b
    n_a = as.numeric(length(a))
    n_b = as.numeric(length(b))
    n = n_a + n_b
    n1 <- n + 1

    #Mean of Intervention
    meana <- mean(a)
    #Mean of Control
    meanb <- mean(b)

    #SD of Intervention
    sda <- sd(a)
    #SD of Control
    sdb <- sd(b)

    #Get estimated functions for a and b;
    #Identify z-scores for each variable;
    faz <- (meana - meanb) / (n1 * sqrt(sda^2 + sdb^2))
    fab <- (meanb - meana) / (n1 * sqrt(sda^2 + sdb^2))

    #Get normal CDF scores for a and b;
    fa <- pnorm(faz)
    fb <- pnorm(fab)

    #Get weighted objective function
    wo <- weighted_objective_cont(sda = sda, sdb = sdb, fa = fa, fb = fb, cme = cme,
                             n = n, n_a = n_a, n_b = n_b, zalpha = zalpha, zbeta = zbeta,
                             weight = weight, nu = nu)

    allocation <- wo$allocation

    #Allocation assignment
    if(allocation == 1){
      a = as.matrix(c(a, rnorm(1, mu[1], sig[1])))
    } else {
      a = a
    }

    if(allocation == 0){
      b = as.matrix(c(b, rnorm(1, mu[2], sig[2])))
    } else {
      b = b
    }

    #a <- ifelse(allocation == 1, c(a, rnorm(1, mu[1], sig[1])), a)
    #b <- ifelse(allocation == 0, c(b, rnorm(1, mu[2], sig[2])), b)6
    n_a <- as.numeric(length(a))
    n_b <- as.numeric(length(b))
    n <- n_a + n_b

    #Get new values
    newval <- list(a, b)
    names(newval) <- c("a", "b")

    #Set what we would like to return as output
    out <- list(length(a), length(b), n, wo$estimated_n)
    names(out) <- c('n_a', "n_b", "n", "est_n")

    if (n >= wo$estimated_n){
      return(out)
    }
  }
}
