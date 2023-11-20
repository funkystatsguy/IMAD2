#' Title
#'
#' @param n A vector allocating initial counts to each group
#' @param p A vector that takes the actual simulated difference between two groups
#' @param p_clin A vector that takes the minimum clinically significant success probabilities of each group
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param alpha The Type 1 error rate
#' @param beta The Type 2 error rate
#'
#' @return The number of individuals allocated to each group, total sample size, and minimum estimated sample size
#' @export
#'
#' @examples IMAD_simulation_bin()
IMAD_simulation_bin <- function(n = c(20, 20), p = c(0.6, 0.5), p_clin = c(0.6, 0.5), weight = 0.5, alpha = 0.05, beta = 0.8){
  #Bring in random binom
  random_bin <- random_binomial(n = n, p = p)

  #Set values from random binomial
  a <- random_bin$a
  b <- random_bin$b

  #Get nu from clinician expected differene;
  nu = p_clin[1] - p_clin[2]

  #Get z values for alpha, beta, and combined
  zalpha <- qnorm(alpha/2)
  zbeta <- qnorm(beta)
  zc <- zalpha + zbeta

  #CME
  cme <- (nu / (zc)) ^ 2

  #Initial Sample sizes
  n <- random_bin$ninit
  n_a <- random_bin$nainit
  n_b <- random_bin$nbinit

  #Get random pi values
  p <- p

  #Set initial odds ratio
  or_nought <- (p_clin[1] * (1 - p_clin[2])) / (p_clin[2] * (1 - p_clin[1]))

  #Set pi_knot;
  #pi_null <- pi_nought(pa_clin = p_clin[1], pb_clin = p_clin[2], weight = weight, cme = cme)
  pi_null <- 0.5

  repeat{
    #Get intervention values
    if(n_a == random_bin$nainit){
      a = random_bin$a
    } else {
      a = newval$a
    }

    #Get control values
    if(n_b == random_bin$nbinit){
      b = random_bin$b
    } else {
      b = newval$b
    }

    #Get n of a and b
    n_a = as.numeric(length(a))
    n_b = as.numeric(length(b))
    n = n_a + n_b

    #Set probability of a
    pa_hat <- sum(a) / length(a)
    pa_hat <- ifelse(pa_hat == 1, .99, pa_hat)
    qa_hat <- 1- pa_hat

    #Set probability of b
    pb_hat <- sum(b) / length(b)
    pb_hat <- ifelse(pb_hat == 1, .99, pb_hat)
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

    #Allocation assignment
    if(allocation == 0){
      a = as.matrix(c(a, rbinom(1, 1, p[1])))
    } else {
      a = a
    }

    if(allocation == 1){
      b = as.matrix(c(b, rbinom(1, 1, p[2])))
    } else {
      b = b
    }

    #Get new length values of a and b
    n_a <- as.numeric(length(a))
    n_b <- as.numeric(length(b))
    n <- n_a + n_b

    #Get new values
    newval <- list(a, b, n , wo$estimated_n)
    names(newval) <- c("a", "b", "n", "est_n")

    #Set what we would like to return as output
    out <- list(length(a), length(b), n, wo$estimated_n)
    names(out) <- c('n_a', "n_b", "n", "est_n")



    if (n >= wo$estimated_n){
      return(out)
    }
  }
}
