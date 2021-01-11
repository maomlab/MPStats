#' create date code of the form YYMMDD
#' inputs:
#'    d (optional): date code in the format of Sys.Date() for which to generate the date code, defaulting to 'today'
#' @export
date_code <- function(d=NA){
  # reference http://www.r-cookbook.com/node/17
  if(is.na(d)) d <- Sys.Date()
  pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
  paste(
    sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
    sep="")
}


#' Variance of bayesian estimator for binomial trial with flat prior
#' 
#' The bayesian estimator for the success probability p from a binomial trial
#' with n successes and m failures and beta prior with rate parameters a and b is
#'    posterior ~ prior * likelihood
#'    posterior ~ Beta(a,b) * Binomial(n, m+n)
#'    posterior ~ p^(a-1) * (1-p)^(b-1) * p^n * (1-p)^m
#'    posterior ~ p^(n+a-1) * (1-p)^(m+b-1)
#'    posterior ~ Beta(n+a,m+b)
#'    
#' A flat prior is a Beta(a=1, b=1), so the estimate of the posterior distribution
#' for p is Beta(n+1, m+1). One way to think of this is that since the beta distribution is
#' the congugate prior to the binomial distribution, the flat prior effectively adds one
#' pseudo positive count and one pseudo negative count to the observation.
#'
#' The variance of a Beta distribution is
#'
#'    var[Beta(a,b)] = a*b/((a+b)^2(a+b+1))
#'    
#' so the variance of the posteriror is
#' 
#'    var[posterior] = var[Beta(n+1,m+1)] = (n+1)(m+1)/((n+m+2)^2 * (n+m+3))
#' 
#' Ref: https://stats.stackexchange.com/questions/185221/binomial-uniform-prior-bayesian-statistics
#' 
#' @param n_positive number of successes in binomial trial
#' @param n_trials number of binomial trials
#' @param prior_positive number of pseudo successes to add as prior
#' @param prior_negative number of pseudo failures to add as prior
#' @return variance of parameter estimate
binomial_variance <- function(n_positive, n_trials, prior_positive=1, prior_negative=1){
  m_negative <- n_trials - n_positive
  n_positive <- n_positive + prior_positive
  m_negative <- m_negative + prior_negative
  numerator <- n_positive*m_negative
  denominator <- (n_positive + m_negative)^2 * (n_positive + m_negative + 1)
  numerator / denominator
}

#' Quantile of bayesian estimator for binomial trial with flat prior
#' 
#' Using the above derivation, the bayesian estimator for the success probability p
#' from a binomial trial with n successes and m failures and flat prior is Beta(n+1, m+1).
#' 
#' @param n_positive number of successes in binomial trial
#' @param n_trials number of binomial trials
#' @param p probability quantile to return
#' @export
binomial_quantile <- function(n_positive, n_trials, p){
  m_negative <- n_trials - n_positive
  qbeta(p=p, shape1=n_positive+1, shape2=m_negative+1)
}

#' Quantile of bayesian estimators for samples of a poisson process
#' 
#' The gamma distribution is the congugate prior for the poisson distribution. Given a prior `Gamma(k,theta)``,
#' and poisson observations `k_1, k_2, ..., k_n`, the bayesian estimator for the posterior is
#' 
#'      Gamma(k + sum k_i, n + theta)
#' 
#' So the prior can be thought of as adding `theta` pseudo observations each with a count value of `k/theta`.
#' 
#' Example:
#' 
#'     Say we have collected a sample of 3 counts with values `counts <- c(120, 130, 125)`
#'     then the 95% credible interval for the rate poisson rate parameter is
#'     
#'          low <- poisson_quantile(counts, .025)
#'          high <- poisson_quantile(counts, .975)
#'          
#'     Now say we have a prior equivilent to a single count with 124, then we would add 124 to the counts vector:
#'     
#'            counts <- c(120, 130, 125) + c(124)
#'            low <- poisson_quantile(counts, .025)
#'            high <- poisson_quantile(counts, .975)
#'     
#' @param count vector of n counts, assumed to be samples from a poisson process
#' @param probability quantile to return
#' @return count value of the given quantile of the posterior estimator for the poisson process
#'
#' @export
poisson_quantile <- function(count, p){
  qgamma(p=p, shape=sum(count), rate=length(count))
}

#' Univrsity of Michigan colors
#' 
#' https://med.umich.edu/branding/color.html
#' @export
umich_colors <- c(
  blue = "#00274C",
  maize = "#FFCB05",
  umma_tan = "#DDC497",
  soft_yellow = "#F2D384",
  barton_beige = "#979467",
  brown = "#7D6A56",
  lsa_orange = "#BA5827",
  cranberry = "#853754",
  matthaei_violet = "#514E86",
  south_u_blue = "#174992",
  grey_rock = "#878A8F",
  teal = "#20657E",
  wave_field_green = "#B3AC17",
  taubman_teal = "#00B5AF",
  rackham_roof_green = "#00B5AF",
  arboretum_blue = "#0174BB",
  dusty_blue = "#9FB6BC",
  grey_blue = "#465E85")
  


#' @importFrom magrittr %>%
NULL
