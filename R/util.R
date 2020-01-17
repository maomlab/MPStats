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
#' @return variance of parameter estimate
binomial_variance <- function(n_positive, n_trials){
  m_negative <- n_trials - n_positive
  (n_positive+1)*(m_negative+1)/((n_positive + m_negative + 2)^2 * (n_positive + m_negative + 3))
}



#' @importFrom magrittr %>%
NULL