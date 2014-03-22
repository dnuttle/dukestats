library(stats)

#' \code{ci.mean} calculates and returns a confidence interval as a 2-element vector
#' @param pe is the point estimate (the sample mean)
#' @param sd is the standard deviation
#' @param n is the sample size
#' @param conflvl is the desired confidence level (0-1)
#' @return vector containing upper and lower bounds of confidence interval
#' @export
#' @examples
#' ci.mean(3, 2, 25, 0.95)
ci.mean <- function(pe,sd,n,conflvl){
  alpha = 1 - conflvl
  zstar = qnorm(1-(alpha/2))
  loci = pe - zstar * sd/sqrt(n)
  upci = pe + zstar * sd/sqrt(n)  	
  return(c(loci , upci))
}

#' \code{test.mean} conducts a hypothesis test based on a sample mean.
#' @param pe (point estimate) is sample mean
#' @param nullval is the null value
#' @param sd is standard deviation
#' @param n is sample size
#' @param conflvl is confidence level (0-1)
#' @param twosided is TRUE if test is twosided, defaults to FALSE
#' @return list contains values for Z score, p-value and boolean for reject/not reject null hypothesis
#' @export
#' @examples
#' test.mean(pe=1.1, nullval=0, sd=4.9, n=51, conflvl=0.95, twosided=F)
test.mean <- function(pe, nullval, sd, n, conflvl, twosided=F) {
  alpha <- 1 - conflvl
  Z = qnorm(1-(alpha/2))
  pvalue = 1 - pnorm(Z)
  pvalue = pvalue.mean(pe, nullval, sd, twosided)
  reject <- (pvalue < alpha)
  result <- list(Z=Z, pvalue=pvalue, reject=reject)
  return(result)
}

#' \code{prob.negbinom} returns the probability for a negative binomial distribution. It returns the probability of k successes in n trials, given p probability on any trial.
#' @param n is the number of trials
#' @param k is the number of successes
#' @param p is the probability of a success on any trial
#' @return a single value, the probability
#' @export
#' @examples
#' prob.negbinom(n=10, k=3, p=0.5)
prob.negbinom <- function(n, k, p) {
  return(choose(n-1, k-1) * p^k * (1-p)^(n-k))
}


#' \code{pvalue.mean} returns the p-value for a point estimate that is a sample mean.
#' @param pe is the point estimate
#' @param nullval is the null value being compared to the point estimate
#' @param sd is the standard deviation
#' @param twosided is a binary value that determines whether to return a value based on a two-sided test (defaults to FALSE)
#' @return p-value
#' @export
#' @examples
#' pvalue.mean(pe=1.1, nullval=0, sd=4.9, n=51, twosided=F) # Returns 0.0544477
pvalue.mean <- function(pe, nullval, sd, n, twosided=F) {
  SE = sd / sqrt(n)
  pv = 1 - pnorm(abs(pe - nullval)/SE)
  if(twosided) {
    pv = pv * 2
  }
  return(pv)
}

#' \code{bayes} accepts a list representing a decision tree of probabilities,
#' as well as two values, which determine the first and second paths to take in the tree.
#' @param data is a list whose elements represent all of the paths in the first level of
#' the decision tree.  Each of those elements must be a list, whose elements in turn
#' must be two values named p1 and p2.  The p1 value is the probability of the path
#' at the first level of the tree, and p2 is the probability of the path on the
#' second level of the tree.
#' The p1 values of all elements in data should add up to 1.
#' @param a is the path taken on the first level of the tree.  This is a number that should be between 1
#' and the number of elements in data.
#' @param b is the path taken on the second level of the tree.  This must be 1 or 2.
#' @return probability value
#' @export
#' @examples
#' bayes(data = list(path1=list(sick=0.05, positive=0.95), path2=list(healthy=0.95, positive=0.35)), a=1, b=1) # probability of disease, given a positive result (returns 0.125)
bayes <- function(data, a, b) {
  p2 = ifelse(b==1, data[[a]][[2]], 1-data[[a]][[2]])
  p1 = data[[a]][[1]]
  numerator = p1 * p2
  denominator = 0
  for(i in 1:length(data)) {
    p2 = ifelse(b==1, data[[i]][[2]], 1-data[[i]][[2]])
    
    denominator = denominator + data[[i]][[1]] * p2
  }
  return(numerator/denominator)
}

#' \code{samplesize} calculates the sample size needed for a hypothesis test, 
#' given the desired confidence level, the standard deviation, and the margin of error.
#' @param conflvl is the desired confidence level (0-1)
#' @param sd is the standard deviation
#' @param moe is the desired margin of error
#' @param return sample size
#' @export
#' @examples
#' samplesize(0.95, 25, 4) # returns 150.057
samplesize <- function(conflvl, sd, moe) {
  alpha = 1 - conflvl
  zstar = qnorm(1 - (alpha/2))
  return(((zstar * sd)/moe)^2)
}

