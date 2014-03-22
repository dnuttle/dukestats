\name{prob.negbinom}
\alias{prob.negbinom}
\title{\code{prob.negbinom} returns the probability for a negative binomial distribution.
It returns the probability of of the first success being the kth trial
in n trials, given p probability on any trial.}
\usage{
prob.negbinom(n, k, p)
}
\arguments{
  \item{n}{is the number of trials}

  \item{k}{is the number of successes}

  \item{p}{is the probability of a success on any trial}
}
\value{
a single value, the probability
}
\description{
\code{prob.negbinom} returns the probability for a negative
binomial distribution. It returns the probability of of the
first success being the kth trial in n trials, given p
probability on any trial.
}
\examples{
prob.negbinom(n=10, k=3, p=0.5)
}

