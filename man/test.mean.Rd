\name{test.mean}
\alias{test.mean}
\title{\code{test.mean} conducts a hypothesis test based on a sample mean.}
\usage{
test.mean(pe, nullval, sd, n, conflvl, twosided = F)
}
\arguments{
  \item{pe}{(point estimate) is sample mean}

  \item{nullval}{is the null value}

  \item{sd}{is standard deviation}

  \item{n}{is sample size}

  \item{conflvl}{is confidence level (0-1)}

  \item{twosided}{is TRUE if test is twosided, defaults to
  FALSE}
}
\value{
list contains values for Z score, p-value and boolean for
reject/not reject null hypothesis
}
\description{
\code{test.mean} conducts a hypothesis test based on a
sample mean.
}
\examples{
test.mean(pe=1.1, nullval=0, sd=4.9, n=51, conflvl=0.95, twosided=F)
}

