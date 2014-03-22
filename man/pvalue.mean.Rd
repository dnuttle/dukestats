\name{pvalue.mean}
\alias{pvalue.mean}
\title{\code{pvalue.mean} returns the p-value for a point estimate that is a sample mean.}
\usage{
pvalue.mean(pe, nullval, sd, n, twosided = F)
}
\arguments{
  \item{pe}{is the point estimate}

  \item{nullval}{is the null value being compared to the
  point estimate}

  \item{sd}{is the standard deviation}

  \item{twosided}{is a binary value that determines whether
  to return a value based on a two-sided test (defaults to
  FALSE)}
}
\value{
p-value
}
\description{
\code{pvalue.mean} returns the p-value for a point estimate
that is a sample mean.
}
\examples{
pvalue.mean(pe=1.1, nullval=0, sd=4.9, n=51, twosided=F) # Returns 0.0544477
}

