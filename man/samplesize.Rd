\name{samplesize}
\alias{samplesize}
\title{\code{samplesize} calculates the sample size needed for a hypothesis test,
given the desired confidence level, the standard deviation, and the margin of error.}
\usage{
samplesize(conflvl, sd, moe)
}
\arguments{
  \item{conflvl}{is the desired confidence level (0-1)}

  \item{sd}{is the standard deviation}

  \item{moe}{is the desired margin of error}

  \item{return}{sample size}
}
\description{
\code{samplesize} calculates the sample size needed for a
hypothesis test, given the desired confidence level, the
standard deviation, and the margin of error.
}
\examples{
samplesize(0.95, 25, 4) # returns 150.057
}

