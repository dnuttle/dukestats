\name{ci.mean}
\alias{ci.mean}
\title{\code{ci.mean} calculates and returns a confidence interval as a 2-element vector}
\usage{
ci.mean(pe, sd, n, conflvl)
}
\arguments{
  \item{pe}{is the point estimate (the sample mean)}

  \item{sd}{is the standard deviation}

  \item{n}{is the sample size}

  \item{conflvl}{is the desired confidence level (0-1)}
}
\value{
vector containing upper and lower bounds of confidence
interval
}
\description{
\code{ci.mean} calculates and returns a confidence interval
as a 2-element vector
}
\examples{
ci.mean(3, 2, 25, 0.95)
}

