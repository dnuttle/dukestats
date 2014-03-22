\name{bayes}
\alias{bayes}
\title{\code{bayes} accepts a list representing a decision tree of probabilities,
as well as two values, which determine the first and second paths to take in the tree.}
\usage{
bayes(data, a, b)
}
\arguments{
  \item{data}{is a list whose elements represent all of the
  paths in the first level of the decision tree.  Each of
  those elements must be a list, whose elements in turn
  must be two values named p1 and p2.  The p1 value is the
  probability of the path at the first level of the tree,
  and p2 is the probability of the path on the second level
  of the tree. The p1 values of all elements in data should
  add up to 1.}

  \item{a}{is the path taken on the first level of the
  tree.  This is a number that should be between 1 and the
  number of elements in data.}

  \item{b}{is the path taken on the second level of the
  tree.  This must be 1 or 2.}
}
\value{
probability value
}
\description{
\code{bayes} accepts a list representing a decision tree of
probabilities, as well as two values, which determine the
first and second paths to take in the tree.
}
\examples{
bayes(data = list(path1=list(sick=0.05, positive=0.95), path2=list(healthy=0.95, positive=0.35)), a=1, b=1) # probability of disease, given a positive result (returns 0.125)
}

