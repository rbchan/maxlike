\name{maxlike-package}
\alias{maxlike-package}
\docType{package}
\title{
  Model species distributions by estimating the probability of occurrence
  using presence-only data.
}
\description{
  A species distribution can be characterized by the probability that it
  occurs at some location in space. This can be easily accomplished
  using presence-absence data, but often researchers only have presence
  locations and environmental data for the study area. MAXENT is a
  popular software program for modeling species distributions, but it does
  not estimate the probability of occurrence. Rather it returns various
  indices that are not easy to interpret (see Royle et al. In
  Press). Package "maxlike" provides a simple
  likelihood-based alternative.
}
\note{
  All presence-only models require a random sample of data from
  locations where the species is present. Unfortunately, random sampling
  is not a feature of most presence-only datasets, and bias in the
  estimated probability surface should
  be expected in such cases. This assumption can be greatly relaxed if
  one has presence-absence data, which will always
  contain more information about a species' distribution.
  }
\details{
\tabular{ll}{
Package: \tab maxlike\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2011-12-20\cr
License: \tab GPL (>=3) \cr
LazyLoad: \tab yes\cr
}
}
\references{
  J. A. Royle, R. B. Chandler, C. Yackulic and J.
  D. Nichols. In Press. Likelihood analysis of species occurrence probability
  from presence-only data for modelling species distributions. Methods
  in Ecology and Evolution.
}

\seealso{ \code{\link{maxlike}}, \code{\link{raster}} }

\keyword{ package }
