\name{maxlike}
\alias{maxlike}
\title{
  Estimate probability of occurrence from presence-only data
}
\description{
  Estimate the probability of occurrence from presence-only data
  (ie, no observed absences) using covariate data
  from a random sample of presence points, and covariate data from
  the entire region of interest.
}
\usage{
  maxlike(formula, raster, points, starts)
}
\arguments{
  \item{formula}{
     A right-hand side formula describing the model
}
  \item{raster}{
    A `raster stack' object created by the \code{\link{stack}} function
    in the \code{\link{raster-package}}.
}
  \item{points}{
    A matrix of coordinates for the `presence' locations.
}
  \item{starts}{
    Parameter starting values
}
}
\details{
  Data from the entire region of interest are contained in the stack of
  rasters. Data from the presence locations are extracted using the
  coordinates in \code{points}. These should both use the same
  coordinate system.
}
\value{
  A list with 3 components
  \item{Est}{data.frame containing the parameter estimates (Ests) and
    standard errors (SE).}
  \item{vcov}{variance-covariance matrix}
  \item{AIC}{AIC}
}
\note{
  In general it is very hard to obtain a random sample of presence
  points. This could be done using radio telemetry, but other sources
  of information such as museum specimens will almost surely include
  some sort of sampling bias, which will invalidate the results of
  presence-only models.
}


\seealso{ \code{\link{maxlike-package}}, \code{\link{raster}} }

\examples{

data(MaungaWhau)

elev <- MaungaWhau$elev
precip <- MaungaWhau$precip
xy <- MaungaWhau$xy

# Elevation, standardized
plot(elev)
points(xy)

# Precipitation, standardized
plot(precip)
points(xy)

# Stack them and make sure they are named
ep <- stack(elev, precip)
layerNames(ep) <- c("elev", "precip")

# Fit a model
fm <- maxlike(~elev + I(elev^2) + precip, ep, xy)
fm

# Expected probability of occurrence
beta <- coef(fm)
E.logitpsi <- beta[1] + beta[2]*elev + beta[3]*(elev^2) + beta[4]*precip
E.psi <- exp(E.logitpsi)/(1+exp(E.logitpsi))
plot(E.psi)
points(xy)

}

