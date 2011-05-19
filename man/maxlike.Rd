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

# Fake data based on R's volcano dataset

data(MaungaWhau)
elev <- raster(MaungaWhau$elev, 0, 61, 0, 87)
precip <- raster(MaungaWhau$precip, 0, 61, 0, 87)
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



# Simulation example

x1 <- sort(rnorm(100))
x1 <- raster(outer(x1, x1), xmn=0, xmx=100, ymn=0, ymx=100)

x2 <- raster(matrix(runif(1e4), 100, 100), 0, 100, 0, 100)

logit.psi <- -1 + 1*x1 + 0*x2
psi <- exp(logit.psi)/(1+exp(logit.psi))
plot(psi)

r <- stack(x1, x2)
r@layernames <- c("x1", "x2")
plot(r)

set.seed(131)
pa <- matrix(NA, 100, 100)
pa[] <- rbinom(1e4, 1, as.matrix(psi))
str(pa)
table(pa)

pa <- raster(pa, 0, 100, 0, 100)
plot(pa)

xy <- xyFromCell(pa, sample(Which(pa==1, cells=TRUE), 1000))

plot(x1)
points(xy)

maxlike(~x1 + x2, raster=r, points=xy)
maxlike(~x1, r, xy)

}

