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
  maxlike(formula, covData, ptData, removeDuplicates = FALSE, starts,
          hessian = TRUE, na.action = "na.omit", ...)
}
\arguments{
  \item{formula}{
    A right-hand side formula describing the model
    }
  \item{covData}{
    Either a `raster stack' created by the \code{\link{stack}} function
    in the \code{\link{raster-package}}, or a \code{\link{data.frame}}
    in which each row is a pixel with covariate data.
    }
  \item{ptData}{
    If covData is a raster \code{\link{stack}}, ptData must be a
    \code{\link{data.frame}} with the X and Y coordinates. If covData is
    a \code{\link{data.frame}}, ptData must be a data.frame or a matrix
    with the same column names as covData. In this case, each row is
    covariate data for the presence locations.
    }
  \item{removeDuplicates}{
    Logical. If there is more than one point in a pixel, should the
    duplicates be removed? This is only applicable if covData is a
    raster \code{\link{stack}}. If it is a \code{\link{data.frame}}, you
    should handle duplicates yourself.
    }
  \item{starts}{
    Parameter starting values
    }
  \item{hessian}{
    Should the hessian be computed and the variance-covariance matrix
    returned?
    }
  \item{na.action}{
    See \code{\link{options}} for choices
  }
  \item{...}{
    Additional arguments passed to \code{\link{optim}}
    }
}
\details{
  Data from the entire region of interest are contained in covData.
  Data from the presence locations are extracted using the
  coordinates in \code{points} when covData is a raster
  \code{\link{stack}}. Otherwise, they should be provided directly as a
  \code{\link{data.frame}}. If using spatial data, covData and ptData
  should both use the same coordinate system.
  }
\value{
  A list with 3 components
  \item{Est}{data.frame containing the parameter estimates (Ests) and
    standard errors (SE).}
  \item{vcov}{variance-covariance matrix}
  \item{AIC}{AIC}
  \item{call}{the original call}
  \item{retained}{the points retained if duplicates were removed}
  \item{pts.removed}{The points removed due to missing values}
  \item{pix.removed}{The pixesls removed due to missing values}
  }
\note{
  In general it is very hard to obtain a random sample of presence
  points. This could be done using radio telemetry, but other sources
  of information such as museum specimens will almost surely include
  some sort of sampling bias, which will invalidate the results of
  presence-only models. Presence-absence data collected using random
  sampling is always preferable.
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

