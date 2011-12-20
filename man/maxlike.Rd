\name{maxlike}
\alias{maxlike}
\title{
  Model probability of occurrence using presence-only data
}
\description{
  Model the probability of occurrence from presence-only data
  (i.e., datasets lacking observed absences) using covariate data
  from a random sample of presence points, and covariate data from
  the entire region of interest.
}
\usage{
  maxlike(formula, rasters, points, starts,
          hessian = TRUE, na.action = "na.omit", ...)
}
\arguments{
  \item{formula}{
    A right-hand side \code{\link{formula}} describing the model.
    At least 1 continuous covariate must be present in the formula.
    }
  \item{rasters}{
    A `raster stack' created by the \code{\link{stack}} function
    in the \code{\link{raster-package}}. It's a good idea to standardize
    these by subtracting the mean and dividing by the standard
    deviation. This will make it easier for \code{\link{optim}} to find
    the maximum-likelihood estimates.
    }
  \item{points}{
    A \code{\link{data.frame}} with the X and Y coordinates of the
    presence locations.
    }
  \item{starts}{
    Starting values for the parameters. This should be a vector with as
    many elements as there are parameters. By default, all starting
    values are 0, which should be adequate if covariates are standardized.
    }
  \item{hessian}{
    Logical. Should the hessian be computed and the variance-covariance matrix
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
  Data from the entire region of interest are contained in rasters.
  Data from the presence locations are extracted using the
  coordinates in \code{points}. Rasters and points
  should both use the same coordinate system. The program does not check
  this so it is up to the user.
  }
\value{
  A list with 3 components
  \item{Est}{data.frame containing the parameter estimates (Ests) and
    standard errors (SE).}
  \item{vcov}{variance-covariance matrix}
  \item{AIC}{AIC}
  \item{call}{the original call}
  \item{pts.removed}{The points removed due to missing values}
  \item{pix.removed}{The pixesls removed due to missing values}
  }
\note{
  In general it is very hard to obtain a random sample of presence
  points, which is a requirement of both the Royle et al. (In press)
  method and of MAXENT. This is one of many reasons why presence-absence
  data are preferable. When presence-absence data are available, they
  can be modeled using functions such as \code{\link{glm}}. Creating
  maps of prediction from \code{\link{glm}} can be accomplished in the
  same way as is shown in the example below.

  Another weakness of models for presence-only data is that they do not
  allow one to model detection probability, which is typically less than
  one in field conditions. If detection probability is affected by the same
  covariates that affect occurrence probability, then bias
  is inevitable. The R package \code{\link{unmarked}}
  (Fiske and Chandler 2011) offers numerous
  methods for jointly modeling both occurrence and detection probability
  when detection/non-detection data are available.
  }
\references{
  Royle, J.A., R.B. Chandler, C. Yackulic and J.
  D. Nichols. In Press. Likelihood analysis of species occurrence probability
  from presence-only data for modelling species distributions. Methods
  in Ecology and Evolution.

  Fiske, I. and R.B. Chandler. 2011. unmarked: An R Package for Fitting
  Hierarchical Models of Wildlife Occurrence and Abundance. Journal of
  Statistical Software 43(10).
}

\seealso{ \code{\link{maxlike-package}}, \code{\link{raster}} }

\examples{

# Fake data based on R's volcano dataset
# See help(raster) for other ways of importing raster data

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
summary(fm)


# Produce species distribution map (i.e., expected probability of occurrence)
beta <- coef(fm)
E.logitpsi <- beta[1] + beta[2]*elev + beta[3]*(elev^2) + beta[4]*precip
E.psi <- exp(E.logitpsi)/(1+exp(E.logitpsi))
plot(E.psi)
points(xy)



# Simulation example

set.seed(131)
x1 <- sort(rnorm(100))
x1 <- raster(outer(x1, x1), xmn=0, xmx=100, ymn=0, ymx=100)

x2 <- raster(matrix(runif(1e4), 100, 100), 0, 100, 0, 100)

# The raster package doesn't handle factors, so you need to
# code them as dummy variables.
# Note, using asFactor(x3) will not help
x3 <- raster(matrix(c(0,1), 100, 100), 0, 100, 0, 100)

logit.psi <- -1 + 1*x1 + 0*x2
psi <- exp(logit.psi)/(1+exp(logit.psi))
plot(psi)

r <- stack(x1, x2, x3)
r@layernames <- c("x1", "x2", "x3")
plot(r)

pa <- matrix(NA, 100, 100)
pa[] <- rbinom(1e4, 1, as.matrix(psi))
str(pa)
table(pa)

pa <- raster(pa, 0, 100, 0, 100)
plot(pa)

xy <- xyFromCell(pa, sample(Which(pa==1, cells=TRUE), 1000))

plot(x1)
points(xy)

fm2 <- maxlike(~x1 + x2 + x3, r, xy)

summary(fm2)
confint(fm2)
AIC(fm2)
logLik(fm2)


}

