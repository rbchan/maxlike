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
  maxlike(formula, rasters, points, starts, hessian = TRUE,
          fixed, removeDuplicates=FALSE, na.action = "na.omit", ...)
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
  \item{fixed}{
    Optional vector for fixing parameters. It must be
    of length equal to the number of parameters in the
    model. If an element of fixed is NA, then the parameter is
    estimated, otherwise if it is a real number, the parameter is fixed
    at this value.
    }
  \item{removeDuplicates}{
    Logical. Should duplicate points be removed. Defaults to FALSE, but
    note that the MAXENT default is TRUE. It is unclear why duplicate
    points should be removed, but this option is provided for
    comparative purposes.}
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
  A list with 8 components
  \item{Est}{data.frame containing the parameter estimates (Ests) and
    standard errors (SE).}
  \item{vcov}{variance-covariance matrix}
  \item{AIC}{AIC}
  \item{call}{the original call}
  \item{pts.removed}{The points removed due to missing values}
  \item{pix.removed}{The pixels removed due to missing values}
  \item{optim}{The object returned by \code{\link{optim}}}
  \item{not.fixed}{A logical vector indicating if a parameter was
    estimated or fixed.}
  }
\note{
  In general it is very hard to obtain a random sample of presence
  points, which is a requirement of both the Royle et al. (2012)
  method and of MAXENT. This is one of many reasons why presence-absence
  data are preferable. When presence-absence data are available, they
  can be modeled using functions such as
  \code{\link[stats]{glm}}. Creating
  maps of prediction from \code{\link[stats]{glm}} can be accomplished
  in the
  same way as is shown in the example below.

  The MAXENT software assumes that species prevalence is known a
  priori. Even if the user does not specify a value for prevalence, a
  value is chosen automatically. In the unlikely scenario that
  prevalence is known a priori, one could use this infomation
  by fixing the
  intercept using the \code{fixed} argument in \code{maxlike}. See
  examples below.

  Another weakness of models for presence-only data is that they do not
  allow one to model detection probability, which is typically less than
  one in field conditions. If detection probability is affected by the same
  covariates that affect occurrence probability, then bias
  is inevitable. The R package \pkg{unmarked}
  (Fiske and Chandler 2011) offers numerous
  methods for jointly modeling both occurrence and detection probability
  when detection/non-detection data are available.
}

\section{Warning}{
  Maximizing the log-likelihood function is achieved using the
  \code{\link{optim}} function, which can fail to find the global optima
  of likelihood function if sensible starting values are not
  supplied. The default starting values are \code{rep(0, npars)}, which will
  often be adequate if the covariates have been standardized.
  }

\references{
  Royle, J.A., R.B. Chandler, C. Yackulic and J.
  D. Nichols. 2012. Likelihood analysis of species occurrence probability
  from presence-only data for modelling species distributions. Methods
  in Ecology and Evolution. doi: 10.1111/j.2041-210X.2011.00182.x

  Fiske, I. and R.B. Chandler. 2011. unmarked: An R Package for Fitting
  Hierarchical Models of Wildlife Occurrence and Abundance. Journal of
  Statistical Software 43(10).
}

\seealso{ \code{\link{maxlike-package}}, \code{\link{raster}},
          \code{\link{carw}} }

\examples{

# Carolina Wren data used in Royle et. al (2012)
data(carw)

# Covert data.frame to a list of rasters
rl <- lapply(carw.data$raster.data, function(x) {
   m <- matrix(x, nrow=carw.data$dim[1], ncol=carw.data$dim[2], byrow=TRUE)
   r <- raster(m)
   extent(r) <- carw.data$ext
   r
})

# Create a raster stack and add layer names
rs <- stack(rl[[1]], rl[[2]], rl[[3]], rl[[4]], rl[[5]], rl[[6]])
layerNames(rs) <- names(carw.data$raster.data)

plot(rs)


# Fit a model
fm <- maxlike(~pcMix + I(pcMix^2) + pcDec + I(pcDec^2)+ pcCon +
    I(pcCon^2) + pcGr + I(pcGr^2) +
    Lat + I(Lat^2) + Lon + I(Lon^2), rs, carw.data$xy1,
    method="BFGS", removeDuplicates=TRUE)

fm


\dontrun{
# Fit model using MAXENT.
# maxent appears to ignore the defaultprevalence argument
# To see it's effect, use the stand-along MAXENT program
library(dismo)
fm.me <- maxent(rs, carw.data$xy1,
    args=c("noproduct", "nothreshold", "nohinge", "noautofeature",
           "defaultprevalence=0.5"))
}

# MAXENT sets "default prevalence" to an arbitrary value, 0.5
# We could do the same by fixing the intercept at logit(0.5)=0
# However, it seems more appropriate to estimate this important parameter

fm.fix <- update(fm, fixed=c(0, rep(NA,length(coef(fm))-1)))



# Produce species distribution map (ie, expected probability of occurrence)
psi.hat <- predict(fm)
plot(psi.hat)
points(carw.data$xy1, pch=16, cex=0.1)






\dontrun{

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

}

