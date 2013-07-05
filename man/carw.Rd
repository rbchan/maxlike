\name{carw}
\alias{carw}
\alias{carw.data}
\docType{data}
\title{
  The Carolina Wren data used by Royle et al. (2012)
}
\description{
  These data come from the North American Breeding Bird Survey. They
  include presence-only and presence-absence data for all 25
  square-kilometer covering the contiguous United States. See
  Royle et al. (2012) for details.
}
\usage{data(carw)}
\format{
  The format is:
List of 5
$ raster.data:'data.frame':	31980 obs. of  6 variables:

$ pa.data    :'data.frame':	4615 obs. of  9 variables:

$ xy1        :'data.frame':	12082 obs. of  2 variables:

$ ext        : num [1:4] -3043880 3106120 114746 3364746

$ dim        : num [1:2] 130 246
}
\details{
  The component \code{raster.data} contains the spatially-referenced
  covariate data that can be coverted to a \code{\link{raster}}
  object. \code{pa.data} is a data.frame of presence-absence data.
  \code{xy1} is a data.frame of coordinates of routes where Carolina
  Wrens were detected. The other components are attibutes necessary for
  formatting \code{raster.data} as a \code{\link{raster}} object.
}
\source{
  \url{http://www.pwrc.usgs.gov/bbs/RawData/Choose-Method.cfm}
}
\references{
  Royle, J.A., R.B. Chandler, C. Yackulic and J.
  D. Nichols. 2012. Likelihood analysis of species occurrence probability
  from presence-only data for modelling species distributions. Methods
  in Ecology and Evolution. doi: 10.1111/j.2041-210X.2011.00182.x
}
\examples{
data(carw)

# Convert data frame to a list of rasters
rl <- lapply(carw.data$raster.data, function(x) {
   m <- matrix(x, nrow=carw.data$dim[1], ncol=carw.data$dim[2],
               byrow=TRUE)
   r <- raster(m)
   extent(r) <- carw.data$ext
   r
})

# Stack and name them
rs <- stack(rl[[1]], rl[[2]], rl[[3]], rl[[4]], rl[[5]], rl[[6]])
names(rs) <- names(carw.data$raster.data)

plot(rs)

}
\keyword{datasets}
