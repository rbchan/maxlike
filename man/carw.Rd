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
List of 6
 $ raster.data:'data.frame':	31980 obs. of  6 variables:
  ..$ pcMix: num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
  ..$ pcDec: num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
  ..$ pcCon: num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
  ..$ pcGr : num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
  ..$ Lat  : num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
  ..$ Lon  : num [1:31980] NA NA NA NA NA NA NA NA NA NA ...
 $ pa.data    :'data.frame':	4615 obs. of  9 variables:
  ..$ pcMix: num [1:4615] 0.577 0.198 1.473 0.847 0.708 ...
  ..$ pcDec: num [1:4615] 1.399 1.33 2.897 2.031 0.836 ...
  ..$ pcCon: num [1:4615] -0.291 -0.249 0.177 -0.168 -0.033 ...
  ..$ pcGr : num [1:4615] -0.356 -0.344 -0.348 -0.337 -0.303 ...
  ..$ Lat  : num [1:4615] -0.413 -0.413 -0.413 -0.466 -0.44 ...
  ..$ Lon  : num [1:4615] 0.415 0.401 0.387 0.458 0.458 ...
  ..$ y    : num [1:4615] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ X    : num [1:4615] 768620 743620 718620 843620 843620 ...
  ..$ Y    : num [1:4615] 1352246 1352246 1352246 1302246 1327246 ...
 $ xy1        :'data.frame':	12082 obs. of  2 variables:
  ..$ X: num [1:12082] 760241 757836 755943 755301 754561 ...
  ..$ Y: num [1:12082] 1346223 1344781 1343874 1344089 1343971 ...
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

# Covert data frame to a list of rasters
rl <- lapply(carw.data$raster.data, function(x) {
   m <- matrix(x, nrow=carw.data$dim[1], ncol=carw.data$dim[2],
               byrow=TRUE)
   r <- raster(m)
   extent(r) <- carw.data$ext
   r
})

# Stack and name them
rs <- stack(rl[[1]], rl[[2]], rl[[3]], rl[[4]], rl[[5]], rl[[6]])
layerNames(rs) <- names(carw.data$raster.data)

plot(rs)

}
\keyword{datasets}
