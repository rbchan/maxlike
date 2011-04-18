\name{MaungaWhau}
\alias{MaungaWhau}
\docType{data}
\title{
Fake data on elevation and precipitation for the Manuga Whau volcano
}
\description{
This is a list with 2 \link{raster} objects and a matrix of x and y
coordinate values for the presence locations.
}
\usage{data(MaungaWhau)}
\format{
  The format is:
List of 3
 $ elev  :Formal class 'RasterLayer' [package "raster"] with 13 slots
 $ precip:Formal class 'RasterLayer' [package "raster"] with 13 slots
 $ xy    : num [1:1000, 1:2] 38.5 42.5 30.5 25.5 16.5 36.5 22.5 56.5 49.5 18.5 ...
}
\details{
Both elev and precip are standardized to have mean 0 and unit variance.
}
\seealso{\link{volcano}}
\examples{
data(MaungaWhau)
plot(MaungaWhau$elev)
plot(ManugaWhau$precip)
}
\keyword{datasets}
