\name{MaungaWhau}
\alias{MaungaWhau}
\docType{data}
\title{
Fake data for the Maunga Whau volcano
}
\description{
Based upon \link{volcano} in the datasets package.
}
\usage{data(MaungaWhau)}
\format{
  The format is:
List of 3
 $ elev  : num [1:87, 1:61] -1.17 -1.13 -1.09 -1.05 -1.01 ...
 $ precip: num [1:87, 1:61] -2.27 -2.25 -2.24 -2.22 -2.2 ...
 $ xy    : num [1:1000, 1:2] 38.5 42.5 30.5 25.5 16.5 36.5 22.5 56.5 49.5 18.5 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "x" "y"
}
\details{
No coordinate system attributes are included
}
\examples{
data(MaungaWhau)
elev <- raster(MaungaWhau$elev, xmn=0, xmx=61, ymn=0, ymx=87)
precip <- raster(MaungaWhau$precip, xmn=0, xmx=61, ymn=0, ymx=87)
}
\keyword{datasets}
