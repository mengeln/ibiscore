\name{IBIlocation}
\alias{IBIlocation}
\title{
SoCal location
}
\description{
Ids ecoregion
}
\usage{
IBIlocation(points)
}

\arguments{
  \item{points}{
A dataframe of coordinates with StationCode.
}
}
\examples{
## The function is currently defined as
function (points) 
{
    if (is.data.frame(points) == F) {
        print("Input must be a data frame")
        stop
    }
    if (colnames(points)[1] != "StationCode") {
        print("First column must be 'StationCode'")
        stop
    }
    load("map.RData")
    require("rgdal")
    coordinates(points) <- ~Longitude + Latitude
    points@proj4string <- map@proj4string
    zone <- cbind(points@data$StationCode, over(points, map))
    colnames(zone)[1] <- "StationCode"
    zone$LEVEL3_NAM <- as.character(zone$LEVEL3_NAM)
    zone$LEVEL3_NAM[which(is.na(zone$LEVEL3_NAM))] <- "No match"
    chaparral <- as.character(zone[zone$LEVEL3_NAM == "Southern and Central California Chaparral and Oak Woodlands", 
        1])
    mountains <- as.character(zone[zone$LEVEL3_NAM == "Southern California Mountains", 
        1])
    outside <- as.character(zone[which(!(zone$StationCode \%in\% 
        c(chaparral, mountains))), 1])
    location <- rep(c("Chaparral", "Mountains", NA), times = c(length(chaparral), 
        length(mountains), length(outside)))
    names(location) <- c(chaparral, mountains, outside)
    location <- location[match(zone$StationCode, names(location))]
    return(location)
  }
}

\keyword{ ibi }
