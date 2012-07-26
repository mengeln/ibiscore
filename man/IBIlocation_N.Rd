\name{IBIlocation_N}
\alias{IBIlocation_N}
\title{
NorCal location
}
\description{
IDs NorCal ecoregion
}
\usage{
IBIlocation_N(points)
}

\arguments{
  \item{points}{
A dataframe with coordinates and StationCode.
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
    load("map2.RData")
    require("rgdal")
    coordinates(points) <- ~Longitude + Latitude
    points@proj4string <- map2@proj4string
    zone <- cbind(points@data$StationCode, over(points, map2))
    colnames(zone)[1] <- "StationCode"
    zone$LEVEL3_NAM <- as.character(zone$LEVEL3_NAM)
    zone$LEVEL3_NAM[which(is.na(zone$LEVEL3_NAM))] <- "No match"
    coast <- as.character(zone[zone$LEVEL3_NAM == "Coast Range", 
        1])
    klamath <- as.character(zone[zone$LEVEL3_NAM == "Klamath Mountains", 
        1])
    wood <- as.character(zone[zone$LEVEL3_NAM == "Southern and Central California Chaparral and Oak Woodlands", 
        1])
    outside <- as.character(zone[which(!(zone$StationCode \%in\% 
        c(coast, klamath, wood))), 1])
    location <- rep(c("Coast Range", "Klamath", "Chaparral", 
        NA), times = c(length(coast), length(klamath), length(wood), 
        length(outside)))
    names(location) <- c(coast, klamath, wood, outside)
    location <- location[match(zone$StationCode, names(location))]
    return(location)
  }
}
\keyword{ ibi }
