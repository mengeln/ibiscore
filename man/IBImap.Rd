\name{IBImap}
\alias{IBImap}
\title{
Build an IBI map
}
\description{
Produces a map based on IBI scores or ecoregion data.
}
\usage{
IBImap(locationinfo, data, attribute, type = "Ecoregion", zone)
}
\arguments{
  \item{locationinfo}{
From function IBIlocation.
}
  \item{data}{
A data frame.
}
  \item{attribute}{
Attribute data to be mapped.
}
  \item{type}{
Ecoregion map or IBI map
}
  \item{zone}{
SoCal or NorCal
}
}

\examples{
## The function is currently defined as
function (locationinfo, data, attribute, type = "Ecoregion", 
    zone) 
{
    library(ggmap)
    load("map.RData")
    datalocation <- unique(data[, c("StationCode", "SampleID")])
    datalocation$attribute <- sapply(1:length(datalocation$StationCode), 
        function(i) attribute[which(names(attribute) == datalocation$StationCode[i])])
    datalocation$StationCode <- as.character(datalocation$StationCode)
    datalocation$lon <- sapply(1:length(datalocation$StationCode), 
        function(i) locationinfo[which(locationinfo$StationCode == 
            datalocation$StationCode[i]), "Longitude"])
    datalocation$lat <- sapply(1:length(datalocation$StationCode), 
        function(i) locationinfo[which(locationinfo$StationCode == 
            datalocation$StationCode[i]), "Latitude"])
    datalocation <- datalocation[which(datalocation$attribute != 
        "character(0)"), ]
    datalocation <- datalocation[which(datalocation$attribute != 
        "numeric(0)"), ]
    datalocation$attribute <- unlist(datalocation$attribute)
    datalocation$lon <- unlist(datalocation$lon)
    datalocation$lat <- unlist(datalocation$lat)
    if (zone == "SoCal") {
        testmap <- get_googlemap(center = c(lon = -118, lat = 33.7), 
            maptype = "satellite", crop = T, style = "feature:road|element:all|visibility:off", 
            zoom = 8)
    }
    if (zone == "NorCal") {
        testmap <- get_googlemap(center = c(lon = -122.454285, 
            lat = 40.501269), maptype = "satellite", crop = T, 
            style = "feature:road|element:all|visibility:off", 
            zoom = 7)
    }
    if (type == "Ecoregion") {
        ibimap <- ggmap(testmap) + geom_point(aes(x = lon, y = lat, 
            colour = attribute), data = datalocation) + labs(colour = attribute) + 
            scale_colour_discrete(name = deparse(substitute(attribute)))
    }
    else {
        datalocation$Status <- sapply(1:length(datalocation$StationCode), 
            function(i) if (datalocation$attribute[i] >= 40) {
                "Acceptable"
            }
            else {
                "Below standards"
            })
        datalocation$Status <- as.factor(datalocation$Status)
        ibimap <- ggmap(testmap) + geom_point(aes(x = lon, y = lat, 
            shape = Status, colour = attribute), size = 3, data = datalocation) + 
            labs(colour = attribute) + scale_colour_gradient(limits = c(0, 
            100), name = type, space = "Lab", low = "yellow", 
            high = "blue")
    }
  }
}

\keyword{ ibi }
