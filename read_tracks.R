

library(rgdal)
library(XML)
library(geosphere)

# https://www.camptocamp.org/outings#act=skitouring&bbox=886188%252C5810119%252C897501%252C5825636&offset=150
# 671 done of 671

track.path <- "D:/InterestingStuff/avalanches/tracks simplon/"
gpx.list <- dir(track.path, pattern = ".gpx$", full.names = FALSE)

gpx.file <- gpx.list[1]
for (gpx.file in gpx.list) {
  # read the xml file
  track <- xmlTreeParse(paste0(track.path, gpx.file), useInternalNodes = T)
  # extract data
  # track name
  track.name <- xpathSApply(track, path = "/gpx//version", xmlValue)
  elevations <- as.numeric(xpathSApply(track, path = "//rtept/ele", xmlValue))

  coords <- xpathSApply(track, path = "//rtept", xmlAttrs)
  lats <- as.numeric(coords['lat',])
  lons <- as.numeric(coords['lon',])
  
  times.str <- xpathSApply(track, path = "//rtept/time", xmlValue)
  
  times <- strptime(times.str, format = "%Y-%m-%dT%H:%M:%SZ", tz = "CET")
  nt <- length(times)
  
  dt.sec <- c(0, (times[2:nt] - times[1:(nt-1)]))
  
  dx.m <- rep(0, times=nt)
  for (i in (2:nt)) {
    dx.m[i] <- distGeo(c(lons[i-1], lats[i-1]), c(lons[i], lats[i])) 
  }
  
  speed.mps <- dx.m / dt.sec
  speed.kmph <- speed.mps * 3.6
  plot(speed.kmph)
  
  plot(lons, lats)
  track <- readOGR(dsn = gpx.file, layer = "tracks")
}
