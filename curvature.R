

setwd("D:/InterestingStuff/avalanches/R")
library(raster)
source("DEMderiv.R")

# read DEM of the Alps
dem.alps <- raster("D:/InterestingStuff/slopeangle/EUD_CP-DEMS_4500025000-AA.tif")

# crop an area
lon.min <- 7.967
lon.max <- 8.127
lat.min <- 46.177
lat.max <- 46.291
corner.points <- SpatialPoints(coords = matrix(c(lon.min, lon.max, lat.min, lat.max), nrow = 2),
                               proj4string = CRS("+init=epsg:4326"))
spTransform(corner.points, CRS("+init=epsg:3035"))
extent.simplon <- extent(spTransform(corner.points, CRS("+init=epsg:3035")))
dem.simplon <- crop(dem.alps, extent.simplon)
aspect.simplon.rad <- terrain(dem.simplon, opt = "aspect", unit = 'radians')
slope.simplon <- terrain(dem.simplon, opt = "slope", unit = 'degrees')
plot(slope.simplon>45)
plot(dem.simplon, col = terrain.colors(20))

pcurv <- DEMderiv(dem.simplon, attr = "plan.curvature", method = "evans")
plot(pcurv)

prof.curv <- DEMderiv(dem.simplon, attr = "prof.curvature", method = "evans")
plot(prof.curv > 0.0035)

gpx.raw <- readOGR(dsn = paste0("Spitzhorli.gpx"), layer = "track_points")
gpx.3035 <- spTransform(gpx.raw, CRS("+init=epsg:3035"))

plot(dem.simplon, col = terrain.colors(20))
points(gpx.3035)
# steepest way up from a point 
track.pt.i <- 10
uphill.line.coords <- matrix(0, nrow = 50, ncol = 2)
#<- Line(matrix(0, ncol = 2, nrow = 10))
for (i.trackpoint in 1:NROW(gpx.3035@coords)) {
  uphill.line.coords <- matrix(0, nrow = 50, ncol = 2)
  uphill.line.coords[1,] <- gpx.3035@coords[i.trackpoint,]
  step <- 25
  for (i in 2:50) {
    previous.pt <- uphill.line.coords[(i-1),]
    aspect.prev <- extract(aspect.simplon.rad, matrix(previous.pt, ncol = 2))
    slope.prev <-
    next.pt <- previous.pt + step * -c(sin(aspect.prev), cos(aspect.prev))
    uphill.line.coords[i,] <- next.pt
  }
  
  uphill.line <- SpatialLines(list(Lines(list(Line(uphill.line.coords)), ID = "a")), CRS("+init=epsg:3035"))
  plot(uphill.line, add=TRUE, col = "red", lwd = 3)
  
}


extract(pcurv, matrix(uphill.line.coords[1,], ncol = 2))

i<-2
