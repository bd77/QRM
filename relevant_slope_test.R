# -----------------------
# test relevant_slope.R
# -----------------------

wd <- "C:/Documenten/QRM/Rscripts"
wd <- "D:/QRM/Rscripts"
setwd(wd)

library(raster)
library(rgdal)
library(rgeos)
source(file = "relevant_slope.R")
source("DEMderiv.R")

# read DEM of the Alps
dem.alps <- raster("../DEM/Alps/EUD_CP-DEMS_4500025000-AA.tif")
# crop Simplon area
lon.min <- 7.967
lon.max <- 8.127
lat.min <- 46.177
lat.max <- 46.291
corner.points <- SpatialPoints(coords = matrix(c(lon.min, lon.max, lat.min, lat.max), nrow = 2),
                               proj4string = CRS("+init=epsg:4326"))
# transform from WGS84 to LAEA
extent.simplon <- extent(spTransform(corner.points, CRS("+init=epsg:3035")))
#elevation
dem.simplon <- crop(dem.alps, extent.simplon)
aspect.simplon.rad <- terrain(dem.simplon, opt = "aspect", unit = 'radians')
slope.simplon <- terrain(dem.simplon, opt = "slope", unit = 'degrees')
plot(slope.simplon>45)
plot(dem.simplon, col = terrain.colors(20))

# plan curvature of the domain
pcurv.simplon <- DEMderiv(dem.simplon, attr = "plan.curvature", method = "evans")
plot(pcurv.simplon)


# read a track
gpx.raw <- readOGR(dsn = paste0("Spitzhorli.gpx"), layer = "track_points")
gpx.3035 <- spTransform(gpx.raw, CRS("+init=epsg:3035"))
jpeg("spitzhorli.jpg", width = 2*480, height = 2*480)
track.extent <- extent(gpx.3035)
track.extent@xmin <- track.extent@xmin-1000
track.extent@ymin <- track.extent@ymin-1000
track.extent@xmax <- track.extent@xmax+1000
track.extent@ymax <- track.extent@ymax+1000
dem.track <- crop(dem.simplon, track.extent)
plot(dem.track, col = terrain.colors(20))
track.line <- SpatialLines(list(Lines(Line(gpx.3035@coords), ID = 1)))

dev.off()
i<-21
source(file = "relevant_slope.R")
for( i in 1:NROW(gpx.3035@coords)) {
  print(i)
  # convert the track point to a 1x2 matrix
  track.point <- matrix(gpx.3035@coords[i,], ncol = 2)
  # determine the relevant slope area
  rsa.points <- rsa(dem.simplon, slope.simplon, aspect.simplon.rad, pcurv.simplon, track.point)
  
  # plot the DEM, the track, the selected point, it's RSA
  jpeg(paste0("spitzhorli_", sprintf("_%03d", i), ".jpg"), width = 2*480, height = 2*480)
  plot(dem.track, col = terrain.colors(20))
  points(track.point, col = 'red', pch=19)
  lines(track.line)
  lines(rsa.points, col="blue")
  dev.off()
  #chull(as.points(rsa.points))
}

rsa.points.matrix <- matrix(ncol = 2, nrow = 0)  
for (j in (1:3)) {
  rsa.points.matrix <- rbind(rsa.points.matrix,
                             rsa.points@lines[[j]]@Lines[[1]]@coords)
}
rsa.hull <- rsa.points.matrix[chull(rsa.points.matrix),]
points(rsa.hull, col="green", pch=19)
rsa.polygon <- SpatialPolygons( list(  Polygons(list(Polygon(rsa.hull)), 1)))

rsa.slopes <- extract(slope.simplon, rsa.polygon)
rsa.slope.quant.80 <- quantile(rsa.slopes[[1]], 0.8)
rsa.slope.mean <- mean(rsa.slopes[[1]])
rsa.slope.max <- max(rsa.slopes[[1]])



plot(rsa.polygon, add=T)
polygon.raster <- crop(dem.track, extent(rsa.polygon))
r.slope.polygon <- crop(slope.simplon, extent(rsa.polygon))
points.in.polygon <- over(SpatialPoints(rasterToPoints(r.slope.polygon)), 
                          rsa.polygon)
points.in.polygon[is.na(points.in.polygon)] <- FALSE
points.in.polygon[points.in.polygon == 1] <- TRUE
sum(values(r.slope.polygon) * points.in.polygon) / sum(points.in.polygon)
max(values(r.slope.polygon) * points.in.polygon)
plot(r.slope.polygon)
