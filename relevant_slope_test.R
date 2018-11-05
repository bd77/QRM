# -----------------------
#  test relevant_slope.R
# -----------------------

# clean up
rm(list = ls())

wd <- "C:/Documenten/QRM/Rscripts"
# wd <- "D:/QRM/Rscripts"
setwd(wd)

library(raster)
library(rgdal)
library(rgeos)
source(file = "relevant_slope.R")
source(file = "DEMderiv.R")

# read DEM of the Alps in LAEA epsg:3035
dem.alps <- raster("../DEM/Alps/EUD_CP-DEMS_4500025000-AA.tif")
# slope.alps <- terrain(dem.alps, opt = "slope", unit = 'degrees')
# writeRaster(slope.alps, "../DEM/Alps/slope_alps.tif")

# Extent of the Simplon area in WGS84, epsg:4326
lon.min <- 7.966
lon.max <- 8.603
lat.min <- 46.160
lat.max <- 46.543
corner.points <- SpatialPoints(coords = matrix(c(lon.min, lon.max, lat.min, lat.max), nrow = 2),
                               proj4string = CRS("+init=epsg:4326"))
# transform Simplon extent from WGS84 to LAEA
extent.simplon <- extent(spTransform(corner.points, CRS("+init=epsg:3035")))

# crop the DEM to the Simplon region
dem.simplon <- crop(dem.alps, extent.simplon)
rm(dem.alps) # save some memory

# -- calculate slope attributes --
# aspect in radians
aspect.simplon.rad <- terrain(dem.simplon, opt = "aspect", unit = 'radians')
# slope angle in degrees
slope.simplon.deg <- terrain(dem.simplon, opt = "slope", unit = 'degrees')
brks <- c(seq(from = 0, to = 30, by =5), 35, 40, 45, 55, 90)
nbrks <- length(brks)-1
slopecols <- c(colorRampPalette(c('white', 'olivedrab1'))(6), 
               'yellow', 'orange', 'firebrick1', 'darkorchid', 'black')
png('Slopes_Simlon_deg.png', height = 4*480, width = 4*480,
     pointsize = 2*12)
plot(slope.simplon.deg, breaks=brks, col=slopecols, asp=1)
dev.off()
# cliffs
plot(slope.simplon.deg > 45)
plot(dem.simplon, col = terrain.colors(20))

# plan curvature of the domain
pcurv.simplon <- DEMderiv(dem.simplon, attr = "plan.curvature", method = "evans")
# plot(pcurv.simplon)

# loop over all tracks
gpx.track.path <- "C:/Documenten/QRM/gpxtracks/"
animations.path <- "C:/Documenten/QRM/Animations/"
track.list <- dir(path = gpx.track.path, pattern = "Daten.gpx$")
gpx.file <- "14495_16908_95538_GPS-Daten.gpx"

track.df <- data.frame()

for (gpx.file in track.list) {
  # read a track
  gpx.file.full.path <- paste0(gpx.track.path, gpx.file)
  gpx.name.no.ext <- sub(".gpx$", "", gpx.file)
  gpx.raw <- readOGR(dsn = gpx.file.full.path, layer = "track_points")
  
  # folder for the animation frames
  track.animation.path <- paste0(animations.path, gpx.name.no.ext)
  if (!dir.exists(track.animation.path)) {
    dir.create(track.animation.path)
  }
  
  # transform it to LAEA, the projection of the DEM
  gpx.track <- spTransform(gpx.raw, proj4string(dem.simplon))
  track.extent <- extent(gBuffer(as(extent(gpx.track), 'SpatialPolygons'), width = 1000))
  dem.track <- crop(dem.simplon, track.extent)
  
  

  i<-21
  # source(file = "relevant_slope.R")
  # loop over the whole track
  np <- NROW(gpx.track@coords)
  for(i in 1:np) {
    print(paste0(round(i/np*100), "%"))
    # convert the track point to a 1x2 matrix
    track.point <- matrix(gpx.track@coords[i,], ncol = 2)
    # determine the relevant slope area
    # source(file = "relevant_slope.R")
    rsa.props <- rsa(dem.simplon, slope.simplon.deg, aspect.simplon.rad, pcurv.simplon, track.point)
    
    track.point.df <- data.frame(gpx.file = gpx.name.no.ext,
                                 point.id = i,
                                 x = track.point[1,1],
                                 y = track.point[1,1],
                                 lon = as.numeric(gpx.raw@coords[i,1]),
                                 lat = as.numeric(gpx.raw@coords[i,2]),
                                 local.slope = rsa.props@localSlope,
                                 local.aspect = rsa.props@localAspect,
                                 max.slope = rsa.props@maxSlope@data$slope)
    track.df <- rbind(track.df, track.point.df)
  
    # plot the DEM, the track, the selected point, it's RSA
    # jpeg(paste0(track.animation.path, "/", sprintf("_%04d", i), ".jpg"), 
    #      width = 4*480, height = 2*480,
    #      pointsize = 24)
    # par(mfrow=c(1,2))
    # plot(dem.track, col = terrain.colors(20), asp=1)
    # points(track.point, col = 'red', pch=19)
    # track.line <- Line(gpx.track@coords)
    # lines(track.line)
    # lines(rsa.props@borderLines, col="blue")
    # 
    # # detailed slope plot
    # local.extent <- extent(gBuffer(rsa.props@relevantSlopePolygon, width = 100))
    # local.slope.deg <- crop(slope.simplon.deg, local.extent)
    # plot(local.slope.deg, breaks=brks, col=slopecols, asp=1 )
    # points(track.point, col = 'red', pch=19)
    # lines(track.line)
    # lines(rsa.props@borderLines, col="blue")
    # plot(rsa.props@relevantSlopePolygon, add=T)
    # points(rsa.props@maxSlope)
    # text(x = rsa.props@maxSlope@coords[1], y = rsa.props@maxSlope@coords[2],
    #      paste0("Max slope ", round(rsa.props@maxSlope@data$slope)), pos=4)
    # dev.off()
    
  } # loop over the track
}

library(geoR)
v1 <- variog(coords = as.matrix(track.df[, c("x", "y")]), data = track.df$local.slope,
             breaks = seq(0, 3000, 50))
# breaks = c(0, 2^(0:12))log="x",
plot(v1)
v1 <- variog(coords = as.matrix(track.df[, c("x", "y")]), data = track.df$max.slope,
             breaks = seq(0, 3000, 50))
# breaks = c(0, 2^(0:12))log="x",
plot(v1)
var(track.df$max.slope)

library(gstat)
coordinates(track.df) = ~x+y
local.slope.variogram <- variogram(local.slope~1, track.df, width = 10)
plot(local.slope.variogram, xlim = c(0, 100))

local.aspect.variogram <- variogram(local.aspect~1, track.df, width = 10, cutoff = 10000)
plot(local.aspect.variogram)


# rsa.points.matrix <- matrix(ncol = 2, nrow = 0)  
# for (j in (1:3)) {
#   rsa.points.matrix <- rbind(rsa.points.matrix,
#                              rsa.points@lines[[j]]@Lines[[1]]@coords)
# }
# rsa.hull <- rsa.points.matrix[chull(rsa.points.matrix),]
# points(rsa.hull, col="green", pch=19)
# rsa.polygon <- SpatialPolygons( list(  Polygons(list(Polygon(rsa.hull)), 1)))
# 
# rsa.slopes <- extract(slope.simplon, rsa.polygon)
# rsa.slope.quant.80 <- quantile(rsa.slopes[[1]], 0.8)
# rsa.slope.mean <- mean(rsa.slopes[[1]])
# rsa.slope.max <- max(rsa.slopes[[1]])
# 
# rsa.pcurv <- extract(pcurv.simplon, rsa.polygon)
# rsa.pcurv.mean <- mean(rsa.pcurv[[1]])
# 
# extract(slope.simplon, rsa.polygon, cellnumbers=T)
# 
# rasterToPoints(slope.simplon)[75656,]
# values(slope.simplon)[75656]
# extract(slope.simplon, matrix(c(4166663, 2573388), ncol = 2)) coordinates(slope.simplon)[75656,])
# 
# plot(rsa.polygon, add=T)
# polygon.raster <- crop(dem.track, extent(rsa.polygon))
# r.slope.polygon <- crop(slope.simplon, extent(rsa.polygon))
# points.in.polygon <- over(SpatialPoints(rasterToPoints(r.slope.polygon)), 
#                           rsa.polygon)
# points.in.polygon[is.na(points.in.polygon)] <- FALSE
# points.in.polygon[points.in.polygon == 1] <- TRUE
# sum(values(r.slope.polygon) * points.in.polygon) / sum(points.in.polygon)
# max(values(r.slope.polygon) * points.in.polygon)
# plot(r.slope.polygon)
