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

# loop over all tracks
gpx.track.path <- "../gpxtracks/"
animations.path <- "../Animations/"
# track.list <- dir(path = gpx.track.path, pattern = ".gpx$")

# read the tracks from the points data frame
track.points.df <- read.table(file = "../gpxtracks/summary/tracks_points.txt",
                              header = TRUE, sep = "\t")
track.list <- unique(track.points.df$gpx.file)
e.tracks <- extent(c(min(track.points.df$lon)-0.1, max(track.points.df$lon)+0.1,
                     min(track.points.df$lat)-0.1, max(track.points.df$lat)+0.1))

# read DEM of the Alps in LAEA epsg:3035
dem.alps <- raster("../DEM/Alps/EUD_CP-DEMS_4500025000-AA.tif")
# slope.alps <- terrain(dem.alps, opt = "slope", unit = 'degrees')
# writeRaster(slope.alps, "../DEM/Alps/slope_alps.tif")

# Extent of the tracks area in WGS84, epsg:4326 as a polygon
corner.points <- SpatialPoints(coords = matrix(c(e.tracks@xmin, e.tracks@xmax, e.tracks@ymin, e.tracks@ymax), nrow = 2),
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

if (file.exists("track_point_properties.txt")) {
  tpp.df <- read.table("track_point_properties.txt", sep = "\t", header = TRUE)
  track.list <- track.list[!(track.list %in% paste0(unique(tpp.df$gpx.file), ".gpx"))]
} else {
  track.df <- data.frame()
}

gpx.file <- "2018-01-28 Canalone Gelato.gpx"

for (gpx.file in track.list) {
  # read a track
  # gpx.file.full.path <- paste0(gpx.track.path, gpx.file)
  gpx.name.no.ext <- sub(".gpx$", "", gpx.file)
  # gpx.raw <- readOGR(dsn = gpx.file.full.path, layer = "track_points")
  
  # get the points data from the points table (reading gpx is slow)
  mytrack.points.df <- track.points.df[track.points.df$gpx.file == gpx.file,]
  gpx.raw <- SpatialPointsDataFrame(coords = mytrack.points.df[,c("lon", "lat")],
                                    data = mytrack.points.df[, c("dt.s", "dx.m", "ele")],
                                    proj4string = CRS("+init=epsg:4326"))
  
  # folder for the animation frames
  track.animation.path <- paste0(animations.path, gpx.name.no.ext)
  if (!dir.exists(track.animation.path)) {
    dir.create(track.animation.path)
  }
  
  # transform it to LAEA, the projection of the DEM
  gpx.track <- spTransform(gpx.raw, proj4string(dem.simplon))
  track.extent <- extent(gBuffer(as(extent(gpx.track), 'SpatialPolygons'), width = 1000))
  dem.track <- crop(dem.simplon, track.extent)

  i<-1400
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
                                 y = track.point[1,2],
                                 lon = as.numeric(gpx.raw@coords[i,1]),
                                 lat = as.numeric(gpx.raw@coords[i,2]),
                                 ele = extract(dem.track, track.point),
                                 local.slope = rsa.props@localSlope,
                                 local.aspect = rsa.props@localAspect,
                                 max.slope = rsa.props@maxSlope@data$slope)
    tpp.df <- rbind(tpp.df, track.point.df)
  
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

write.table(track.df, file = "track_point_properties.txt",
            row.names = FALSE, sep = "\t")

library(geoR)
v1 <- variog(coords = as.matrix(tpp.df[, c("x", "y")]), data = track.df$local.slope,
             breaks = seq(0, 1000, 20))
# breaks = c(0, 2^(0:12))log="x",
png("Variogram local slope.png")
plot(v1, main = "Variogram for local slope")
dev.off()

v1 <- variog(coords = as.matrix(tpp.df[, c("x", "y")]), data = track.df$max.slope,
             breaks = seq(0, 2000, 20))
# breaks = c(0, 2^(0:12))log="x",
png("Variogram max slope.png")
plot(v1, main = "Variogram for maximum slope")
dev.off()

var(track.df$max.slope)

library(gstat)
coordinates(tpp.df) = ~x+y
local.slope.variogram <- variogram(local.slope~1, tpp.df, width = 20)
plot(local.slope.variogram, xlim = c(0, 1000))

local.aspect.variogram <- variogram(local.aspect~1, tpp.df, width = 10, cutoff = 2000)
png("Variogram aspect.png")
plot(local.aspect.variogram, main = "Variogram for aspect")
dev.off()

# check for colinearity
png("Collinearity local and max slope.png")
plot(tpp.df$local.slope, tpp.df$max.slope,
     xlab = "Local slope (deg)", ylab = "Maximum relevant slope (deg)")
dev.off()
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
