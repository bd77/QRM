# -------------------------------------------
# read GPS tracks and add spatial attributes
# -------------------------------------------

# load libraries
library(rgdal)
library(raster)
library(geosphere)

# clean up
rm(list = ls())

wd <- "D:/InterestingStuff/avalanches/"
setwd(wd)

# digital elevation models
# dem.alps <- raster("D:/InterestingStuff/slopeangle/EUD_CP-DEMS_4500025000-AA.tif") 
dem.ch <- raster("D:/InterestingStuff/avalanches/data/swissdem/switzerland_dem.tif")
slope.ch <- raster("D:/InterestingStuff/avalanches/data/swissdem/switzerland_slope.tif") 
aspect.ch <- terrain(dem.ch, opt='aspect', unit='degrees', neighbors=4)
# plot(aspect.ch)

# folder with gpx tracks from camptocamp.org
# Simplon area
# https://www.camptocamp.org/outings#act=skitouring&bbox=886188%252C5810119%252C897501%252C5825636&offset=150
# 671 done of 671
gpx.folder <- "gpxtracks/"
gpx.list <- list.files(gpx.folder, pattern = ".gpx$", full.names = FALSE)

# data frames for tour summary
tours.summary.df <- data.frame()
tours.df <- data.frame()

for (gpx.file in gpx.list) {
  # check if the gpx is a track or a route
  if (length(ogrFIDs(dsn = paste0(gpx.folder, gpx.file), layer = "track_points")) == 0) {
    print(paste0(gpx.file, " is not a track"))
  } else {
    c2c.id <- sub('.gpx', '', gpx.file)
    # read the track
    gpx.raw <- readOGR(dsn = paste0(gpx.folder, gpx.file), layer = "track_points")
    np <- NROW(gpx.raw@data)
    gpx.raw@coords
    gpx.raw@data <- cbind(c2c.id, 
                          gpx.raw@data[, c("track_seg_point_id", "ele", "time")])
    
    # calculate distances between points
    dx.m <- rep(0, times=np)
    for (i in (2:np)) {
      dx.m[i] <- distGeo(c(gpx.raw@coords[i-1, 1], gpx.raw@coords[i-1, 2]),
                         c(gpx.raw@coords[i, 1], gpx.raw@coords[i, 2])) 
    }
    # calculate tour length in km
    tour.length.km <- sum(dx.m) / 1000
    
    # process the time stamps
    datetimes <- strptime(gpx.raw@data$time, format = "%Y/%m/%d %H:%M:%S")
    time.s <- c(0, datetimes[2:np] - datetimes[1:(np-1)])
    # route time in hours
    tour.time.h <- sum(time.s) / 3600
    
    # slope in degrees
    slope.deg <- extract(slope.ch, gpx.raw@coords)

    # aspect in degrees from North, clockwise
    aspect.deg <- extract(aspect.ch, gpx.raw@coords)
    
    # calculate speed
    # speed.kmph <- 
    
    # calculate uphill/downhill/flat
    
    # add the tour summary data to the summary df
    tours.summary.df <- rbind(tours.summary.df,
                              data.frame(c2c.id = c2c.id,
                                         tour.date = strptime(datetimes[1], format = "%Y-%m-%d"),
                                         tour.length.km = tour.length.km,
                                         tour.time.h = tour.time.h,
                                         n.points = length(datetimes)))
    tour.df <- data.frame(gpx.raw@data,
                          lon = gpx.raw@coords[, 1],
                          lat = gpx.raw@coords[, 2],
                          time.s = time.s,
                          dx.m = dx.m,
                          slope.deg = slope.deg,
                          aspect.deg = aspect.deg)
    
    tours.df <- rbind(tours.df, tour.df)
  }
}  

acf(tours.df$slope.deg, lag.max = NULL,
    type = "correlation",
    plot = TRUE)

acf(tours.df$aspect.deg, lag.max = NULL,
    type = "correlation",
    plot = TRUE)

# spatial variogram
library(geoR)

slope.data <- as.geodata(tours.df[, c("lon", "lat", "slope.deg")],
                         coords.col = 1:2, data.col = 3)

slope.variog <- variog(slope.data, breaks = seq(0, 0.1, 0.005))
plot(slope.variog)
mean(tours.df$dx.m)
sd(tours.df$dx.m)

#   ogrListLayers(dsn = paste0(gpx.folder, gpx.file))
#   ogrFIDs(dsn = paste0(gpx.folder, gpx.file), layer = "track_points")
#   ogrInfo(dsn = paste0(gpx.folder, gpx.file))
#   gpx.raw <- readOGR(dsn = paste0(gpx.folder, gpx.file), layer = "route_points")
#   gpx.raw@data <- cbind(c2c.id, 
#                         gpx.raw@data[, c("route_point_id", "ele", "time")])
#   
# }
# sr.wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# gpx.ext <- as(extent(gpx.raw), 'SpatialPolygons')
# proj4string(gpx.ext) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# 
# gpx.ext.3035 <- extent(spTransform(gpx.ext, CRS("+init=epsg:3035")))
# dem.gpx <- crop(dem.alps, gpx.ext.3035)
# 
# 
# slope.gpx <- terrain(dem.gpx, opt='slope', unit='degrees', neighbors=4)
# aspect.gpx <- terrain(dem.gpx, opt='aspect', unit='degrees', neighbors=4)
# 
# gpx.3035 <- spTransform(gpx.raw, CRS("+init=epsg:3035"))
# 
# plot(slope.gpx)
# plot(spTransform(gpx.raw, CRS("+init=epsg:3035")) , add=T, col="red")
# 
# plot(aspect.gpx)
# plot(spTransform(gpx.raw, CRS("+init=epsg:3035")) , add=T, col="red")
# plot(rasterToContour(dem.gpx), levels = seq(from = 1700, to = 3000, by = 20), add=T)
# 
# gpx.slope <- extract(slope.gpx, gpx.3035@coords)
# gpx.aspect <- extract(aspect.gpx, gpx.3035@coords)