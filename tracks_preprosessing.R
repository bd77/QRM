# ---------------------------
# Track preprocessing
# - Filter the tracks out
# - Plot all the tracks 
# - Make summary table
#----------------------------

# clean up
rm(list = ls())

library(ggmap)

# path of the gpx tracks and file list
gpx.track.path <- "C:/Documenten/QRM/gpxtracks/"
gpx.plots.path <- "C:/Documenten/QRM/gpxtracks/plots/"
gpx.summary.path <- "C:/Documenten/QRM/gpxtracks/summary/"
track.list <- dir(path = gpx.track.path, pattern = ".gpx$")

# summary dataframe for the tracks, skip tracks already included
summary.file <- paste0(gpx.summary.path, "tracks_summary.txt")
if (file.exists(summary.file)) {
  tracks.summary.df <- read.table(summary.file, header = TRUE, sep = "\t")
                                  #as.is = c("start.time"))
  # exclude existing tracks
  track.list <- track.list[!(track.list %in% tracks.summary.df$gpx.file)]
} else {
  tracks.summary.df <- data.frame()
}

# function calculating track properties
calculate.track.props <- function(gpx.spdf, gpx.plots.path, gpx.file) {
  # it is a proper track
  np <- NROW(gpx.spdf@coords)
  # convert datetime string in a time object
  # %OS for seconds because sometimes they have 3 digits after the comma
  track.time <- strptime(paste0(gpx.spdf@data$time, "00"), "%Y/%m/%d %H:%M:%OS%z")
  # time intervals between points
  track.dt.s <- c(0, track.time[2:np] - track.time[1:(np-1)])
  # calculate track properties
  start.time <- track.time[1]
  track.date <- format(start.time, "%Y/%m/%d")
  finish.time <- track.time[np]
  # gpx.sl <- SpatialLines(list(Lines(list(Line(gpx.spdf@coords)), ID=1)), 
  #              proj4string = CRS(proj4string(gpx.spdf)))
  ele.min <- min(gpx.spdf@data$ele, na.rm = TRUE)
  ele.max <- max(gpx.spdf@data$ele, na.rm = TRUE)
  track.length.m <- geosphere::lengthLine(gpx.spdf@coords)
  track.time.h <- as.numeric(finish.time - start.time)
  avg.speed.kmph <- track.length.m / 1000 / track.time.h
  track.dh <- c(0, gpx.spdf@data$ele[2:np] - gpx.spdf@data$ele[1:(np-1)])
  uphill.m <- sum(track.dh[track.dh > 0], na.rm = TRUE)
  time.up.s <- sum(track.dt.s[track.dh > 0], na.rm = TRUE)
  speed.up.mph <- uphill.m / time.up.s * 3600
  downhill.m <- sum(track.dh[track.dh < 0], na.rm = TRUE)
  time.down.s <- sum(track.dt.s[track.dh < 0], na.rm = TRUE)
  speed.down.mph <- downhill.m / time.down.s * 3600
  
  # speed up and downhill: many points have 0 elevation gain.
  # These have to be assigned to the uphill or downhill
  
  #(speed.up.mph + time.down.s)
  # store track properties in a data frame
  track.summary.df <- data.frame(gpx.file = gpx.file,
                                 track.date = track.date,
                                 start.time = start.time,
                                 finish.time = finish.time,
                                 ele.min = ele.min,
                                 ele.max = ele.max,
                                 uphill.m = uphill.m,
                                 downhill.m = downhill.m,
                                 track.length.m = track.length.m,
                                 track.time.h = track.time.h,
                                 avg.speed.kmph = avg.speed.kmph)

  # # plot the track
  # e.track <- extent(gpx.spdf)
  # e.track.width <- e.track@xmax-e.track@xmin
  # e.track.height <- e.track@ymax-e.track@ymin
  # # enlarged extent as bottom left/top right box
  # e.ggmap <- c(e.track@xmin-e.track.width/2, e.track@ymin-e.track.height/2,
  #              e.track@xmax+e.track.width/2, e.track@ymax+e.track.height/2)
  # # get a google maps map
  # map <- get_map(location = e.ggmap)
  # m <- ggmap(map) 
  # m2 <- m + geom_path(data = as.data.frame(gpx.spdf@coords),
  #                     aes(x = coords.x1, y = coords.x2),
  #                     col = "red", lwd = 1)
  # # plot the map with the track
  # jpeg(filename = paste0(gpx.plots.path, sub(pattern = ".gpx", ".jpg", gpx.file)),
  #      width = 4*480, height = 4*480, pointsize = 2*12)
  # print(m2)
  # dev.off()
  
  return(track.summary.df)
}

# gpx.file <- "14495_16908_95538_GPS-Daten.gpx"
gpx.file <- "215850.gpx" # route that is not a track
# gpx.file <- "969507.gpx"
gpx.file <- "516053.gpx" # route that is a track
gpx.file <- "163200.gpx" # route without date
gpx.file <- "4793_18786_92471_GPS-Daten.gpx" # not read correctly
gpx.file <- "486660.gpx" # no elevation

for (gpx.file in track.list) {
  print(gpx.file)
  gpx.file.full.path <- paste0(gpx.track.path, gpx.file)
  
  # try to open the gpx file as a track
  try.track <- try(gpx.track.spdf<- readOGR(dsn = gpx.file.full.path, layer = "track_points"))
  if (!inherits(try.track, "try-error")) {
    # It reads as a track
    track.summary.df <- calculate.track.props(gpx.track.spdf, gpx.plots.path, gpx.file)  
    # write.table(file = paste0(gpx.summary.path, sub(pattern = ".gpx", ".txt", gpx.file)),
    #             track.summary.df, row.names = FALSE, sep = '\t')
    # if (is.na(track.summary.df$start.time)) {
    #   browser()
    # }
    tracks.summary.df <- rbind(tracks.summary.df, track.summary.df)
    # if (is.na(sum(is.na(tracks.summary.df$start.time)))) {
    #   browser()
    # }
    
  } else {
    # try to read as a route
    try.route <- try(gpx.route.spdf<- readOGR(dsn = gpx.file.full.path, layer = "route_points"))
    if (!inherits(try.route, "try-error")) {
      # It reads as a track
      # Check if the time field contains only NA
      pct.NA <- sum(is.na(gpx.route.spdf@data$time)) / NROW(gpx.route.spdf@data) * 100
      if (pct.NA == 100) {
        # no time field, this is a route. Copy to routes
        print(paste0("Neither route nor track: ", gpx.file, ". Moved to /routes"))
        if (file.copy(gpx.file.full.path, paste0(gpx.track.path, "routes/", gpx.file))) {
          file.remove(gpx.file.full.path) }
      } else {
        track.summary.df <- calculate.track.props(gpx.route.spdf, gpx.plots.path, gpx.file) 
        # write.table(file = paste0(gpx.summary.path, sub(pattern = ".gpx", ".txt", gpx.file)),
        #             track.summary.df, row.names = FALSE, sep = '\t')
        # if (is.na(track.summary.df$start.time)) {
        #   browser()
        # }
        tracks.summary.df <- rbind(tracks.summary.df, track.summary.df)
        # if (is.na(sum(is.na(tracks.summary.df$start.time)))) {
        #   browser()
        # }
        
      }
    } else {
      # dump it
      print(paste0("Neither route nor track: ", gpx.file, ". Moved to /undefined."))
      if (file.copy(gpx.file.full.path, paste0(gpx.track.path, "undefined/", gpx.file))) {
        file.remove(gpx.file.full.path) }
    }
  }
} # loop over all gpx files
      
    # # check if it's a route instead of a track
    #   print(paste0("Not a track but a route: ", gpx.file, ". Moved to /routes"))
    #   if (file.copy(gpx.file.full.path, paste0(gpx.track.path, "routes/", gpx.file))) {
    #     file.remove(gpx.file.full.path) }

write.table(file = paste0(gpx.summary.path, "tracks_summary.txt"), sep = "\t",
                          tracks.summary.df, row.names = FALSE, quote = FALSE)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# code to make a map with slope background. This takes time
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# map <- get_map(location = e.ggmap)
# map
# str(map)
#         #+ geom_raster(simplon.slope.wgs84)
# lines(gpx.spdf@coords)
# 
# # crop slope map to the extent of the google map of the track
# corner.points.wgs84 <- SpatialPoints(coords = t(matrix(e.ggmap, ncol = 2)),
#                                      proj4string = CRS("+init=epsg:4326"))
# corner.points.3035 <- spTransform(corner.points.wgs84, CRS("+init=epsg:3035"))
# r.slope.track <- crop(slope.simplon.deg, extent(corner.points.3035))
# plot(r.slope.track)
# r.slope.track.wgs84 <- projectRaster(r.slope.track, crs = proj4string(gpx.spdf))
# plot(r.slope.track.wgs84)
# # convert raster into data.frame with x, y and slope
# slope.df <- data.frame(rasterToPoints(r.slope.track.wgs84))
# brks <- c(seq(from = 0, to = 30, by =5), 35, 40, 45, 55, 90)
# nbrks <- length(brks)-1
# slopecols <- c(colorRampPalette(c('white', 'olivedrab1'))(6), 
#                'yellow', 'orange', 'firebrick1', 'darkorchid', 'black')
# 
# g <- ggplot() + geom_raster(data = slope.df, aes(x=x, y=y, fill = slope, alpha = 0.3), interpolate = FALSE) + scale_fill_gradientn(colours = slopecols, breaks = brks)
# g
# 
# m3 <- m2 + geom_tile(data = slope.df, aes(x=x, y=y, fill = slope, alpha = 0.1)) + scale_fill_gradientn(colours = slopecols, breaks = brks)
# m3
# 
# spTransform(rasterToPoints(), proj4string(gpx.spdf))
# simplon.slope.wgs84 <- projectRaster(slope.simplon.deg, 
#                                      crs = CRS(projargs = proj4string(gpx.spdf)))
# 
# g <- ggplot(faithfuld, aes(waiting, eruptions)) + geom_raster(aes(fill = density))
# g
# g <- ggplot() + geom_raster(data = faithfuld,
#                                    aes(x=waiting, y=eruptions, fill = density))
# g
