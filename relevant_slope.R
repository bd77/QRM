# ---------------------------------------------------
# function to determine the RSA (relevant slope area)
# ---------------------------------------------------

slopeProps <- setClass("slopeProps",
                       slots = c(borderLines = "SpatialLines",
                                 relevantSlopePolygon = "SpatialPolygons",
                                 localSlope = "numeric",
                                 localAspect = "numeric",
                                 maxSlope = "SpatialPointsDataFrame",
                                 meanSlope = "numeric"))

# auxiliary functions
run.away <- function(start.point, r.aspect, step, max.dist, step.dir) {
  # input:
  # start.point: point from which a step is made
  # r.aspect: raster of aspect angle in radians clockwise from North
  # step: size of the step in meters
  # max.dist: maximum distance from start.point in meters
  # output: a Line object of points running from the start.point towards the specified
  # direction
  
  # define the correct function to convert the aspect in a direction vector
  if (step.dir == "up") { 
    aspect2vector <- function(aspect) { c(-sin(aspect), -cos(aspect)) }
  } else if (step.dir == "down") { 
    aspect2vector <- function(aspect) { c(sin(aspect), cos(aspect)) }
  } else if (step.dir == "right") { 
    aspect2vector <- function(aspect) { c(cos(aspect), -sin(aspect)) }  
  } else if (step.dir == "left") { 
    aspect2vector <- function(aspect) { c(-cos(aspect), sin(aspect)) }
  } else { 
    print("Invalid direction for run.away function")
  }
  
  # determine the point step meters from the input point in the direction
  dist2start <- 0
  run.away.line <- start.point
  current.point <- start.point
  # limit the number of steps because the line can get stuck on a ridge
  # or sink closer then max.dist
  n.iter <- 0
  max.iter <- 2*max.dist/step
  # create a line of points running in the desired direction
  while (dist2start <= max.dist & n.iter < max.iter) {
    # retrieve the aspect of the current position
    aspect.cp <- extract(r.aspect, current.point)
    # step to the next point and add it to the line
    next.point <- current.point + step * aspect2vector(aspect.cp)
    run.away.line <- rbind(run.away.line, next.point)
    # recalculate the distance to the starting point
    dist2start <- sqrt(sum((start.point - current.point)**2))
    # preparing the next step
    current.point <- next.point
    n.iter <- n.iter+1
  }
  return(Line(run.away.line))
}

rsa <- function(r.elev, r.slope, r.aspect, r.pcurv, track.point) {
  # r.elev: elevation raster
  # r.slope: raster with the slope angle in degrees
  # r.aspect: raster with the aspect in radians, N=0, E=pi/2, S=pi, W=3*pi/2
  # r.pcurv: plan curvature, <0 = concave (valley), >0 = convex (ridge)
  # track.point: 1x2 matrix with coordinates
  # the coordinate system should be orthogonal
  
  # use convex hull function "chull" to make an area out of the points
  
  # this should be in a config file
  max.dist.up <- 200
  max.dist.down <- 50
  max.dist.side <- 50
  step <- 12.5
  
  # Line object with steepest way and down up from the track point 
  uphill.line.from.track.point <- run.away(track.point, r.aspect, step, max.dist.up, "up")
  downhill.line.from.track.point <- run.away(track.point, r.aspect, step, max.dist.down, "down")
  
  # make a step left and right
  left.line.from.track.point <- run.away(track.point, r.aspect, step, max.dist.side, "left")
  npl <- NROW(left.line.from.track.point@coords)
  point.left <- matrix(left.line.from.track.point@coords[npl,], ncol = 2)
  right.line.from.track.point <- run.away(track.point, r.aspect, step, max.dist.side, "right")
  npr <- NROW(right.line.from.track.point@coords)
  point.right <- matrix(right.line.from.track.point@coords[npr,], ncol = 2)
  
  # steepest way up from the right and left point 
  uphill.line.from.left <- run.away(point.left, r.aspect, step, max.dist.up, "up")
  uphill.line.from.right <- run.away(point.right, r.aspect, step, max.dist.up, "up")

  # steepest way down from the right and left points 
  downhill.line.from.left <- run.away(point.left, r.aspect, step, max.dist.down, "down")
  downhill.line.from.right <- run.away(point.right, r.aspect, step, max.dist.down, "down")
  
  # put everything in a spatialLines object
  borderLines <- SpatialLines(list(Lines(uphill.line.from.track.point, ID="centre_up"),
                                 Lines(downhill.line.from.track.point, ID="centre_down"),
                                 Lines(uphill.line.from.left, ID="left_up"),
                                 Lines(downhill.line.from.left, ID="left_down"),
                                 Lines(uphill.line.from.right, ID="right_up"),
                                 Lines(downhill.line.from.right, ID="right_down")),
                              proj4string = CRS(proj4string(r.slope)))

  # Construct a polygon of the relevant slope
  # convert the lines into a polygon
  rsa.points.matrix <- matrix(ncol = 2, nrow = 0)
  for (j in (1:length(borderLines@lines))) {
    rsa.points.matrix <- rbind(rsa.points.matrix,
                               borderLines@lines[[j]]@Lines[[1]]@coords)
  }
  # convex hull, maybe not ideal. 'chull' returns the row numbers of the points
  # forming the hull
  rsa.hull <- rsa.points.matrix[chull(rsa.points.matrix),]
  # convert the hull to a polygon
  
  relevantSlopePolygon <- SpatialPolygons(list(Polygons(list(Polygon(rsa.hull)), ID=1)))
  
  # # select raster points in the polygon
  local.slope <- crop(r.slope, extent(gBuffer(relevantSlopePolygon, width = 50)))
  slope.points <- rasterToPoints(local.slope) # gives a matrix with x, y, and slope
  # convert slope points into SpatialPointsDataFrame
  slope.at.raster.points <- SpatialPointsDataFrame(coords = slope.points[,c('x', 'y')], 
                                                   data = data.frame(slope = slope.points[,'slope']))
  # determine the points inside the polygon of the relevant slope
  slope.points.in.rsa <- intersect(slope.at.raster.points, relevantSlopePolygon)
  max.slope.in.rsa <- max(slope.points.in.rsa@data$slope)
  max.slope.point.in.rsa <- slope.points.in.rsa[which.max(slope.points.in.rsa@data$slope),]
  mean.slope.in.rsa <- mean(slope.points.in.rsa@data$slope)
  # inverse distance weighted slope angle
  dist2track.point <- sqrt(rowSums((slope.points.in.rsa@coords - 
                         matrix(1, nrow = NROW(slope.points.in.rsa@coords)) %*% track.point)^2))
  dist2track.point[dist2track.point < 6] <- 6
  sum(slope.points.in.rsa@data$slope * (1/dist2track.point / sum(1/dist2track.point)))
  #extract(ras_temp.new, polygon, fun=mean, na.rm=TRUE, sp = T)
  
  return(slopeProps(borderLines = borderLines,
                    relevantSlopePolygon = relevantSlopePolygon,
                    localSlope = extract(r.slope, track.point),
                    localAspect = extract(r.aspect, track.point),
                    maxSlope = max.slope.point.in.rsa,
                    meanSlope = mean.slope.in.rsa))
}

