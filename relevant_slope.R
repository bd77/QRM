# ---------------------------------------------------
# function to determine the RSA (relevant slope area)
# ---------------------------------------------------


# auxiliary functions
run.up <- function(start.point, r.aspect, step, max.dist) {
  # determine the point step meters uphill from the input point
  dist2start <- 0
  uphill.line <- start.point
  current.point <- start.point
  n.iter <- 0
  max.iter <- 2*max.dist/step
  while (dist2start <= max.dist & n.iter < max.iter) {
    aspect.cp <- extract(r.aspect, current.point)
    next.point <- current.point + step * c(-sin(aspect.cp), -cos(aspect.cp))
    dist2start <- sqrt(sum((start.point - current.point)**2))
    uphill.line <- rbind(uphill.line, next.point)
    current.point <- next.point
    n.iter <- n.iter+1
  }
  return(uphill.line)
}

rsa <- function(r.elev, r.slope, r.aspect, r.pcurv, track.point) {
  # r.elev: elevation raster
  # r.slope: raster with the slope angle in degrees
  # r.aspect: raster with the aspect in radians, N=0, E=pi/2, S=pi, W=3*pi/2
  # r.pcurv: plan curvature, <0 = concave (valley), >0 = convex (ridge)
  # track.point: 1x2 matrix with coordinates
  # the coordinate system should be orthogonal
  
  # use convex hull function "chull" to make an area out of the points
  
  max.dist.up <- 500
  max.dist.down <- 50
  max.dist.side <- 50
  step <- 25
  
  # steepest way up from the track point 
  uphill.line.from.track.point <- run.up(track.point, r.aspect, step, max.dist.up)
  
  # make a step left and right
  aspect.tp <- extract(r.aspect, track.point)
  point.left <- track.point + max.dist.side * c(-cos(aspect.tp), sin(aspect.tp))
  point.right <- track.point + max.dist.side * c(cos(aspect.tp), -sin(aspect.tp))
  
  # steepest way up from the right and left point 
  uphill.line.from.left <- run.up(point.left, r.aspect, step, max.dist.up)
  uphill.line.from.right <- run.up(point.right, r.aspect, step, max.dist.up)
  
  all.lines <- SpatialLines(list(Lines(Line(uphill.line.from.track.point), ID="centre"),
                                 Lines(Line(uphill.line.from.left), ID="left"),
                                 Lines(Line(uphill.line.from.right), ID="right")))
  
  return(all.lines)
}

