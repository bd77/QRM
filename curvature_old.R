# curvature

setwd("C:/Documenten/QRM/Rscripts")
library(raster)
source("DEMderiv.R")
library(rasterVis)
library(RColorBrewer)

dem.alps <- raster("C:/Documenten/slopeangle/EUD_CP-DEMS_4500025000-AA.tif") 

simplon.extent <- extent(c(4164000, 4170000, 2571000, 2575000))
dem.simplon <- crop(dem.alps, simplon.extent)
plot(dem.simplon)

curvature.simplon <- curvature(dem.simplon, s=5, type="bolstad")
curvature.simplon <- DEMderiv(dem.simplon, attr = "plan.curvature", method = "evans")

quantile(curvature.simplon)

levelplot(curvature.simplon, breaks =)
breaks <- as.vector(quantile(curvature.simplon))

levelplot(curvature.simplon, at = breaks, 
          margin=FALSE)
