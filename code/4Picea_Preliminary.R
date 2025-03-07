dev.off()
rm(list=ls())


#install.packages("elevatr")
library(terra)
library(raster)
library(elevatr)
#library(rgdal)
library(sf)

#shape <- read_sf(dsn = "~/Google Drive/Shared Drives/4Picea/4Picea/gsp/", layer = "Spruce_plantation")
shape <- read_sf(dsn = "G:/Shared Drives/4Picea/4Picea/gsp/", layer = "Spruce_plantation")

plot(shape)

pie.ft <- get_elev_raster(shape,z=14,override_size_check = TRUE)
plot(pie.ft)
plot(shape,add=TRUE)

#pie.ft
#f <- drawExtent()

picea <- crop(pie.ft,shape)
plot(picea)
plot(shape,add=TRUE)

getwd()

#writeRaster(picea,"~/Google Drive/Shared Drives/4Picea/4Picea/gsp/",overwrite=TRUE)
writeRaster(picea,"C:/Users/ashley.lynn.carter/Documents/4PICEA_DEM.tif",overwrite=TRUE)

#test <- raster("~/Google Drive/Shared Drives/4Picea/4Picea/gsp/")
test <- raster("C:/Users/ashley.lynn.carter/Documents/4PICEA_DEM.tif")
plot(test)

demSt <- terrain(test,opt=c("slope","aspect","TPI","TRI","roughness","flowdir"),unit='degrees',neighbors=8)
plot(demSt)

demSt
plot(demSt$slope)
plot(shape,add=TRUE)

demSt
plot(demSt$aspect)
plot(shape,add=TRUE)

library(topmodel)
library(MEForLab)
bob <- create_layers(test)
plot(bob$twi)
plot(shape,add=TRUE)

picea.stack <- raster::stack(picea,demSt,bob)
plot(picea.stack)


#writeRaster(new.cov,"C:/Users/premerm/Desktop/batch/PNW_Covariates_v4.tif",options="INTERLEAVE=BAND",overwrite=TRUE)

#terra::writeRaster(picea.stack,"~/Desktop/picea_stack.tif",options="INTERLEAVE=BAND",overwrite=TRUE)
terra::writeRaster(picea.stack,"C:/Users/ashley.lynn.carter/Documents/4PICEA_DEM.tif",options="INTERLEAVE=BAND",overwrite=TRUE)

#pull from google drive

ft <-raster::stack("G:/Shared drives/4Picea/4Picea/gsp/picea_stack.tif")
ft
plot(ft)     
names(ft) <- c("elevation","tri","tpi","roughness","slope","aspect",
               "flowdir","filled.elevations","upslope.area","twi")
plot(ft)
plot(ft$aspect)
plot(shape,add=TRUE,col="darkred",lwd=5)


