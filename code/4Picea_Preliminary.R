dev.off()
rm(list=ls())


#install.packages("elevatr")
library(elevatr)
library(rgdal)
library(sf)

shape <- read_sf(dsn = "~/Google Drive/Shared Drives/4Picea/4Picea/gsp/", layer = "Spruce_plantation")
plot(shape)

pie.ft <- get_elev_raster(shape,z=14,override_size_check = TRUE)
plot(pie.ft)
plot(shape,add=TRUE)

pie.ft
f <- drawExtent()

picea <- crop(pie.ft,f)
plot(picea)
plot(shape,add=TRUE)

getwd()

writeRaster(picea,"~/Desktop/4PICEA_DEM.tif",overwrite=TRUE)

test <- raster("~/Desktop/4PICEA_DEM.tif")
plot(test)

demSt <- terrain(test,opt=c("slope","aspect","TPI","TRI","roughness","flowdir"),unit='degrees',neighbors=8)
plot(demSt)

demSt
plot(demSt$slope)
plot(shape,add=TRUE)

library(topmodel)
bob <- create_layers(test)
plot(bob$twi)
plot(shape,add=TRUE)

picea.stack <- raster::stack(picea,demSt,bob)
plot(picea.stack)

#writeRaster(new.cov,"C:/Users/premerm/Desktop/batch/PNW_Covariates_v4.tif",options="INTERLEAVE=BAND",overwrite=TRUE)

terra::writeRaster(picea.stack,"~/Desktop/picea_stack.tif",options="INTERLEAVE=BAND",overwrite=TRUE)



ft <- raster::stack("~/Desktop/picea_stack.tif")


plot(ft)     
names(ft) <- c("elevation","tri","tpi","roughness","slope","aspect",
               "flowdir","filled.elevations","upslope.area","twi")
plot(ft)
plot(ft$aspect)
plot(shape,add=TRUE,col="darkred",lwd=5)
