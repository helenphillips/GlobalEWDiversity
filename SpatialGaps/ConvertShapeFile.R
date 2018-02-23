library(raster)
library(rgdal)

## Temperature
p <- shapefile('I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialGaps\\GlobalLayers\\2080-2099.temp.anom.a2\\pcmdi_long_anom.tas.sresa2.2080-2099.median.shp')
p
r <- raster(res(acc))
ext <- extent(p)
r <- raster(ext, res= 0.0083)  
projection(r) <- proj4string(p)
rr <- rasterize(p, field="annual", r, filename = 'I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialGaps\\GlobalLayers\\rasters2\\temp.tif')

## Precipitation

p <- shapefile('I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialGaps\\GlobalLayers\\2080-2099.precip.anom.a2\\futureB.ppt.totals.median.shp')
p
ext <- extent(p)
r <- raster(ext, res= 0.0083)  
projection(r) <- proj4string(p)
rr <- rasterize(p, field="ANNUAL", r, filename = 'I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialGaps\\GlobalLayers\\rasters2\\precip.tif')
