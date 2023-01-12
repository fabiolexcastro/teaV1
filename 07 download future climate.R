
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(rmapshaper, rnaturalearthdata, rnaturalearth, dismo, ggspatial, terra, colourpicker, fs, sf, tidyverse, gtools, rgeos, raster, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
zone <- st_read('shpf/base/uga_ken.shp')
zone <- mutate(zone, gid = 1)
limt <- ms_dissolve(input = zone, field = 'gid')
mdls <- c( "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")
sspe <- '370'
time <- c('2021-2040', '2041-2060')

# To download -------------------------------------------------------------
download <- function(mdl, ssp, tme){
   
  # mdl <- mdls[1]
  # ssp <- sspe
  # tme <- time[1]
   
  cat(mdl, '\n')
  
  # Download
  prec_1 <- geodata::cmip6_tile(lon = 35, lat = 0, model = mdl, ssp = '370', time = tme, var = 'prec', path = 'tmpr')
  tmax_1 <- geodata::cmip6_tile(lon = 35, lat = 0, model = mdl, ssp = '370', time = tme, var = 'tmax', path = 'tmpr')
  tmin_1 <- geodata::cmip6_tile(lon = 35, lat = 0, model = mdl, ssp = '370', time = tme, var = 'tmin', path = 'tmpr')
  
  prec_2 <- geodata::cmip6_tile(lon = 35, lat = 2, model = mdl, ssp = '370', time = tme, var = 'prec', path = 'tmpr')
  tmax_2 <- geodata::cmip6_tile(lon = 35, lat = 2, model = mdl, ssp = '370', time = tme, var = 'tmax', path = 'tmpr')
  tmin_2 <- geodata::cmip6_tile(lon = 35, lat = 2, model = mdl, ssp = '370', time = tme, var = 'tmin', path = 'tmpr')
  
  prec_3 <- geodata::cmip6_tile(lon = 25, lat = 2, model = mdl, ssp = '370', time = tme, var = 'prec', path = 'tmpr')
  tmax_3 <- geodata::cmip6_tile(lon = 25, lat = 2, model = mdl, ssp = '370', time = tme, var = 'tmax', path = 'tmpr')
  tmin_3 <- geodata::cmip6_tile(lon = 25, lat = 2, model = mdl, ssp = '370', time = tme, var = 'tmin', path = 'tmpr')
  
  # Extract by mask
  prec_1 <- terra::crop(prec_1, vect(limt)) %>% terra::mask(., vect(limt))
  prec_2 <- terra::crop(prec_2, vect(limt)) %>% terra::mask(., vect(limt))
  prec_3 <- terra::crop(prec_3, vect(limt)) %>% terra::mask(., vect(limt))
  
  tmax_1 <- terra::crop(tmax_1, vect(limt)) %>% terra::mask(., vect(limt))
  tmax_2 <- terra::crop(tmax_2, vect(limt)) %>% terra::mask(., vect(limt))
  tmax_3 <- terra::crop(tmax_3, vect(limt)) %>% terra::mask(., vect(limt))
  
  tmin_1 <- terra::crop(tmin_1, vect(limt)) %>% terra::mask(., vect(limt))
  tmin_2 <- terra::crop(tmin_2, vect(limt)) %>% terra::mask(., vect(limt))
  tmin_3 <- terra::crop(tmin_3, vect(limt)) %>% terra::mask(., vect(limt))
  
  prec <- list(prec_1, prec_2, prec_3) %>% sprc() %>% mosaic()
  tmax <- list(tmax_1, tmax_2, tmax_3) %>% sprc() %>% mosaic()
  tmin <- list(tmin_1, tmin_2, tmin_3) %>% sprc() %>% mosaic()
  
  rm(prec_1, prec_2, prec_3, tmax_1, tmax_2, tmax_3, tmin_1, tmin_2, tmin_3)
  gc()
  
  # To create biovars
  prec <- stack(prec); tmax <- stack(tmax); tmin <- stack(tmin)
  # bioc <- dismo::biovars(prec = prec, tmin = tmin, tmax = tmax)
  
  # To write 
  dout <- glue('tif/climate/c6/{ssp}/{tme}')
  # raster::writeRaster(x = bioc, filename = glue('{dout}/{mdl}_bioc.tif'))
  raster::writeRaster(x = prec, filename = glue('{dout}/{mdl}_prec.tif'), overwrite = T)
  raster::writeRaster(x = tmax, filename = glue('{dout}/{mdl}_tmax.tif'), overwrite = T)
  raster::writeRaster(x = tmin, filename = glue('{dout}/{mdl}_tmin.tif'), overwrite = T)
  cat('Done!\n')
  
}

purrr::map(.x = 6:10, .f = function(i){download(mdl = mdls[i], ssp = '370', tme = '2021-2040')})
purrr::map(.x = 1:5, .f = function(i){download(mdl = mdls[i], ssp = '370', tme = '2041-2060')})

dir_ls('tif/climate/c6/370/2021-2040') %>% as.character %>% grep('prec', ., value = T)
