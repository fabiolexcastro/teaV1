
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, rgeos, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
fles <- dir_ls('rf/output/run_1/results/process/ssp370/2021_2040') %>% as.character()
mdls <- fles %>% str_split(pattern = '_') %>% map(5) %>% unlist() %>% unique() %>% gsub('.tif$', '', .)
zone <- rbind(vect('shpf/base/kenya1.shp'), vect('shpf/base/uganda1.shp'))

# Function mixed ----------------------------------------------------------
load('rData/run_1/threshold_prob.rData')
thr_prb <- thrs
load('rData/run_1/threshold_unc.rData')
thr_unc <- thrs

# Function ----------------------------------------------------------------
calcMixed <- function(mdl){
  
  cat(mdl, '\n')
  lim <- grep(paste0(mdl, '.tif'), fles, value = T) %>% grep('lim', ., value = T) %>% raster()
  unc <- dir_ls('rf/output/run_1/results/raw/ssp370/2021-2040') %>% grep(paste0(mdl, '.tif'), ., value = T) %>% as.character() %>% grep('uncr', ., value = T) %>% raster()
  prb <- dir_ls('rf/output/run_1/results/raw/ssp370/2021-2040') %>% grep(paste0(mdl, '.tif'), ., value = T) %>% as.character() %>% grep('prob', ., value = T) %>% raster()
  rslt <- lim
  rslt[which(unc[] < thr_unc & prb[] > thr_prb)] <- 7
  raster::writeRaster(x = rslt, filename = glue('rf/output/run_1/results/process/ssp370/2021_2040/rf_unc_{mdl}.tif'), overwrite = T)
  
}

map(mdls[1:length(mdls)], calcMixed)

# Calc modal --------------------------------------------------------------
fles <- dir_ls('rf/output/run_1/results/process/ssp370/2021_2040') %>% grep('unc', ., value = T) %>% as.character()
rstr <- terra::rast(fles)
plot(rstr)

mdal <- terra::modal(rstr)
terra::writeRaster(x = mdal, filename = 'rf/output/run_1/results/process/ssp370/rf_modal_2021_2040.tif', overwrite = TRUE)

# Add the NA --------------------------------------------------------------
which.lyr(is.na(mdal))
mdal <- raster(mdal)
mdal[is.na(mdal)] <- 1
mdal <- terra::rast(mdal) %>% terra::crop(., zone) %>% terra::mask(., zone)

terra::writeRaster(x = mdal, filename = 'rf/output/run_1/results/process/ssp370/rf_modal_2021_2040.tif', overwrite = TRUE)


