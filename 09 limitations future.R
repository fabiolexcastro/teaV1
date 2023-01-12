
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, rgeos, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
fles <- dir_ls('rf/output/run_1/results/raw/ssp370/2021-2040') %>% as.character()
mdls <- fles %>% str_split(pattern = '_') %>% map(4) %>% unlist() %>% unique() %>% gsub('.tif$', '', .)
zone <- rbind(vect('shpf/base/kenya1.shp'), vect('shpf/base/uganda1.shp'))

# Function limitations ----------------------------------------------------
load('rData/run_1/threshold_prob.rData')

mtx_prb <- matrix(c(0, thrs, 0, thrs, 1, 2), ncol = 3, byrow = T)
mtx_cls <- matrix(c(0.5, 2 + 0.5, 0, 2 + 0.5, 5 + 0.5, 1), nrow = 2, byrow = T)

calcLimitations <- function(mdl){
  
  cat(mdl, '\n')
  fls <- grep(paste0(mdl, '.tif'), fles, value = T)
  cls <- grep('class', fls, value = T) %>% raster()
  prb <- grep('prob', fls, value = T) %>% raster()
  
  cls.rcl <- raster::reclassify(x = cls, rcl = mtx_cls)
  prb.rcl <- raster::reclassify(x = prb, rcl = mtx_prb)
  
  diff <- prb.rcl - cls.rcl
  rslt <- raster(cls)
  rslt[which(diff[] == -1)] <- 2 + 3 + 1
  rslt[which(diff[] == 2)] <- 2 + 3 + 1
  
  writeRaster(x = rslt, filename = glue('rf/output/run_1/results/process/ssp370/2021_2040/rf_lim_{mdl}.tif'), overwrite = T)
  
}

map(mdls[5:length(mdls)], calcLimitations)


