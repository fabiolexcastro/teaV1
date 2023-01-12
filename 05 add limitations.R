
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, rgeos, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
fles <- dir_ls('rf/output/run_1/results/raw', regexp = '.asc$') %>% as.character()
prob <- grep('Prob', fles, value = T) %>% terra::rast()
clst <- grep('Clust', fles, value = T) %>% terra::rast()

load('rData/run_1/clustereddata.rData')

occ
zone <- rbind(vect('shpf/base/kenya1.shp'), vect('shpf/base/uganda1.shp'))

# Extract the values ------------------------------------------------------
vles <- terra::extract(prob, occ[,c('x', 'y')]) %>% drop_na() %>% pull(2)
qntl <- quantile(vles, seq(0, 1, 0.01))
qntl <- as.data.frame(qntl)
qntl <- rownames_to_column(qntl)
thrs <- filter(qntl, rowname == '5%')
thrs <- pull(thrs, qntl)

save(thrs, file = 'rData/run_1/threshold_prob.rData')

# Limitations -------------------------------------------------------------
mtx_prb <- matrix(c(0, thrs, 0, thrs, 1, 2), ncol = 3, byrow = T)
mtx_cls <- matrix(c(0.5, 2 + 0.5, 0, 2 + 0.5, 5 + 0.5, 1), nrow = 2, byrow = T)

prob_rclf <- raster::reclassify(raster(prob), mtx_prb)
clst_rclf <- raster::reclassify(raster(clst), mtx_cls)

diff <- prob_rclf - clst_rclf
rslt <- raster(clst)
rslt[which(diff[] == -1)] <- 2 + 3 + 1
rslt[which(diff[] == 2)] <- 2 + 3 + 1

raster::writeRaster(x = rslt, filename = 'rf/output/run_1/results/process/rf_limitations.tif')


