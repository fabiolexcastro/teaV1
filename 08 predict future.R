
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, glue, rgeos, gtools, stringr, fs, outliers, Hmisc, cclust, sf, randomForest, multcomp, dismo, magrittr, ggpubr, corrplot)

rm(list = ls())
source('FunctionsRF.R')
run <- 'run_1'
myproj <- CRS('+proj=longlat +datum=WGS84')
options(stringsAsFactors = FALSE)

# Load data ---------------------------------------------------------------
load('rData/run_1/rflist_3.rdata')
rff <- do.call(randomForest::combine, rflist)
mdls <- c( "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")
fles <- dir_ls('tif/climate/c6/370/2021-2040', regexp = '.tif$') %>% as.character()

# Function ----------------------------------------------------------------
myPredict <- function(mdl){
  
  mdl <- mdls[18]
  mdl <- 'EC-Earth3-Veg-LR'
  fls <- grep(mdl, fles, value = T)
  ppt <- grep('prec', fls, value = T) %>% stack()
  tmx <- grep('tmax', fls, value = T) %>% stack()
  tmn <- grep('tmin', fls, value = T) %>% stack()
  
  cat('To calculate the bioclimatic variables\n')
  bio <- dismo::biovars(prec = ppt, tmin = tmn, tmax = tmx)
  writeRaster(x = bio, filename = glue('tif/climate/c6/370/2021-2040/{mdl}_bioc.tif'), overwrite = T)
  
  # bio <- stack(grep('bioc', fls, value = T))
  names(bio) <- glue('bioc_{1:19}')
  
  climatevalues <- data.frame(getValues(bio))
  
  cat('Climate values\n')
  rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
  NumberOfClusters <- 3
  rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
  uncertainty <- apply(rasterProbs, 1, max)  
  
  rasterRFprob <- bio[[1]]
  values(rasterRFprob) <- rasterRF 
  
  rasterRFuncertainty <- bio[[1]]
  values(rasterRFuncertainty) <- uncertainty 
  
  rasterRF <- max.col(rasterProbs, 'first')
  rasterRFclass <- bio[[1]]
  values(rasterRFclass) <- rasterRF
  
  plot(rasterRFclass)
  plot(rasterRFprob)
  plot(rasterRFuncertainty)
  
  dir <- 'rf/output/run_1/results/raw/ssp370/2021-2040/'; dir_create(dir)
  terra::writeRaster(x = rasterRFclass, filename = glue('{dir}/rf_class_{mdl}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFprob, filename = glue('{dir}/rf_prob_{mdl}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFuncertainty, filename = glue('{dir}/rf_uncr_{mdl}.tif'), overwrite = TRUE)
  
  
}


