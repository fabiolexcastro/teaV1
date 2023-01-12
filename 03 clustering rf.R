
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, hrbrthemes, ggspatial, gtools, rgeos, geodata, stringr, rgbif, rnaturalearthdata, rnaturalearth, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('FunctionsRF.R')

# Functions ---------------------------------------------------------------
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  
  datRF_presences <- occ[,3:ncol(occ)]
  print(nrow(datRF))
  
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  
  return(list(labelRF, clusterdata))
  
}


# Load data ---------------------------------------------------------------
bioc <- terra::rast('tif/climate/wc/ken_uga/country/bioc_all.tif')
tble <- read_csv('tble/points/points_v2.csv')
tble <- tble[,c(1, 2, 7:25)]
occ <- tble
names(bioc) <- glue('bioc_{1:19}')

# Clustering --------------------------------------------------------------
env_values <- as.matrix(occ[,3:ncol(occ)]); nrow(env_values)
datRF <- as.data.frame(occ[,3:ncol(occ)]); nrow(datRF) 
d <- dist(datRF, method = "euclidean")  
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 3)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,3:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 3

dir.create('rData/run_1', recursive = T)
run <- 'run_1'
save(datRF, file = paste0('rData/', run, '/datRF.rData'))
save(clusterdata, file = paste0('rData/', run, '/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('rData/', run, '/clustereddata.rData'))

# Make mosaic  ------------------------------------------------------------
fles <- dir_ls('tif/climate/wc/country') %>% grep('bioc', ., value = T) %>% as.character()
bioc <- map(fles, rast)
bioc <- sprc(bioc)
bioc <- mosaic(bioc)
wrld <- ne_countries(returnclass = 'sf', scale = 50)
wrld <- filter(wrld, sov_a3 %in% c('MWI', 'UGA', 'IND', 'KEN'))
bioc <- terra::crop(bioc, vect(wrld)) %>% terra::mask(., vect(wrld))
terra::writeRaster(x = bioc, filename = glue('tif/climate/wc/country/bioc_all.tif'), overwrite = T)
