

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, outliers, Hmisc, cclust, sf, randomForest, multcomp, dismo, magrittr, ggpubr, corrplot)

rm(list = ls())
source('FunctionsRF.R')
run <- 'run_1'
myproj <- CRS('+proj=longlat +datum=WGS84')
options(stringsAsFactors = FALSE)

# Functions ---------------------------------------------------------------
rf.clust <-function(occ, nforest, ntrees, nVars, nclasses){
  # occ = back_swd; nforest = 50; ntrees = 500; nVars = 8; nclasses = 2
  datRF_presences <- occ[,3:ncol(occ)] %>% as.data.frame()
  print(nrow(datRF_presences))
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
bio <- raster::stack('tif/climate/wc/country/bioc_all.tif')
msk <- raster('tif/climate/wc/country/bioc_all.tif') * 0
load(file = 'rData/run_1/clustereddata.rData')
occ
wrl <- ne_countries(returnclass = 'sf', scale = 50) %>% vect()
sub <- mutate(occ, country = pull(terra::extract(wrl, occ[,c('x', 'y')]), 'sov_a3')) %>% filter(country %in% c('KEN', 'UGA') )
table(occ$country)

# Bias raster - process ---------------------------------------------------
SPspecies <- SpatialPoints(occ[,1:2]) 
crs(SPspecies) <- myproj
back_raster <- msk
speciescell <- raster::extract(msk, SPspecies, cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA #remove the cell with presences
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- max(clusteredpresdata$cluster) 
ratio <- NumberOfClusters/1
numberofpresences <- nrow(clusteredpresdata) 
crs(back_raster) <- myproj
back <- randomPoints(back_raster, 1*numberofpresences) %>% # Cambiar al 2
  as_data_frame()
coordinates(back) <- ~ x + y
lyr <- bio
names(lyr) <- glue('bioc_{1:19}')
back_swd <- raster::extract(lyr, back) %>% 
  cbind(coordinates(back), .)
nrow(back_swd) == nrow(back_swd[complete.cases(back_swd),])
write.csv(back_swd, 'tble/points/back_swd.csv', row.names = FALSE)
write.csv(occ, 'tble/points/occ_swd.csv', row.names = FALSE)

back_swd <- as.data.frame(back_swd)

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500 # Cambari el 8 al default
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
detach(datRF)

# classdata <- cbind(pb = as.factor(c(rep(1, 2000), rep(2, 1148))), back_swd[,3:ncol(back_swd)])
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, bioc_1:bioc_19)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background

dim(classdata_2); dim(presvalue_swd)
presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
write.csv(allclasses_swd, '../tbl/points/run5/all_classes_swd_v3.csv', row.names = FALSE)
write.csv(allclasses_swd, 'tble/points/all_classes_swd_v1.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
vrs <- gsub('.tif', '', vrs) 
vrs <- gsub('\\$', '', vrs)
vrs <- paste0('bioc_', 1:19)
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)
library(pROC)

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) # nodesize = 1
  
  save(rfmodel, file = paste('rf/output/run5/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}


print('Load libraries')
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, outliers, Hmisc, cclust, sf, randomForest, multcomp, dismo, magrittr, ggpubr, corrplot)


rm(list = ls())
source('FunctionsRFclustering.R')
run <- 'run5'
myproj <- CRS('+proj=longlat +datum=WGS84')
options(stringsAsFactors = FALSE)

print('Functions to use')
rf.clust <-function(occ, nforest, ntrees, nVars, nclasses){
  # occ = back_swd; nforest = 50; ntrees = 500; nVars = 8; nclasses = 2
  datRF_presences <- occ[,3:ncol(occ)] %>% as.data.frame()
  print(nrow(datRF_presences))
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
shp <- shapefile('../shp/base/dptos.shp')
# shp <- shp[shp@data$NOMBRE_DPT %in% 'RISARALDA',]
msk <- raster('../climate/current/wc_v20_v2/bio_1.tif') * 0
vrs <- paste0('bio_', 1:19, '.tif')

load(file = '../rData/run5/clustereddata.rData')

vrs <- colnames(occ[grep('bio', colnames(occ), value = F)]) %>% paste0(., '.asc$')
lyr <- list.files('../climate/current/wc_v20_v2', full.names = TRUE, pattern = '.tif$') %>%
  mixedsort() %>% 
  # grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>%
  stack()

msk <- lyr[[1]] * 0

# Bias raster - process ---------------------------------------------------
SPspecies <- SpatialPoints(occ[,1:2]) 
crs(SPspecies) <- myproj
back_raster <- msk
speciescell <- raster::extract(msk, SPspecies, cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA #remove the cell with presences
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- max(clusteredpresdata$cluster) 
ratio <- NumberOfClusters/1
numberofpresences <- nrow(clusteredpresdata) 
crs(back_raster) <- myproj
back <- randomPoints(back_raster, 2*numberofpresences) %>% # Cambiar al 2
  as_data_frame()
coordinates(back) <- ~ x + y
back_swd <- raster::extract(lyr, back) %>% 
  cbind(coordinates(back), .)
nrow(back_swd) == nrow(back_swd[complete.cases(back_swd),])
write.csv(back_swd, '../tbl/points/run5/back_swd.csv', row.names = FALSE)
write.csv(occ, '../tbl/points/run5/occ_swd.csv', row.names = FALSE)

plot(shp, border = 'grey')
points(back_swd[,1], back_swd[,2], pch = 16, col = 'red')
back_swd <- back_swd[complete.cases(back_swd),]
back_swd <- as.data.frame(back_swd)

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500 # Cambari el 8 al default
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, bio_1:bio_19)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background

dim(classdata_2); dim(presvalue_swd)
presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
write.csv(allclasses_swd, '../tbl/points/run4/all_classes_swd.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
vrs <- gsub('.tif', '', vrs) 
vrs <- gsub('\\$', '', vrs)
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)
library(pROC)

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) # nodesize = 1
  
  save(rfmodel, file = paste('rf/output/run_1/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

dir.create('rData/run_1')

save(rflist, file = paste('rData/run_1/', '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('rData/run_1', '/importanceRF.rData'))
save(auc, file = paste0('rData/run_1', '/aucRF_dist.rData'))
save(rff, file = paste0('rData/run_1', '/rff_dist.rData'))

# Predict modell
zone <- rbind(terra::vect('shpf/base/kenya1.shp'), terra::vect('shpf/base/uganda1.shp'))
lyr <- raster::crop(lyr, as(zone, 'Spatial')) %>% raster::mask(., as(zone, 'Spatial'))
climatevalues  <- data.frame(getValues(lyr))
NumberOfClusters <- 3

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- lyr[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- lyr[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- lyr[[1]]
values(rasterRFclass) <- rasterRF

plot(rasterRF)

dir.create('rf/output/run5/results/raw', recursive = TRUE)
writeRaster(rasterRFclass, paste0('../rf/output/run5/results/raw/RF_5Clust_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0('../rf/output/run5/results/raw/RF_5Prob_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste0('../rf/output/run5/results/raw/RF_5Unc_current.asc'), format = 'ascii', overwrite = T)

print('To make the map')
crn <- rasterToPoints(rasterRFclass, xy = TRUE) %>% 
  as_tibble() %>% 
  gather(var, value, -x, -y) %>% 
  inner_join(., data.frame(var = 1:5, type = c('Unsuitable', 'Unsuitable', 'Type 1', 'Type 2', 'Type 3')), by = c('value' = 'var'))

ext <- extent(zmb[zmb@data$NAME_1 %in% 'Northern',])
nth <- zmb[zmb@data$NAME_1 %in% 'Northern',]

g1 <- ggplot(data = crn) +
  geom_tile(aes(x = x, y = y, fill = factor(value))) + 
  geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
  coord_equal() +
  scale_fill_manual(values = c('grey', 'grey', 'brown', '#868A08', '#0B3B17'),
                    labels = c('Unsuitable', 'Unsuitable', paste0('Cluster ', 1:3)),
                    name = 'AEZ') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom') +
  labs(x = 'Longitude',
       y = 'Latitude', 
       caption = 'Baseline')
ggsave(plot = g1, filename = '../png/maps/run4/crn_v1.png', units = 'in', width = 9, height = 6, dpi = 300)


auc <- unlist(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

save(rflist, file = paste('../rData/run4/', '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('../rData/run4/', '/importanceRF.rData'))
save(auc, file = paste0('../rData/run4/', '/aucRF_dist.rData'))
save(rff, file = paste0('../rData/run4/', '/rff_dist.rData'))

# Predict modell
climatevalues  <- data.frame(getValues(lyr))
NumberOfClusters <- 3

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- lyr[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- lyr[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- lyr[[1]]
values(rasterRFclass) <- rasterRF

dir_create('rf/output/run_1/results/raw')
writeRaster(rasterRFclass, paste0('rf/output/run_1/results/raw/RF_3Clust_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0(' rf/output/run_1/results/raw/RF_3Prob_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste0('rf/output/run_1/results/raw/RF_3Unc_current.asc'), format = 'ascii', overwrite = T)

print('To make the map')
crn <- rasterToPoints(rasterRFclass, xy = TRUE) %>% 
  as_tibble() %>% 
  gather(var, value, -x, -y) %>% 
  inner_join(., data.frame(var = 1:5, type = c('Unsuitable', 'Unsuitable', 'Type 1', 'Type 2', 'Type 3')), by = c('value' = 'var'))

table(crn$type)

g1 <- ggplot() +
  geom_tile(data = crn, aes(x = x, y = y, fill = factor(value))) + 
  geom_sf(data = st_as_sf(zone), color = 'grey', fill = 'NA') +
  coord_sf() +
  scale_fill_manual(values = c('grey', 'brown', '#868A08', '#0B3B17'),
                    labels = c('Unsuitable', paste0('Cluster ', 1:3)),
                    name = 'AEZ') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom') +
  labs(x = 'Longitude',
       y = 'Latitude', 
       caption = 'Baseline')
ggsave(plot = g1, filename = 'png/maps/run_1/crn_v1.png', units = 'in', width = 9, height = 6, dpi = 300)


