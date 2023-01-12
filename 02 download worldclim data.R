
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, hrbrthemes, ggspatial, gtools, rgeos, geodata, stringr, rgbif, rnaturalearthdata, rnaturalearth, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
pnts <- read_csv('tble/points/malawi_rgbif.csv')
mlwi <- read_csv('tble/points/malawi.csv') %>% mutate(scientificName = 'Camelia sinensis')
pnts <- rbind(pnts, mlwi)

shpf <- vect('shpf/points1.shp')
isos <- c('UGA', 'KEN', 'MWI', 'IND')

indi <- gadm(country = 'IND', level = 0, path = 'tmpr')
keny <- gadm(country = 'KEN', level = 1, path = 'tmpr')

# Download ----------------------------------------------------------------
purrr::map(.x = 1:length(isos), .f = function(i){
  
  cat(isos[i], '\n')
  iso <- isos[i]
  tmx <- geodata::worldclim_country(country = iso, var = 'tmax', path = 'tmpr', version = '2.1')
  tmn <- geodata::worldclim_country(country = iso, var = 'tmin', path = 'tmpr', version = '2.1')
  ppt <- geodata::worldclim_country(country = iso, var = 'prec', path = 'tmpr', version = '2.1')
  bio <- geodata::worldclim_country(country = iso, var = 'bioc', path = 'tmpr', version = '2.1')
  
  dir <- 'tif/climate/wc/country'
  terra::writeRaster(x = tmx, filename = glue('{dir}/tmax_{iso}.tif'), overwrite = T)
  terra::writeRaster(x = tmn, filename = glue('{dir}/tmin_{iso}.tif'), overwrite = T)
  terra::writeRaster(x = ppt, filename = glue('{dir}/tmax_{iso}.tif'), overwrite = T)
  terra::writeRaster(x = bio, filename = glue('{dir}/bioc_{iso}.tif'), overwrite = T)
  cat('Done!\n')
  
})

# Shapefile to points -----------------------------------------------------
tbl1 <- shpf %>% crds() %>% as_tibble() %>% mutate(source = 'Google Earth')
tbl2 <- pnts %>% mutate(source = 'GBIF') %>% dplyr::select(x = Lon, y = Lat, source)
tble <- rbind(tbl1, tbl2)

# Extract the country ------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)
wrld <- vect(wrld)
cntr <- terra::extract(wrld, tble[,c('x', 'y')]) %>% dplyr::select(admin, adm0_a3) %>% as_tibble()
tble <- mutate(tble, country = cntr$admin, iso = cntr$adm0_a3)

mlwi <- geodata::gadm(country = 'MWI', level = 1, path = 'tmpr')
mlwi_pnts <- tble %>% filter(iso == 'MWI')

gpn <- ggplot() + 
  geom_sf(data = mlwi, fill = NA, col = 'grey30') + 
  geom_point(data = tble %>% filter(iso == 'MWI'), aes(x = x, y = y), col = 'darkred') + 
  coord_sf() + 
  theme_minimal() + 
  theme()

freq <- as.data.frame(table(tble$iso)) %>% arrange(desc(Freq)) %>% filter(Var1 %in% c('MWI', 'IND', 'UGA', 'KEN'))
tble <- tble %>% filter(iso %in% c('MWI', 'IND', 'UGA', 'KEN'))

zone <- wrld[wrld$adm0_a3 %in% c('MWI', 'IND', 'UGA', 'KEN'),]
plot(zone)
points(tble$x, tble$y, pch = 16, col = 'red')

# Filtering uganda, malawi, india -----------------------------------------
zone <- wrld[wrld$sov_a3 %in% c('MWI', 'KEN', 'UGA', 'IND')]
tble <- filter(tble, iso %in% c('MWI', 'KEN', 'UGA', 'IND'))
table(tble$iso)

# Join uganda and kenya ---------------------------------------------------
fles <- dir_ls('tif/climate/wc/country', regexp = '.tif') %>% grep(paste0(c('UGA', 'KEN'), collapse = '|'), ., value = T)
bioc <- grep('bioc', fles, value = T)
bioc <- as.character(bioc)
bioc <- map(bioc, rast)
bioc <- terra::sprc(bioc)
bioc <- terra::mosaic(bioc)
area <- zone[zone$sov_a3 %in% c('KEN', 'UGA'),]
bioc <- terra::crop(bioc, area) 
bioc <- terra::mask(bioc, area)
names(bioc) <- glue('bioc_{1:19}')

terra::writeRaster(x = bioc, filename = 'tif/climate/wc/ken_uga/bioc.tif', overwrite = TRUE)
write.csv(tble, 'tble/points/points_v1.csv', row.names = FALSE)

# Remove duplicated by cell  ----------------------------------------------
tble <- distinct(tble)

# Map by each country  ----------------------------------------------------
purrr::map(.x = 1:4, .f = function(i){
  
  cat('Start\n')
  tbl <- filter(tble, iso == isos[i])
  lim <- geodata::gadm(country = isos[i], level = 1, path = 'tmpr')
  lim <- st_as_sf(lim)
  
  gpn <- ggplot() + 
    geom_sf(data = lim, fill = NA, col = 'grey60', lwd = 0.3) + 
    geom_point(data = tbl, aes(x = x, y = y), col = 'darkred', size = 0.9) + 
    coord_sf() + 
    labs(x = 'Lon', y = 'Lat') +
    theme_minimal() + 
    theme()
  
  ggsave(plot = gpn, filename = glue('png/maps/location_tea_{isos[i]}.png'), units = 'in', width = 9, height = 7, dpi = 300)
  cat('Done!\n')
  
})

write.csv(tble, 'tble/points/points_v2.csv', row.names = FALSE)

# To extract the values ---------------------------------------------------
vles <- purrr::map_dfr(.x = 1:length(isos), .f = function(i){
  
  cat(isos[i], '\n')
  iso <- isos[i]
  fls <- dir_ls('tif/climate/wc/country', regexp = '.tif$')
  fls <- as.character(fls)
  bio <- grep('bioc', fls, value = T) 
  bio <- grep(iso, bio, value = T)
  bio <- terra::rast(bio)
  names(bio) <- glue('bioc_{1:19}')
  tbl <- filter(tble, iso == isos[i])
  vls <- terra::extract(bio, tbl[,c('x', 'y')])
  tbl <- cbind(tbl, vls)
  tbl <- as_tibble(tbl)
  cat('Done!\n')
  return(tbl)
  
})

write.csv(vles, 'tble/points/points_v2.csv', row.names = FALSE)


