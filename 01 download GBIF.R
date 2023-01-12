
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, hrbrthemes, ggspatial, gtools, rgeos, geodata, stringr, rgbif, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
spce <- 'Camellia sinensis'
tble <- occ_search(scientificName = spce, hasCoordinate = TRUE)
tble <- tble$data

wrld <- ne_countries(returnclass = 'sf', type = 'countries', scale = 50)

# A simple plot -----------------------------------------------------------
plot(st_geometry(wrld), border = 'grey60')
points(tble$decimalLongitude, tble$decimalLatitude, pch = 16, col = 'darkred')
write.csv(tble, 'tble/points/rgbif_raw.csv', row.names = FALSE)

# Malawi presences --------------------------------------------------------
pnts <- read_csv('tble/points/malawi.csv')
malw <- geodata::gadm(country = 'MWI', level = 0, path = 'tmpr')

# A simple plot
plot(malw, border = 'grey60')
points(pnts$Lon, pnts$Lat, pch = 16, col = 'brown')

tble <- tble %>% dplyr::select(2, decimalLongitude, decimalLatitude)
pnts <- mutate(pnts, scientificName = 'Camelia sinensis')
colnames(tble) <- c('scientificName', 'Lon', 'Lat')
tble <- relocate(tble, scientificName)
pnts <- relocate(pnts, scientificName)
pnts <- rbind(pnts, tble)

write.csv(pnts, 'tble/points/malawi_rgbif.csv', row.names = FALSE)

# To make the general map  ------------------------------------------------

gmain <- ggplot() + 
  geom_sf(data = wrld, fill = NA, col = 'grey60', lwd = 0.5) + 
  geom_point(data = pnts, aes(x = Lon, y = Lat), col = 'brown', size = 0.6) + 
  coord_sf() + 
  ggtitle(label = 'Tea location (GBIF)') +
  labs(x = '', y = '') +
  theme_ipsum_ps() + 
  theme(plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, color = 'grey50')) +
  annotation_scale(location =  "br", width_hint = 0.5, text_col = 'grey70', bar_cols= c('white', 'grey80')) +
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_col = 'grey70'))
gmain

ggsave(plot = gmain, filename = 'png/maps/location_tea_world.png', units = 'in', width = 8, height = 5.5, dpi = 300)

# Malawi points as kml file -----------------------------------------------
pnts <- st_as_sf(pnts, coords = c('Lon', 'Lat'), crs = st_crs(4326))
pnts
pnts
plot(st_geometry(pnts))
dir_create('kml')
st_write(pnts, 'kml/malawi_points.kml')

# Rest of points as kml  --------------------------------------------------
tble
shpf <- st_as_sf(tble, coords = c('Lon', 'Lat'), crs = st_crs(4326))
shpf
st_write(shpf, 'kml/world_points.kml')

st_write(shpf, 'shpf/world_points.shp')


ugan <- geodata::gadm(country = 'UGA', level = 1, path = 'tmpr')
dir_create('shpf/base')
writeVector(ugan, filename = 'shpf/base/uganda1.shp')
