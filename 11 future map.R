
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, rnaturalearth, rmapshaper, ggspatial, gtools, rgeos, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
rstr <- terra::rast('rf/output/run_1/results/process/ssp370/rf_modal_2021_2040.tif')
zone <- rbind(vect('shpf/base/kenya1.shp'), vect('shpf/base/uganda1.shp'))
tble <- terra::as.data.frame(rstr, xy = T)
tble <- as_tibble(tble)
colnames(tble)[3] <- 'layer'

# To make the map ---------------------------------------------------------
lbls <- tibble(layer = 1:7, class = c('Unsuitable', 'Unsuitable', 'Type 1', 'Type 2', 'Type 3', 'Limitations', 'Mixed'))
tble <- inner_join(tble, lbls, by = 'layer')
tble <- mutate(tble, class = factor(class, levels = unique(lbls$class)))

wrld <- ne_countries(returnclass = 'sf', scale = 50)
cntr <- ms_dissolve(input = st_as_sf(zone), field = 'COUNTRY')

lkes <- 'D:/DATA/world/ne_50m_lakes/ne_50m_lakes.shp' %>% st_read() 
library(colourpicker); colourWidget()

tble$class %>% unique()
gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = class)) + 
  scale_fill_manual(values = c('white', '#E6E6E6', '#F2EFC2')) +
  geom_sf(data = lkes, fill = '#4A869E', col = 'grey30') +
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey60', lwd = 0.4) +
  geom_sf(data = wrld, fill = NA, col = 'grey70') +
  geom_sf(data = cntr, fill = NA, col = 'grey30', lwd = 10.8) +
  labs(x = 'Lon', y = 'Lat', fill = 'AEZ') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) + 
  theme_bw() + 
  theme(legend.position = c(0.1, 0.13),
        legend.title = element_text(face = 'bold'),
        legend.background = element_rect(fill = 'white'),
        text = element_text(family = 'georg', color = 'grey50'),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  annotation_scale(location =  "br", width_hint = 0.5, text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = gmap, filename = glue('png/maps/run_1/suitability_future.png'), units = 'in', width = 9, height = 7, dpi = 300)


