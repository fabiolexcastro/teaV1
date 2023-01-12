
# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, magrittr)

# Initial Setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run_1'

# Load Data
crn <- raster('rf/output/run_1/results/process/rf_unc.tif')
f50 <- raster('rf/output/run_1/results/process/ssp370/rf_modal_2021_2040.tif')

all_options <- read_csv('tble/classesImpGraLimMix.csv')
unique(all_options$category) 
labelss <- data.frame(value = c(0, 1, 2, 3, 4, 5), category = c('Unsuit', 'cope', 'adjust', 'transform', 'opportunity', 'resilience'))

# Function to use
impGra <- function(crn, ftr){
  
  crn <- crn
  ftr <- f50
  
  msk <- crn * 0
  crd_df <- coordinates(crn)
  
  x <- raster::extract(crn, crd_df, cellnumbers = TRUE) %>% as_data_frame()
  ncell <- dplyr::select(x, cells)
  x <- select_(x, names(crn))
  colnames(x) <- 'current'
  
  y <- raster::extract(ftr, crd_df[,c('x', 'y')], cellnumbers = TRUE) %>% as_data_frame()
  y <- select_(y, names(ftr))
  colnames(y) <- 'future'
  
  z <- data.frame(x, y, ncell) %>% as_tibble()
  
  print('To Results')
  rslts <- left_join(z, all_options, by = c('current', 'future'))
  labls <- as_tibble(labelss) %>% mutate(category = as.character(category))
  
  final <- full_join(rslts, labls, by = 'category') %>%
    dplyr::select(value) %>%
    pull(1)
  
  final <- left_join(rslts, labls, by = 'category') %>%
    dplyr::select(value) %>%
    pull(1)
  
  length(final)
  length(msk)
  hist(final)
  
  rst <- raster::setValues(msk, final)
  return(rst)
  
}

hist(rst[])

writeRaster(rst, 'rf/output/run_1/results/process/rf_impact_30.tif', overwrite = TRUE) 

# Impact gradient map -----------------------------------------------------
tbl <- rasterToPoints(rst, spatial = FALSE) %>% as_tibble()
colnames(tbl)[3] <- 'value'

unique(tbl$value)

lbl <- tibble(value = c(0, 3, 5), class = c('Unsuitable', 'Transformation', 'Systemic resilience'))
tbl <- inner_join(tbl, lbl, by = 'value')
tbl <- mutate(tbl, class = factor(class, levels = c('Unsuitable', 'Transformation', 'Systemic resilience')))

wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- rbind(vect('shpf/base/kenya1.shp'), vect('shpf/base/uganda1.shp'))
cntr <- ms_dissolve(input = st_as_sf(zone), field = 'COUNTRY')
lkes <- 'D:/DATA/world/ne_50m_lakes/ne_50m_lakes.shp' %>% st_read() 

gmap <- ggplot() + 
  geom_tile(data = tbl, aes(x = x, y = y, fill = class)) + 
  scale_fill_manual(values = c('#FFFFFF', '#9E4C4C', '#E0CD3D')) +
  geom_sf(data = lkes, fill = '#4A869E', col = 'grey30') +
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey60', lwd = 0.4) +
  geom_sf(data = wrld, fill = NA, col = 'grey70') +
  geom_sf(data = cntr, fill = NA, col = 'grey30', lwd = 10.8) +
  labs(x = 'Lon', y = 'Lat', fill = 'Gradient') +
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

ggsave(plot = gmap, filename = glue('png/maps/run_1/impact_gradient_2021_2040.png'), units = 'in', width = 9, height = 7, dpi = 300)


# Models ------------------------------------------------------------------

mdls <- dir_ls('rf/output/run_1/results/process/ssp370/2021_2040') %>% as.character %>% grep('lim', ., value = T) %>% basename() %>% str_split(pattern = '_') %>% map(3) %>% gsub('.tif', '', .) 
write.csv(mdls, 'test.csv')

