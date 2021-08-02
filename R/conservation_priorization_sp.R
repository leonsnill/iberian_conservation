# ====================================================================================================
#
# IBERIAN CONSERVATION
# Prioritization
#
# ====================================================================================================


# ----------------------------------------------------------------------------------------------------
# IMPORT
# ----------------------------------------------------------------------------------------------------
library(sf) # Simple Features for R
library(rnaturalearth) # World Map Data from Natural Earth
library(here) # A Simpler Way to Find Your Files
library(stars) # Spatiotemporal Arrays, Raster and Vector Data Cubes
library(dplyr) # A Grammar of Data Manipulation
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(ggnewscale) # Multiple Fill and Color Scales in 'ggplot2'
library(scico) # Colour Palettes Based on the Scientific Colour-Maps
library(geobgu) # install from GitHub ("michaeldorman/geobgu")
library(ggrepel) # Automatically Position Non-Overlapping Text Labels with 'ggplot2'
library(raster)
library(tidyr)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(rgeos) 
library(dismo) 
library(rnaturalearthdata)

#species <- "ursusarctos"
species <- "lynxpardinus"

# CRS of the project
crs_project <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"



# We are going to follow four different strategies to spatial prioritirize the conservation of both species
#1.	Identify top-ranking areas with highest suitability per species.
#Assess how much of this is already in protected areas Natura 2000
#2. Select those areas that will likely be colonized under the low dispersal scenarios (dispersal without much human pressure)
# Assess how much of this is adjacent or closeby protected areas 
# Identify areas that are potentially large enough for hosting populations if colonized (candidate sites for new protected areas?)
#3. Select areas that are colonized under the high human pressure and identify 
#Which of these cells fall into areas with high human density. These are areas where conflicts will become likely in the future (or already are)
#4. Select areas that are colonized under low dispersal limitations but not under high

# ----------------------------------------------------------------------------------------------------
# 1. PRIORITISE BASED ON SDM ####
# ----------------------------------------------------------------------------------------------------
r_ens <- raster(paste0("Data/SDMs/", species, "_SDM_ensemble_mn_md_wmn_cav_sdp.tif"))
r_ens_bin <- raster(paste0("Data/SDMs/", species, "_SDM_ensembles_binary_mn_md_wmn.tif"))

r_ens_thresh <- r_ens
values(r_ens_thresh)[values(r_ens_bin)<1] <- NA
mask <- raster("Data/Mask/IBERIA_MASK_10km.tif")

r_ens_thresh <- projectRaster(r_ens_thresh, mask)

# Select top 15% as target areas
top15_suit <- quantile(r_ens_thresh, probs=0.85)
plot(r_ens_thresh  >= quantile(r_ens_thresh, probs=0.85))

r_suit_top15 <- r_ens_thresh >= top15_suit
plot(r_suit_top15)
values(r_suit_top15)[values(r_suit_top15)<1] <- NA
poly_suit_top15 <- rasterToPolygons(r_suit_top15)
poly_suit_top15 <-gUnaryUnion(poly_suit_top15)
poly_suit_top15 <- disaggregate(poly_suit_top15)
poly_suit_top15 <- st_as_sf(poly_suit_top15)
poly_suit_top15 <- st_transform(poly_suit_top15, crs= crs_project)

st_write(poly_suit_top15, "top15_suitability.shp", overwrite=TRUE)  

# Suitability area inside the network 2000
network2000 <- st_read("Data/rednatura2000/Natura2000_end2019_epsg3035.shp")
network2000 = st_transform(network2000, crs = crs_project)
network2000_mask <- st_crop(network2000, poly_suit_top15 )
#exclude A: SPAs (Special Protection Areas - sites designated under the Birds Directive);
network2000_mask <- network2000_mask %>% 
  filter(!SITETYPE== "A")
network2000_mask$area <- st_area(network2000_mask) #Take care of units
network2000_mask$area1 <- as.numeric(network2000_mask$area) #Take care of units
plot(network2000_mask$geometry)
pi <- st_intersection(poly_suit_top15, network2000_mask)
plot(st_geometry(pi))
pi$area <- st_area(pi) #Take care of units
poly_suit_top15$area <- st_area(poly_suit_top15)
plot(st_geometry(poly_suit_top15))
area_inside_conservation <- sum(pi$area)
area_total <- sum(poly_suit_top15$area)
suitability15_area <- area_inside_conservation/area_total

#################
# Amount of the current distribution protected by Natura 2000

current_dist <- raster("Data/RedList/lynxpardinus_rasterized_10km.tif")
plot(current_dist)
poly_current_dist <- rasterToPolygons(current_dist)
poly_current_dist <-gUnaryUnion(poly_current_dist)
poly_current_dist <- disaggregate(poly_current_dist)
poly_current_dist <- st_as_sf(poly_current_dist)
poly_current_dist = st_transform(poly_current_dist, crs = crs_project)
pi_dist <- st_intersection(poly_current_dist, network2000_mask)
pi_dist$area <- st_area(pi_dist) #Take care of units
poly_current_dist$area <- st_area(poly_current_dist)
plot(st_geometry(poly_current_dist))
#st_write(poly_suit_top15, "R/new_analysis_05_2021/results_bear/poly_suit_top15.shp", overwrite=TRUE)  
area_inside_conservation <- sum(pi_dist$area)
area_total <- sum(poly_current_dist$area)
dist_area <- area_inside_conservation/area_total


# ----------------------------------------------------------------------------------------------------
# 2. PRIORITISE areas that will likely be colonized under the low dispersal scenarios ####
# ----------------------------------------------------------------------------------------------------
# Assess how much of this is adjacent or closeby protected areas 
# Identify areas that are potentially large enough for hosting populations if colonized ( candidate sites for new protected areas?)
r_disp_025 <- raster(paste0("Data/Dispersal/", species,"_th0.025/", species, "_th0.025.tif"))
r_disp_05 <- raster(paste0("Data/Dispersal/", species, "_th0.05/", species, "_th0.05.tif"))
r_disp_1 <- raster(paste0("Data/Dispersal/", species, "_th0.1/", species, "_th0.1.tif"))
r_disp_25 <- raster(paste0("Data/Dispersal/", species, "_th0.25/", species, "_th0.25.tif"))

values(r_disp_025)[values(r_disp_025) < -1 & !is.na(values(r_disp_025))] <- NA
values(r_disp_05)[values(r_disp_05) < -1 & !is.na(values(r_disp_05))] <- NA
values(r_disp_1)[values(r_disp_1) < -1 & !is.na(values(r_disp_1))] <- NA
values(r_disp_25)[values(r_disp_25) < -1 & !is.na(values(r_disp_25))] <- NA
# Work with the dispersal scenario in low human activity
top15_low_disp <- quantile(r_disp_025, probs=0.85)
plot(r_disp_025  >= quantile(r_disp_025, probs=0.85))

r_disp_025_top15 <- r_disp_025 >= top15_low_disp

values(r_disp_025_top15)[values(r_disp_025_top15)<1] <- NA
plot(r_disp_025_top15)
poly_disp_025_top15 <- rasterToPolygons(r_disp_025_top15)
poly_disp_025_top15 <-gUnaryUnion(poly_disp_025_top15)
poly_disp_025_top15 <- disaggregate(poly_disp_025_top15)
poly_disp_025_top15 <- st_as_sf(poly_disp_025_top15)
poly_disp_025_top15 = st_transform(poly_disp_025_top15, crs_project)
plot(st_geometry(poly_disp_025_top15))
st_write(poly_disp_025_top15, "poly_disp_025_top15.shp", overwrite=TRUE)  
# How much of the dispersal area under the low human pressure scenario are Protected
network2000 = st_transform(network2000, crs_project)
network2000_mask <- st_crop(network2000, poly_disp_025_top15 )
pi_dist <- st_intersection(poly_disp_025_top15, network2000_mask)
pi_dist$area <- st_area(pi_dist) #Take care of units
poly_disp_025_top15$area <- st_area(poly_disp_025_top15)
plot(st_geometry(poly_disp_025_top15))
area_inside_conservation <- sum(pi_dist$area)
area_total <- sum(poly_disp_025_top15$area)
dispersal_area <- area_inside_conservation/area_total

#### a.	Assess how much of this is adjacent or closeby protected areas#####

# Minimum distance between the dispersal area and the conservation areas ####
# create an index of the nearest feature
index <- st_nearest_feature(x =network2000_mask, y = poly_disp_025_top15)

# slice based on the index
poly_disp_025_top15 <- poly_disp_025_top15 %>% slice(index)

plot(st_geometry(poly_disp_025_top15))

# calculate distance between polygons
poly_dist <- st_distance(x = network2000_mask, y= poly_disp_025_top15, by_element = TRUE)

# add the distance calculations to the conservation polygons
network2000_mask$distance_lowdisp <- as.vector(poly_dist)
poly_disp_025_top15$distance_lowdisp <- as.vector(poly_dist)
# Filter out PA that are over 50 Km and overlapping with PA
poly_disp_025_top15_filt <- poly_disp_025_top15 %>% 
  filter(! distance_lowdisp > 50000)
plot(st_geometry(poly_disp_025_top15_filt))

# Plot
world <- ne_countries(scale = "medium", returnclass = "sf")
#world <- st_as_sf(rnaturalearth::countries110)
IP <- dplyr::filter(world, name=="Spain" | name== "Portugal")
# A bounding box for continental Europe.
ip.bbox <- st_polygon(list(
  matrix(c(-10,35,45,35,45,45,-10,45,-10,35),byrow = T,ncol = 2)))

ip.clipped <- suppressWarnings(st_intersection(IP, st_sfc(ip.bbox, crs=st_crs(IP))))

ggplot(ip.clipped) +
  geom_sf(alpha=0.8,col='grey') +
  geom_sf(data = network2000_mask) +
  geom_sf(data = poly_disp_025_top15_filt, colour = "red", fill = NA) + 
  theme_bw()
poly_disp_025_top15_filt <- st_transform(poly_disp_025_top15_filt, crs= crs_project)

st_write(poly_disp_025_top15_filt, "poly_disp_025_top15_filter_distance.shp", overwrite=TRUE)  

#### b.	Identify areas that are potentially large enough for hosting populations if colonised (candidate sites for new protected areas?)

poly_disp_025_top15$area <- st_area(poly_disp_025_top15)
poly_disp_025_top15 <- poly_disp_025_top15 %>% 
  mutate(area2= as.numeric(area))
plot(st_geometry(poly_disp_025_top15))
plot(st_geometry(poly_disp_025_top15 %>% filter(area2 > 650000000)))

dispersal_low_big <- poly_disp_025_top15 %>% 
  filter(! area2 < 650000000)

st_write(dispersal_low_big, "dispersal_low_big.shp", overwrite=TRUE)  

## Calculate differences between dispersal low human pressure scenario and PA

#' st_erase
#'
#' @description Removes intersection from two spacial objects. This is useful when removing water areas from TIGER county files
#' @param x sf dataframe
#' @param y sf datafram
#'
#' @return
#' @export
#'
#'
st_erase <- function(x, y) {
  #Check coordinate reference system
  if(st_crs(x) != st_crs(y)){
    y <- st_transform(y, st_crs(x))
  }
  sf::st_difference(x, st_union(y))
}


diff_PA_dispersal <- st_erase(poly_disp_025_top15, network2000_mask)
plot(st_geometry(diff_PA_dispersal))
# Select areas larger enough to keep carnivores
diff_PA_dispersal_large <- diff_PA_dispersal
diff_PA_dispersal_large$area <- st_area(diff_PA_dispersal_large)
diff_PA_dispersal_large <- diff_PA_dispersal_large %>% 
  mutate(area2= as.numeric(area))
plot(st_geometry(diff_PA_dispersal_large %>% filter(area2 > 650000000)))

diff_PA_dispersal_large_filt <- diff_PA_dispersal_large %>% 
  filter(! area2 < 650000000)
st_write(diff_PA_dispersal_large_filt, "dispersal_big_noPA_lynx.shp", overwrite=TRUE)  


# ----------------------------------------------------------------------------------------------------
# 3. PRIORITISE areas that will likely be colonized under high human pressure ####
# ----------------------------------------------------------------------------------------------------
#Select areas that are colonized under the high human pressure and identify 
#Which of these cells fall into areas with high human density

plot(r_disp_25)
top15_high_disp <- quantile(r_disp_25, probs=0.85)
plot(r_disp_25  >= quantile(r_disp_25, probs=0.85))

r_disp_25_top15 <- r_disp_25 >= top15_high_disp

values(r_disp_25_top15)[values(r_disp_25_top15)<1] <- NA

r_disp_25_top15 <- projectRaster(r_disp_25_top15,
                                 crs = crs_project)

# intersect human density and top 15% of the high human pressure dispersal scenario
human_density <- raster("Data/human_density/pop_density_europe/popu01clcv5.tif")

r_disp_25_top15_pro <- projectRaster(r_disp_25_top15,
                                     crs = crs(human_density))
human_density <- crop(human_density, r_disp_25_top15_pro)
human_density<-resample(human_density, r_disp_25_top15_pro, method="bilinear")

top70_high_hd <- quantile(human_density, probs=0.70)
plot(human_density  >= quantile(human_density, probs=0.70))
human_density_top <- human_density >= top70_high_hd
values(human_density_top)[values(human_density_top)<1] <- NA

# Calculate overlap
overlap_f <- mask(r_disp_25_top15_pro, human_density_top)
overlap_f <- projectRaster(overlap_f,
                           crs = crs_project)
human_density_top<- projectRaster(human_density_top,
                                  crs = crs_project)

# save results
writeRaster(r_disp_25_top15, "r_disp_25_top15.asc", overwrite=T)
writeRaster(overlap_f, "overlap_density_disp25_bis.asc", overwrite=T)
writeRaster(human_density_top, "human_density_top.asc", overwrite=T)


# ----------------------------------------------------------------------------------------------------
# 4. Identify areas that will likely be colonized under low human pressure but no high ####
# ----------------------------------------------------------------------------------------------------
###### Use erase to see what areas are colonized by low human disturbance and not in high human disturbance
# Low human pressure
# Select the top 15%
top_disp_025 <- quantile(r_disp_025, probs=0.85)
r_top_disp_025 <- r_disp_025 >= top_disp_025
values(r_top_disp_025)[values(r_top_disp_025)<1] <- NA
plot(r_top_disp_025)
poly_disp_025 <- rasterToPolygons(r_top_disp_025)
poly_disp_025 <-gUnaryUnion(poly_disp_025)
poly_disp_025 <- disaggregate(poly_disp_025)
poly_disp_025 <- st_as_sf(poly_disp_025)
poly_disp_025 = st_transform(poly_disp_025, crs_project)
plot(st_geometry(poly_disp_025))
st_write(poly_disp_025, "poly_disp_025.shp", overwrite=TRUE)  

#High human pressure
top_disp_25 <- quantile(r_disp_25, probs=0.85)
r_top_disp_25 <- r_disp_25 >= top_disp_25
values(r_top_disp_25)[values(r_top_disp_25)<1] <- NA
plot(r_top_disp_25)
poly_disp_25 <- rasterToPolygons(r_top_disp_25)
poly_disp_25 <-gUnaryUnion(poly_disp_25)
poly_disp_25 <- disaggregate(poly_disp_25)
poly_disp_25 <- st_as_sf(poly_disp_25)
poly_disp_25 = st_transform(poly_disp_25, crs_project)
plot(st_geometry(poly_disp_25))
st_write(poly_disp_25, "poly_disp_25.shp", overwrite=TRUE)  

# How much is protected the dispersal area by high human pressure?
pi_dist <- st_intersection(poly_disp_25, network2000_mask)
pi_dist$area <- st_area(pi_dist) #Take care of units
poly_disp_25$area <- st_area(poly_disp_25)
plot(st_geometry(poly_disp_25))
area_inside_conservation <- sum(pi_dist$area)
area_total <- sum(poly_disp_25$area)
dispersal_h_area <- area_inside_conservation/area_total

# Areas that will be colonized only by low human pressure
low_colonize <- st_erase(poly_disp_025, poly_disp_25)
plot(st_geometry(low_colonize))
st_write(low_colonize, "low_colonize.shp", overwrite=TRUE)  
