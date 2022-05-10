# ====================================================================================================
#
# IBERIAN CONSERVATION
# Spatial priorizitation
# Guillermo Fandos
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

#remotes::install_github("michaeldorman/geobgu")

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(rgeos) #OPERACIONES GEOMÉTRICAS CON INFO GEOGRÁFICA
library(dismo) #LIBRERIA PARA MODELOS DE DISTRIBUCION
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)
library(viridis)

source("scripts/st_erase.R") # Function to delete the intersection between two shapefiles
source("scripts/regional_seas.R") # Create buffer divided by closest region
species <- "ursusarctos"
#species <- "lynxpardinus"

# CRS of the project
crs_project <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs" # CRS EPSG:25830

# We are going to follow four different strategies to spatial prioritirize the conservation of both species
#1. Identify top-ranking areas (15% percent) with highest suitability per species, regardless of whether species can get there or not.
#a.	Assess how much of this is already in protected areas
#2.	Select those areas that will likely be colonized under the “vague” dispersal scenarios (human pressure provide low dispersal cost) 25% Threshold (few human barriers acting as dispersal limitation). Only big cities act as a barrier. 
#a.	Assess how much of this is adjacent or closeby protected areas (candidate sites for protected area expansion or increasing protection level)
#b.	Identify areas that are potentially large enough for hosting populations if colonized ( candidate sites for new protected areas?
#3.	Select areas that are colonized under the restrictive dispersal scenario ( 2.5% scenario where small settlements and roads act as dispersal barrier). Large proportion of barriers and very high cost for crossing the landscape. 
#These are areas where conflicts will become likely in the future (or already are)
#4.	Select areas that are colonized under vague dispersal (25%) scenario but not under restrictive dispersal scenario (2.5%) 
#a.	These are areas where human pressure, and likely conflict, are higih and prevent range expansion – two possible pathways from here
#steer animals away by making other, less conflict prone ares more attractive
#work with local people to increase tolerance (particularly good if this is a critical bottleneck)
#Reduce dispersal limitation by wildlife road bridges

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

st_write(poly_suit_top15, "results/bear/bear_top15_suitability.shp", overwrite=TRUE)  

# Suitability area inside the network 2000
network2000 <- st_read("Data/rednatura2000/Natura2000_end2019_epsg3035.shp")
network2000 = st_simplify(network2000, dTolerance = 1000)  # 1000 m
network2000 <- st_as_sf(network2000)
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
# The current distribution
current_dist <- raster("Data/RedList/ursusarctos_rasterized_10km.tif")
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

# -------------------------------------------------------------------------------------------------------------------
# 2. PRIORITISE areas that will likely be colonized under the the “vague” dispersal scenarios 25% Threshold ######## 
#(human pressure provide low dispersal cost)   ####
# -------------------------------------------------------------------------------------------------------------------
# Assess how much of this is adjacent or closeby protected areas 
# Identify areas that are potentially large enough for hosting populations if colonized ( candidate sites for new protected areas?)
r_disp_restrictive <- raster(paste0("Data/Dispersal/", species,"_th0.025/", species, "_th0.025.tif"))
r_disp_vague <- raster(paste0("Data/Dispersal/", species, "_th0.25/", species, "_th0.25.tif"))

# These raster are classified with:
#-999 cell never occupied
#-1:  cell decolonized 
# 0: stable: cells that were initially occupied and remained occupied
# #1-10: colonization
# 20: potentially suitable cells that remained unoccupied
# We select only the colonize cells. Because we want to see what affect this colonization

values(r_disp_restrictive)[values(r_disp_restrictive) < 1 & !is.na(values(r_disp_restrictive))] <- NA
values(r_disp_vague)[values(r_disp_vague) < 1 & !is.na(values(r_disp_vague))] <- NA
values(r_disp_restrictive)[values(r_disp_restrictive) == 20] <- NA
values(r_disp_vague)[values(r_disp_vague) == 20] <- NA
values(r_disp_restrictive)[values(r_disp_restrictive) < 11] <- 1
values(r_disp_vague)[values(r_disp_vague) < 11] <- 1

plot(r_disp_vague)

# Transform to polygon. 
poly_disp_vague <- rasterToPolygons(r_disp_vague)
poly_disp_vague <-gUnaryUnion(poly_disp_vague)
poly_disp_vague <- disaggregate(poly_disp_vague)
poly_disp_vague <- st_as_sf(poly_disp_vague)
poly_disp_vague = st_transform(poly_disp_vague, crs_project)
plot(st_geometry(poly_disp_vague))
st_write(poly_disp_vague, "results/bear/poly_disp_vague.shp", overwrite=TRUE)  
##################
# How much of the dispersal area under the low human pressure scenario are Protected
network2000 = st_transform(network2000, crs_project)
network2000_mask <- st_crop(network2000, poly_disp_vague )
st_write(network2000_mask, "results/bear/network_2000_mask.shp", overwrite=TRUE)  
pi_dist <- st_intersection(poly_disp_vague, network2000_mask)
pi_dist$area <- st_area(pi_dist) #Take care of units
poly_disp_vague$area <- st_area(poly_disp_vague)
plot(st_geometry(poly_disp_vague))
area_inside_conservation <- sum(pi_dist$area)
area_total <- sum(poly_disp_vague$area)
dispersal_area_vague <- area_inside_conservation/area_total

#### a.	Assess how much of this is adjacent or closeby protected areas#####
#network2000_mask_s = st_simplify(network2000_mask, dTolerance = 1000)  # 2000 m
network2000_mask_s <- st_as_sf(network2000_mask)
plot(st_geometry(network2000_mask_s))
ggplot() +
  geom_sf(data = network2000_mask_s, colour = "red", fill = "cyan")
diff_PA_dispersal <- st_erase(poly_disp_vague, network2000_mask_s)
plot(st_geometry(diff_PA_dispersal))

ggplot() +
  geom_sf(data = diff_PA_dispersal, colour = "red", fill = "cyan")

st_write(diff_PA_dispersal, "results/bear/dispersal_vague_noPA.shp", overwrite=TRUE)  


# Buffer around Network 2000


network2000_mask_buffer <- regional_seas(
  x = network2000_mask_s,
  group = "SITECODE",
  dist = units::set_units(50, km), # buffer distance
  density = units::set_units(0.5, 1/km) # density of points (the higher, the more precise the region attribution)
)


#ggplot() +
 # geom_sf(data = network2000_mask_buffer,
  #        aes(colour = NOM_DEPT, fill = NOM_DEPT),
   #       alpha = 0.25) +
  #geom_sf(data = network2000_mask_s,
   #       aes(fill = NOM_DEPT),
    #      colour = "grey20",
     #     alpha = 0.5) +
  #scale_fill_viridis_d() +
  #scale_color_viridis_d() +
  #theme_bw() + theme(legend.position="none")

network2000_mask_buffer_s <-  st_simplify(network2000_mask_buffer, dTolerance = 1000)  # 2000 m

diff_buffer <- st_erase(network2000_mask_buffer_s, network2000_mask_s)

#ggplot() +
 # geom_sf(data = diff_buffer, colour = "red", fill = "cyan")

poly_disp_vague_s <-  st_simplify(poly_disp_vague, dTolerance = 1000)  # 2000 m

int_buffer <- st_intersection(poly_disp_vague_s, diff_buffer)
#ggplot() +
 # geom_sf(data = int_buffer, colour = "red", fill = "cyan")

int_buffer_2 <- st_collection_extract(int_buffer, "POLYGON")
st_write(int_buffer_2, "results/bear/int_buffer.shp", overwrite=TRUE)  


# -----------------------------------------------------------------------------------------------------
# 3. PRIORITISE areas that will likely be colonized under the restrictive dispersal scenario ####
#2.5% scenario where small settlements and roads act as dispersal barrier
# -----------------------------------------------------------------------------------------------------
#Select areas that are colonized under the high human pressure and identify 

plot(r_disp_restrictive)
r_disp_restrictive <- projectRaster(r_disp_restrictive,
                                    crs = crs_project)
# Transform to polygon 
poly_disp_restrictive <- rasterToPolygons(r_disp_restrictive)
poly_disp_restrictive <-gUnaryUnion(poly_disp_restrictive)
poly_disp_restrictive <- disaggregate(poly_disp_restrictive)
poly_disp_restrictive <- st_as_sf(poly_disp_restrictive)
poly_disp_restrictive = st_transform(poly_disp_restrictive, crs_project)
plot(st_geometry(poly_disp_restrictive))
st_write(poly_disp_restrictive, "results/bear/poly_disp_restrictive.shp", overwrite=TRUE)  

plot(st_geometry(network2000_mask_s))

# How much is from the dispersal restrictive scenario is protected
pi_dist_restrictive <- st_intersection(poly_disp_restrictive, network2000_mask_s)
pi_dist_restrictive$area <- st_area(pi_dist_restrictive) #Take care of units
poly_disp_restrictive$area <- st_area(poly_disp_restrictive)
plot(st_geometry(poly_disp_restrictive))
area_inside_conservation <- sum(pi_dist_restrictive$area)
area_total <- sum(poly_disp_restrictive$area)
dispersal_area_restrictive <- area_inside_conservation/area_total

# ----------------------------------------------------------------------------------------------------
# 4. Identify areas that will likely be colonized under low human pressure but no high ####
# ----------------------------------------------------------------------------------------------------

low_colonize <- st_erase(poly_disp_vague, poly_disp_restrictive)
plot(st_geometry(low_colonize))
st_write(low_colonize, "results/bear/low_colonize.shp", overwrite=TRUE)  


