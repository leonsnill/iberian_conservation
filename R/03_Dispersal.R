# ====================================================================================================
#
# IBERIAN CONSERVATION
# Dispersal Models
#
# ====================================================================================================


# ----------------------------------------------------------------------------------------------------
# IMPORT
# ----------------------------------------------------------------------------------------------------
library(raster)
library(dplyr)
library(tidyverse)
library(MigClim)

setwd("C:\\Users\\Leon\\Google Drive\\03_LSNRS\\Projects\\Iberian_Conservation\\iberian_conservation")
#setwd("~/Documents/UNI/Master/3.Semester/GCIB/Publishing/processing/iberian_conservation")

species <- "ursusarctos"
#species <- "lynxpardinus"
mask <- raster("Data/Mask/IBERIA_MASK_10km.tif")
sdm <- brick(paste0("Data/SDMs/", species, "_SDM_ensemble_mn_md_wmn_cav_sdp.tif"))
sdm <- projectRaster(sdm, mask)
sdm <- sdm[[3]] * 1000  # layer 3 = weighted mean prediction
sdm <- as.integer(sdm)
sdm[is.na(sdm)] <- 0

# get coordinate system
coordsys <- crs(sdm)

# mask
artificial <- raster("Data/Mask/roads-train-artificial_iberia_fraction_10km.tif")/1000
artificial <- mask(artificial, mask, maskvalue = 0)

# Choose mask
bar_mask <- Which(artificial > 0.05)  # Scenario 

bar_mask <- as.data.frame(bar_mask, xy=TRUE)
bar_mask <- drop_na(bar_mask)
bar_mask <- rasterFromXYZ(bar_mask, crs=coordsys)
bar_mask[is.na(bar_mask)] <- 0

# ----------------------------------------------------------------------------------------------------
# DEFINE PARAMS
# ----------------------------------------------------------------------------------------------------
# SDM threshold for suitable habitat
if(species == "lynxpardinus"){
  threshold <- 0.28 # optimal threshold for lynx (weighted mean ensemble) 
} else {
  threshold <- 0.46 # optimal threshold for bear (weighted mean ensemble)
}

# kernel (in cell units)
if(species == "lynxpardinus"){
  dispersal <- 1.47 #(Ferreras et al. 2004)
} else {
  dispersal <- 2
}
distances <- 1:9  # given the resolution of 10x10km -> e.g. 1 = 10km
dispn <- exp(-(distances/dispersal))

# age of sexual maturaty
if(species == "lynxpardinus"){
  matAge <- 3  #(Palomares et al. 2005)
} else {
  matAge <- 4  # sexual maturity of females between 4-8 years
}


# ----------------------------------------------------------------------------------------------------
# PREPARE INPUT
# ----------------------------------------------------------------------------------------------------
sdm[sdm < threshold*1000] <- 0
sdm <- mask(sdm, mask, maskvalue = 0)
sdm[is.na(sdm)] <- 0  # unsuitable must be set to zero

init_dist <- raster(paste0("Data/RedList/", species, "_rasterized_10km.tif"))
init_dist[is.na(init_dist)] <- 0  # unoccupied must be set to zero

# input mask
artificial <- raster("Data/Mask/roads-train-artificial_iberia_fraction_10km.tif")/1000
artificial <- mask(artificial, mask, maskvalue = 0)

# define scenario based on mask
bar_mask <- as.integer(Which(artificial > 0.05))  # Scenario 
bar_mask[is.na(bar_mask)] <- 0


# ----------------------------------------------------------------------------------------------------
# RUN MODEL
# ----------------------------------------------------------------------------------------------------
# MigClim simulation
MigClim.migrate(iniDist = as.data.frame(init_dist, xy = TRUE),
                hsMap = as.data.frame(sdm),
                rcThreshold = 0,
                envChgSteps=1,
                dispSteps=10,  # 10 years
                replicateNb=1,
                iniMatAge = matAge,
                #propaguleProd = c(0.7, 0.7, 0.7, 0.5, 0.5),
                barrier = as.data.frame(bar_mask),
                barrierType = "weak",
                dispKernel = dispn,
                simulName=paste0(species, "_disp"),
                overWrite=T)


# ----------------------------------------------------------------------------------------------------
# CAST TO RASTER
# ----------------------------------------------------------------------------------------------------
# Map colonisation ability
res <- raster(paste0(species, "_disp/", species, "_disp_raster.asc"))
res[res == 30000] <- 111
res[res < 0] <- -1
res[res == 1] <- 1
res[res > 1] <- res[res > 1]-99
res[res == 0] <- NA
crs(res) <- coordsys

writeRaster(res, paste0("dispersal_", species, "_10a_onlylandbarrier.tif"), format="GTiff", overwrite = TRUE)












############
#init_dist <- projectRaster(init_dist, mask)
#init_dist <- mask(init_dist, mask, maskvalue = 0)
#init_dist <- as.data.frame(init_dist, xy=TRUE)
#init_dist <- drop_na(init_dist)
#init_dist <- rasterFromXYZ(init_dist, crs=coordsys)
#init_dist[is.na(init_dist)] <- 0
#init_dist[init_dist == 0] <- 0
#init_dist[init_dist != 0] <- 1

#bar_mask <- as.data.frame(bar_mask, xy=TRUE)
#bar_mask <- drop_na(bar_mask)
#bar_mask <- rasterFromXYZ(bar_mask, crs=coordsys)
#bar_mask[is.na(bar_mask)] <- 0


# crop to Iberian Peninsula
iberia <- raster("Data/Mask/IBERIA_MASK_10km.tif")
res_crop <- crop(res, iberia)
res_crop <- res_crop-1
res_crop[res_crop == -2] <- -1