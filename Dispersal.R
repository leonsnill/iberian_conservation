# ====================================================================================================
#
# 3 GLOBAL CHANGE IMPACTS ON BIODIVERSITY:
# Dispersal
#
# ====================================================================================================

# ----------------------------------------------------------------------------------------------------
# IMPORT and DEFINE
# ----------------------------------------------------------------------------------------------------
library(raster)
library(dplyr)
library(tidyverse)
library(MigClim)

setwd("/Users/leonnill/Documents/MSc_GCG/MSc_GCIB")

sdm <- brick("SDMs/SDM_ensembles_mn_md_wmn_cav_sdp.tif")
sdm <- sdm[[3]] * 1000
sdm <- as.integer(sdm)
sdm[is.na(sdm)] <- 0

suit <- sdm
suit[suit < 0.52*1000] <- 0

init_dist <- raster("GBIF/GBIF_ursusarctos_europe_1990-2020_presence_10km.tif")
init_dist <- as.data.frame(init_dist, xy=TRUE)
init_dist <- drop_na(init_dist)
init_dist <- rasterFromXYZ(init_dist, crs=crs(sdm))
init_dist[is.na(init_dist)] <- 0


# mask
lcf_a <- raster("LC_Fractions/europe_lc_fraction_artificial_10km.tif")/1000
lcf_h <- raster("LC_Fractions/europe_lc_fraction_highveg_10km.tif")/1000
lcf_l <- raster("LC_Fractions/europe_lc_fraction_lowveg_10km.tif")/1000
motorways <- raster("Watermask/EUROPE_MASK_10km_MOTORWAYS.tif")

# Choose mask
mask <- Which(lcf_a > 0.05 | motorways == 2)  # Scenario Artficial + Motorways
mask <- Which(lcf_a > 0.05)  # Scenario Artifical

mask <- as.data.frame(mask, xy=TRUE)
mask <- drop_na(mask)
mask <- rasterFromXYZ(mask, crs=crs(sdm))
mask[is.na(mask)] <- 0
mask <- projectRaster(mask, init_dist, method="ngb")


# kernel (in cell units)
dispersal <- 2
distances <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)  # given the resolution of 10x10km -> e.g. 1 = 10km
dispn <- exp(-(distances/dispersal))


# MigClim simulation
MigClim.migrate(iniDist = as.data.frame(init_dist, xy = TRUE),
                hsMap = as.data.frame(suit),
                rcThreshold = 0,
                envChgSteps=1,
                dispSteps=10,
                replicateNb=1,
                iniMatAge = 4,  # sexual maturity of females between 4-8 years
                propaguleProd = c(0.7, 0.7, 0.7, 0.5, 0.5),
                barrier = as.data.frame(mask),
                barrierType = "weak",
                dispKernel = dispn,
                simulName="ursus_arctos_disp",
                overWrite=T)


# Map colonisation ability
res <- raster("ursus_arctos_disp/ursus_arctos_disp_raster.asc")
res[res == 30000] <- 111
res[res < 0] <- -1
res[res == 1] <- 1
res[res > 1] <- res[res > 1]-99
res[res == 0] <- NA
crs(res) <- crs(lcf_a)

# crop to Iberian Peninsula
iberia <- raster("Watermask/IBERIA_MASK_10km.tif")
res_crop <- crop(res, iberia)
res_crop <- res_crop-1
res_crop[res_crop == -2] <- -1

writeRaster(res_crop, "dispersal_ursus_arctos_10a_onlylandbarrier.tif", format="GTiff", overwrite = TRUE)

