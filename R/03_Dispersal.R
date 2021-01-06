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

#species <- "ursusarctos"
species <- "lynxpardinus"
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
sel_th <- 0.05
bar_mask <- as.integer(Which(artificial > sel_th))  # Scenario (0.025, 0.05, 0.1, 0.25)
bar_mask[is.na(bar_mask)] <- 0


# ----------------------------------------------------------------------------------------------------
# RUN MODEL
# ----------------------------------------------------------------------------------------------------
setwd(paste0(getwd(), "/Data/Dispersal"))
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
                simulName=paste0(species, "_th", toString(sel_th)),
                overWrite=T)


# ----------------------------------------------------------------------------------------------------
# CAST TO RASTER
# ----------------------------------------------------------------------------------------------------
# Map colonisation ability
res <- raster(paste0(species, "_th", toString(sel_th), "/",species, "_th", toString(sel_th), "_raster.asc"))
# reclassify
rules <- c(0,-999, #cell has never been occupied
           -101,-1, #cell decolonized
           1,0, #stable: cells that were initially occupied and remained occupied
           101,1, #cell colonized in year 1
           102,2, #cell colonized in year 2
           103,3, #cell colonized in year 3
           104,4, #cell colonized in year 4
           105,5, #cell colonized in year 5
           106,6, #cell colonized in year 5
           107,7, #cell colonized in year 5
           108,8, #cell colonized in year 5
           109,9, #cell colonized in year 5
           110,10, #cell colonized in year 5
           30000,20) #potentially suitable cells that remained unoccupied
rules <- matrix(rules, ncol=2, byrow=T)
res <- reclassify(res, rules)
crs(res) <- coordsys
res[(res != 0) & (init_dist == 1)] <- -1
res <- mask(res, mask, maskvalue=0)

writeRaster(res, paste0(species, "_th", toString(sel_th), "/",species, "_th", toString(sel_th), ".tif"), 
            format="GTiff", overwrite = TRUE)


# ---- Plot
cls1 <- c("grey60", 'darkred',
         heat.colors(10), 'pink')
legend_colon1 <- c("Never Occupied", "Initial",
                  "1 year", "2 years", "3 years", "4 years", "5 years",
                  "6 years", "7 years", "8 years", "9 years", "10 years",
                  "Suitable, but unoccupied")
cls2 <- c("grey60", "blue", 'darkred',
          heat.colors(10), 'pink')
legend_colon2 <- c("Never Occupied", "Decolonised", "Initial",
                   "1 year", "2 years", "3 years", "4 years", "5 years",
                   "6 years", "7 years", "8 years", "9 years", "10 years",
                   "Suitable, but unoccupied")

cls <- cls2
legend_colon <- legend_colon2
res <- as.factor(res)
rat <- levels(res)[[1]]
rat[["Colonisation"]] <- legend_colon
levels(res) <- rat
library(rasterVis)
png(file=paste0(species, "_th", toString(sel_th), "/", "plot_th", toString(sel_th), ".png"),
    width=1000, height=600)
levelplot(res, margin=F, scales=list(draw=FALSE), col.regions=cls)
dev.off()

