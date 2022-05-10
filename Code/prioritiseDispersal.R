# ====================================================================================================
#
# IBERIAN CONSERVATION
# Prioritization
#
# ====================================================================================================


# ----------------------------------------------------------------------------------------------------
# IMPORT
# ----------------------------------------------------------------------------------------------------
library(raster)

#setwd('data/Projects/Emmy_BIOPIC/WPX_restoration/Pratzer_Nill_IberianConservation')
setwd("~/Documents/UNI/Master/3.Semester/GCIB/Publishing/processing/iberian_conservation")

species <- "ursusarctos"
species <- "lynxpardinus"

# ----------------------------------------------------------------------------------------------------
# PRIORITISE BASED ON SDM
# ----------------------------------------------------------------------------------------------------

r_ens <- raster(paste0("Data/SDMs/", species, "_SDM_ensemble_mn_md_wmn_cav_sdp.tif"))
r_ens_bin <- raster(paste0("Data/SDMs/", species, "_SDM_ensembles_binary_mn_md_wmn.tif"))

r_ens_thresh <- r_ens
values(r_ens_thresh)[values(r_ens_bin)<1] <- NA
r_ens_thresh <- projectRaster(r_ens_thresh, mask)

# Select top 15% as target areas
plot(r_ens_thresh  >= quantile(r_ens_thresh, probs=0.85))


# ----------------------------------------------------------------------------------------------------
# PRIORITISE BASED ON DISPERSAL
# ----------------------------------------------------------------------------------------------------

r_disp_025 <- raster(paste0("Data/Dispersal/", species,"_th0.025/", species, "_th0.025.tif"))
r_disp_05 <- raster(paste0("Data/Dispersal/", species, "_th0.05/", species, "_th0.05.tif"))
r_disp_1 <- raster(paste0("Data/Dispersal/", species, "_th0.1/", species, "_th0.1.tif"))
r_disp_25 <- raster(paste0("Data/Dispersal/", species, "_th0.25/", species, "_th0.25.tif"))

values(r_disp_025)[values(r_disp_025) < -1 & !is.na(values(r_disp_025))] <- NA
values(r_disp_05)[values(r_disp_05) < -1 & !is.na(values(r_disp_05))] <- NA
values(r_disp_1)[values(r_disp_1) < -1 & !is.na(values(r_disp_1))] <- NA
values(r_disp_25)[values(r_disp_25) < -1 & !is.na(values(r_disp_25))] <- NA

r_disp <- stack(r_disp_025, r_disp_05, r_disp_1, r_disp_25)

r_disp_mean <-  calc(r_disp, fun=mean, na.rm=T)
r_disp_sd <-  calc(r_disp, fun=sd, na.rm=T)

# Select top 15% as target areas
plot(r_disp_mean  >= quantile(r_disp_mean, probs=0.85))
plot(r_disp_sd  >= quantile(r_disp_sd, probs=0.85))

