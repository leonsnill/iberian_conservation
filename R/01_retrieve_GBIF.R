# ====================================================================================================
#
# 1 GLOBAL CHANGE IMPACTS ON BIODIVERSITY:
# GBIF Data Download and Pre-Processing
#
# ====================================================================================================

# with some code snippets borrowed from https://damariszurell.github.io/HU-GCIB/

# ----------------------------------------------------------------------------------------------------
# IMPORT and DEFINE
# ----------------------------------------------------------------------------------------------------
library(rgbif)
library(raster)
library(dplyr)
library(tidyr)
library(rgdal)
library(CoordinateCleaner)

setwd("C:/Users/Leon/Google Drive/03_LSNRS/Projects/Iberian_Conservation/iberian_conservation")

species_name <- "Ursus arctos"
scientificName <- "Ursus arctos Linnaeus, 1758"
output_name <- "ursusarctos"
roi <- 'POLYGON((-9.843 35.679, 32.344 35.679, 32.344 71.414, -9.843 71.414, -9.843 35.679))'

species_name <- "Lynx pardinus"
scientificName <- "Lynx pardinus (Temminck, 1827)"
output_name <- "lynxpardinus"
roi <- 'POLYGON((-9.747 35.844, 3.392 35.844, 3.392 43.858, -9.747 43.858, -9.747 35.844))'

# retrieve data from GBIF database based on certain conditions
df <- occ_search(scientificName = species_name, return='data', limit=100000,
                        hasCoordinate = TRUE,
                        basisOfRecord = "HUMAN_OBSERVATION",
                        geometry = roi,
                        year = '1980, 2020')
df <- df$data

# write to .csv
#write.csv(df, paste0("Data/GBIF/", output_name, "_europe_1980-2020.csv"))

# read 
df <- read_csv(paste0("Data/GBIF/", output_name, "_europe_1980-2020.csv"))

output_name <- "lynxpardinus15"

# ----------------------------------------------------------------------------------------------------
# VISUALIZE DATA
# ----------------------------------------------------------------------------------------------------
library(maptools)

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-9.8, 30), ylim=c(35, 70))
points(df$decimalLongitude, df$decimalLatitude, col='red',  pch=4)


# ----------------------------------------------------------------------------------------------------
# CLEAN DATA
# ----------------------------------------------------------------------------------------------------
# FILTER GENUS
df <- df[df$scientificName == scientificName,]

# COORDINATES / DUPLICATES
# remove potential NAs (should alread be taken care of by "hasCoordinate = TRUE" in occ_search)
df <- drop_na(df, c("decimalLatitude", "decimalLongitude"))

df_clean <- clean_coordinates(df, lon="decimalLongitude", lat="decimalLatitude",
                                     countries="countryCode",
                                     tests=c("centroids","outliers"),
                                     value = "clean")
df_clean <- cc_dupl(df_clean, lon="decimalLongitude", lat="decimalLatitude")

# coordinate uncertainty
hist(df_clean$coordinateUncertaintyInMeters/1000, breaks = 20)
df_clean <- df_clean %>%
  filter(coordinateUncertaintyInMeters/1000 <= 10 | is.na(coordinateUncertaintyInMeters))

# PRESENCE / ABSENCE
table(df_clean$individualCount)
# 0     1     2     3     4 
# 60 12463    27    17    11 
df_clean <- df_clean %>% 
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 10 | is.na(individualCount))

table(df_clean$taxonRank)

presence_points <- SpatialPoints(coords = df_clean[,c("decimalLongitude", "decimalLatitude")], proj4string = CRS("+proj=longlat +datum=WGS84"))

# ----------------------------------------------------------------------------------------------------
# IUCN RANGE MAPS
# ----------------------------------------------------------------------------------------------------

# NEW: not filtered by IUCN range maps to consider the entire niche, not only the realised one

#iucn <- readOGR("IUCN/iucn_ursus_arctos/data_0.shp")
#plot(wrld_simpl, xlim=c(-9.5, 32.5), ylim=c(35.5,72))
#plot(wrld_simpl, xlim=c(-9.8, 4.5), ylim=c(35, 45))
#plot(iucn[iucn$PRESENCE == "1",], col='blue', add=T)
#plot(iucn[iucn$PRESENCE == "5",], col='yellow', add=T)
#points(df_clean$decimalLongitude, df_clean$decimalLatitude, col='red',  pch=4)

# filter GBIF based on IUCN extant range
#presence_points <- SpatialPoints(coords = df_clean[,c(5,4)], proj4string = cdfrs(iucn))
#presence_points <- presence_points[iucn,]


# ----------------------------------------------------------------------------------------------------
# CREATE TEMPORARY BINARY MASK
# ----------------------------------------------------------------------------------------------------
# reproject presence points to target crs
ref_img <- raster("Data/Mask/EUROPE_MASK_10km.tif")
# adjust reference image to Iberia for the lynx
if(species_name == "Lynx pardinus"){
  ref_img <- raster("Data/Mask/IBERIA_MASK_10km.tif")
}
presence_points <- spTransform(presence_points, crs(ref_img))

# rasterize points and create binary mask of presence
pres_count <- rasterize(presence_points, ref_img, fun='count')
writeRaster(pres_count, paste0("Data/GBIF/", output_name, "_europe_1990-2020_count.tif"), driver="GTiff", overwrite = TRUE)
pres_bin <- Which(pres_count > 0)
pres_bin <- as.integer(pres_bin)
names(pres_bin) <- "presence"

# mask by water mask and write to disc
pres_bin <- mask(pres_bin, ref_img, maskvalue=1, inverse=TRUE)
writeRaster(pres_bin, paste0("Data/GBIF/", output_name, "_europe_1990-2020_presence.tif"), driver="GTiff", overwrite = TRUE)

# 
df_pres <- as.data.frame(pres_bin, xy = TRUE)
dpdf_pres <- SpatialPointsDataFrame(coords = df_pres[,c(1,2)], data = df_pres, proj4string = crs(pres_bin))
dpdf_pres <- spTransform(dpdf_pres, CRS("+proj=longlat +datum=WGS84"))  # to lat lon
df_pres <- as.data.frame(dpdf_pres)  # to df
colnames(df_pres) <- c("x", "y", "presence", "lon", "lat")


# ----------------------------------------------------------------------------------------------------
# THIN DATASET AND CREATE PSEUDO ABSENCES
# ----------------------------------------------------------------------------------------------------
library(spThin)

# thin presence (1) & absence (0)
l_thinned <- list()

for (i in c(0, 1)){
  pres_temp <- df_pres %>%
    dplyr::select(presence, lon, lat) %>%
    filter(presence == i)
  
  if (i == 0) {
    if (species_name == "Ursus arctos") {
      pres_temp <- sample_n(pres_temp, 20000)
    }
    pres_temp$presence <- 1
  }
  
  temp_thin <- thin(pres_temp,
                    lat.col='lat',
                    long.col='lon',
                    spec.col='presence',
                    thin.par=11,
                    reps=5,
                    write.files=F,
                    locs.thinned.list.return=T)
  
  # maximum presence records
  max_pres <- which.max(sapply(temp_thin, nrow))
  temp_thin <- temp_thin[[max_pres]]
  
  # merge
  temp_thin <- merge(pres_temp, data.frame(lon=temp_thin$Longitude, lat=temp_thin$Latitude),
                     by=c('lon', 'lat'))  # to presence df
  
  temp_thin <- merge(df_pres[,c(1,2,4,5)], temp_thin, by = c("lon", "lat"), all = TRUE)
  #temp_thin[is.na(temp_thin)] <- 0
  temp_thin <- dplyr::select(temp_thin, x, y, presence)
  l_thinned[[i+1]] <- rasterFromXYZ(temp_thin, res = c(10000, 10000), crs=crs(ref_img))
}

img_absence <- l_thinned[[1]]
img_presence <- l_thinned[[2]]

# write to disc
writeRaster(img_presence, paste0("Data/GBIF/", output_name, "_europe_1990-2020_presence_thinned.tif"), driver="GTiff", overwrite = TRUE)
writeRaster(img_absence, paste0("Data/GBIF/", output_name, "_europe_1990-2020_absence_thinned.tif"), driver="GTiff", overwrite = TRUE)


# ----------------------------------------------------------------------------------------------------
# CREATE BUFFER AROUND PRESENCES
# ----------------------------------------------------------------------------------------------------
buffer_presence <- raster::buffer(img_presence, width=50000)
writeRaster(buffer_presence, paste0("Data/GBIF/", output_name, "_europe_1990-2020_presence_thinned_buffer_50km.tif"), driver="GTiff", overwrite = TRUE)




