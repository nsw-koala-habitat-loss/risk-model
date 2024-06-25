# THIS CODE READS IN THE REQUIRED DATA AND ORGANISES THE DATA FOR MODEL FITTING

# load libraries
library(sf)
library(sp)
library(raster)
library(exactextractr)
library(terra)
library(tidyverse)
library(spdep)
library(readxl)
library(qs)

# load functions
source("functions.R")

# load spatial property units for each KMR
SUs <- list(CC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Coast"),
            CST = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Southern_Tablelands"),
            DRP = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Darling_Riverine_Plains"),
            FW = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Far_West"),
            NC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "North_Coast"),
            NT = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northern_Tablelands"),
            NS = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northwest_Slopes"),
            R = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Riverina"),
            SC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "South_Coast"))

# add area to the spatial units attribute table
for (i in names(SUs)) {
  SUs[[i]]$Shape_Area <- as.numeric(st_area(SUs[[i]]))
  SUs[[i]]$Area = SUs[[i]]$Shape_Area / 1e4
}

# Check for Shape_Area < 0 
All_SUs <- do.call(rbind, SUs)
All_SUs[All_SUs$Shape_Area < 0,]

# save processes spatial units
saveRDS(SUs, file = "output/spatial_units/sus.rds")

# load SA1s spatial layer
# note that these are the SA1s for the whole study area, not for each KMR
SA1s_All <- st_read("input/spatial_units/sa1s.gdb", layer = "sa1s")

# dissaggregate SA1s into individual KMRs
SA1s <- list(CC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$CC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             CST = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$CST$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             DRP = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$DRP$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             FW = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$FW$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             NC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             NT = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NT$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             NS = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NS$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             R = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$R$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
             SC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$SC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),])

# save processed SA1s
saveRDS(SA1s, file = "output/spatial_units/sa1s.rds")

# get woody cover in 2011
Woody <- terra::rast("input/woody_cover/woody_nsw.tif") %>% round()

# load woody loss rasters and create raster stack
# get all tif files
WLFiles <- list.files("./input/clearing_data", pattern = "\\.tif$", full.names = TRUE)
WLStack <- terra::rast(WLFiles) %>% round()

# slpit into the different types of clearing
Ag_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 1, 0, 0))) %>% max(na.rm = TRUE) %>% crop(Woody)
In_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 0, 1, 0))) %>% max(na.rm = TRUE) %>% crop(Woody)
Fo_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 0, 0, 1))) %>% max(na.rm = TRUE) %>% crop(Woody)

# create raster stack of woody extent and loss
Stack <- terra::rast(list(aloss = Ag_Loss, iloss = In_Loss, floss = Fo_Loss, woody = Woody))

# save raster stack
saveRDS(Stack, file = "output/raster_stacks/woodyextloss.rds")

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)

# save data
saveRDS(ZStats_Woody, file = "output/data/ZStats_Woody.rds")


rm(list = setdiff(ls(all.names = TRUE), "SUs"))
gc()

# Produce Woody loss vector ----
SUs <- readRDS("output/spatial_units/sus.rds")
ZStats_Woody <- readRDS("output/data/ZStats_Woody.rds")
SUs_Woody <- list()
for (i in 1:length(names(ZStats_Woody))){
  SUs_Woody[[i]] <- bind_cols(SUs[[i]], ZStats_Woody[[i]])
}
SUs_Woody <- do.call(rbind, SUs_Woody)

saveRDS(SUs_Woody, file = "output/data/SUs_Woody.rds")

# Preprocess covariates ----
SUs <- readRDS("output/spatial_units/sus.rds")

# # combined drought indicator
# Drought <- terra::rast("input/covariates/drought.tif")
# # elevation
# Elev <- terra::rast("input/covariates/elev.tif")
# # fire history (2011-2019)
# Fire <- terra::rast("input/covariates/fire.tif")
# # forest code
# ForCode <- terra::rast("input/covariates/forest_code.tif")
# # forest tenure
# ForTen <- terra::rast("input/covariates/forest_tenure.tif")
# # forest tenure type
# ForTenType <- terra::rast("input/covariates/forest_tenure_type.tif")
# # 2016 median household income weekly
# Income <- terra::rast("input/covariates/income.tif")
# # land use (NSW landuse 2013)
# LandUse <- terra::rast("input/covariates/landuse.tif")
# # population density
# PopDen <- terra::rast("input/covariates/pop_den.tif")
# # precipitation
# Precip <- terra::rast("input/covariates/prec.tif")
# # remoteness structure (ABS)
# Remote <- terra::rast("input/covariates/remoteness.tif")
# # slope
# Slope <- terra::rast("input/covariates/slope.tif")
# # soil fertility
# SoilFert <- terra::rast("input/covariates/soil_fert.tif")
# # soil nitrogen
# SoilNit <- terra::rast("input/covariates/soil_nitrogen.tif")
# # soil type
# SoilType <- terra::rast("input/covariates/soil_type.tif")
# # temperature
# Temp <- terra::rast("input/covariates/temp.tif")
# # # Property size
# PropSize <- terra::rast("input/covariates/prop_size.tif")
# # # Property value
# PropVal <- terra::rast("input/covariates/prop_value.tif")


# Load proposed covariates based on workshops from lookup xlsx
CovLookup <- readxl::read_xlsx("Input/covariates/covariate_description.xlsx", sheet = "AllLyr")

# load covariates
Cov_file_bn <- tools::file_path_sans_ext(basename(list.files("input/covariates", pattern = "\\.tif$", full.names = TRUE)))
for(i in 1:nrow(CovLookup)){
  if (CovLookup$`Covariate Code`[i] %in% Cov_file_bn){
    print(paste0("Load ... ", CovLookup$Type[i] ," ... ", CovLookup$`Covariate Code`[i], " = ",  CovLookup$`Covariate Description`[i]))
    assign(CovLookup$`Covariate Code`[i], terra::rast(paste0("input/covariates/", CovLookup$`Covariate Code`[i], ".tif")))
  }
}    

# create raster stack of continuous covariates
StackCovsC <- terra::rast(list(PopDen16, PopGro16, SocioEcon16_PC, DistRoad, DistCity, prop_value, AgProf, TSoilPC, elev, slope, prec, temp , EcolCond))
names(StackCovsC)[c(2, 10, 11, 15, 17, 18)] <- c("PopGro", "PropVal", "AgProf", "Elev", "Precip", "Temp")

# save raster stack
saveRDS(StackCovsC, file = "output/raster_stacks/cont_covs.rds")

# crop continuous covariate rasters by spatial units for each KMR
CropRastCovsC <- map(.x = SUs, .f = get_crop, Raster = StackCovsC)

rm(list = setdiff(ls(all.names = TRUE), c(ls(all.names = TRUE)[sapply(ls(all.names = TRUE), function(x) is.function(get(x)))], "SUs", "CropRastCovsC")))
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
gc()

# get continuous covariate values
ZStats_CovsC <- map2(.x = SUs, .y = CropRastCovsC, .f = get_zonal2, Stat = "mean")

# remove "mean" label from continuous covariates names
for (i in names(ZStats_CovsC)) {
  names(ZStats_CovsC[[i]]) <- names(ZStats_CovsC[[i]]) %>% str_remove("mean.")
}

# save data
saveRDS(ZStats_CovsC, file = "output/data/ZStats_CovsC.rds")

# Clearing stroage and memory space before processing Discrete variables.
rm(list = c(CropRastCovsC, ZStats_CovsC))
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
gc()

SUs <- readRDS("output/spatial_units/sus.rds")

# Load proposed covariates based on workshops from lookup xlsx
CovLookup <- readxl::read_xlsx("Input/covariates/covariate_description.xlsx", sheet = "AllLyr")

# load covariates
Cov_file_bn <- tools::file_path_sans_ext(basename(list.files("input/covariates", pattern = "\\.tif$", full.names = TRUE)))
for(i in 1:nrow(CovLookup)){
  if (CovLookup$`Covariate Code`[i] %in% Cov_file_bn){
    print(paste0("Load ... ", CovLookup$Type[i] ," ... ", CovLookup$`Covariate Code`[i], " = ",  CovLookup$`Covariate Description`[i]))
    assign(CovLookup$`Covariate Code`[i], terra::rast(paste0("input/covariates/", CovLookup$`Covariate Code`[i], ".tif")))
  }
}    

# create raster stack of discrete covariates
StackCovsD <- terra::rast(list(PolPref , LandTen, NSW_forten18_ForTen, NSW_forten18_ForType, NatVegReg , PlanZone, LandUse, drought , Fire   ))

# save raster stack
saveRDS(StackCovsD, file = "output/raster_stacks/disc_covs.rds")

# crop discrete covariate rasters by spatial units for each KMR
CropRastCovsD <- map(.x = SUs, .f = get_crop, Raster = StackCovsD)

# get discrete covariate values
ZStats_CovsD <- map2(.x = SUs, .y = CropRastCovsD, .f = get_zonal2, Stat = "mode") # note here could use Stat = "frac" to get the fraction of each discrete type in each property

# remove "mode" label from discrete covariates names and covert to factors
for (i in names(ZStats_CovsD)) {
  names(ZStats_CovsD[[i]]) <- names(ZStats_CovsD[[i]]) %>% str_remove("mode.")
  ZStats_CovsD[[i]] <- ZStats_CovsD[[i]] %>% mutate_all(~as.factor(.))
}

# save data
saveRDS(ZStats_CovsD, file = "output/data/ZStats_CovsD.rds")

# Clearing stroage and memory space before other processess if necessary
rm(ZStats_CovsD, CropRastCovsD, StackCovsD)
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE) # clear past and current temp file if necessary
gc()

# reclassify discrete covariates as appropriate

# FEI TO WRITE CODE TO DO THIS - E.G., SEE LINES 35-37

## This part onwards is for checking NAs and other summary stats ####
# Combine SUs and Woody
ZStats_Woody <- readRDS("output/data/ZStats_Woody.rds")
lapply(ZStats_Woody, summary)

ZStats_Woody_sf <- ZStats_Woody
for (i in names(ZStats_Woody_sf)) {
  ZStats_Woody_sf[[i]] <- ZStats_Woody_sf[[i]] %>% mutate(Shape = SUs[[i]]$Shape)
}

# combine SUs and Continuous covariates
ZStats_CovsC <- readRDS("output/data/ZStats_CovsC.rds")

lapply(ZStats_CovsC, summary)

ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC) 
summary(ZStats_CovsC_all)

ZStats_CovsC_sf <- ZStats_CovsC
for (i in names(ZStats_CovsC_sf)) {
  ZStats_CovsC_sf[[i]] <- ZStats_CovsC_sf[[i]] %>% mutate(Shape = SUs[[i]]$Shape)
}
ZStats_CovsC_sf <- do.call(rbind, ZStats_CovsC_sf) 

# Check for NA values in continuous covariates 
summary(ZStats_CovsC_sf)

# Replace NA values with 9999 for continuous covariates for plotting in ArcGIS
# ZStats_CovsC_sf <- ZStats_CovsC_sf %>% 
#   mutate(Elev = replace_na(Elev, 9999),
#          Income = replace_na(Income, 9999),
#          PopDen = replace_na(PopDen, 999999),
#          Precip = replace_na(Precip, 9999),
#          Slope = replace_na(Slope, 999),
#          SoilNit = replace_na(SoilNit, 999),
#          Temp = replace_na(Temp, 9999))

# st_write(ZStats_CovsC_sf, "input/ZStats_CovsC_sf.shp", append = FALSE)


# combine SUs and Discrete covariates
ZStats_CovsD <- readRDS("output/data/ZStats_CovsD.rds")
lapply(ZStats_CovsD, summary)

ZStats_CovsD_all <- do.call(rbind, ZStats_CovsD)
summary(ZStats_CovsD_all)

ZStats_CovsD_sf <- ZStats_CovsD
for (i in names(ZStats_CovsD_sf)) {
  ZStats_CovsD_sf[[i]] <- ZStats_CovsD_sf[[i]] %>% mutate(Shape = SUs[[i]]$Shape)
}
ZStats_CovsD_sf <- do.call(rbind, ZStats_CovsD_sf) 

# Check for NA values in Discrete covariates 
summary(ZStats_CovsD_sf)

# st_write(ZStats_CovsD_sf, "input/ZStats_CovsD_sf.shp", append = FALSE)
