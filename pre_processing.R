# THIS CODE READS IN THE REQUIRED DATA AND ORGANISES THE DATA FOR MODEL FITTING

# load libraries
library(sf)
library(sp)
library(raster)
library(exactextractr)
library(terra)
library(tidyverse)
library(spdep)
if (!require("INLA")) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

# load functions
source("functions.R")

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

# proposed covariates based on workshops

#Agriculture:
#Land use
#Combined Drought Indicator
#Property size
#Property value
#Distance to nearest SUA (significant urban area)
#Ecological condition
#Precipitation
#Temperature
#Slope
#Soil fertility

#Infrastructure:
#Land use
#Property size
#Property value
#Distance to nearest SUA
#Population density
#Ecological condition
#Slope

#Forestry:
#Land use
#Combined Drought Indicator
#Forest tenure
#Property size
#Property value
#Household income
#Ecological condition
#Fire
#Slope

# load covariates
# PROBABLY NEEDS TO CHANGE - FEI TO LOOK AT

# combined drought indicator
Drought <- terra::rast("input/covariates/drought.tif")
# elevation
Elev <- terra::rast("input/covariates/elev.tif")
# fire history (2011-2019)
Fire <- terra::rast("input/covariates/fire.tif")
# forest code
ForCode <- terra::rast("input/covariates/forest_code.tif")
# forest tenure
ForTen <- terra::rast("input/covariates/forest_tenure.tif")
# forest tenure type
ForTenType <- terra::rast("input/covariates/forest_tenure_type.tif")
# 2016 median household income weekly
Income <- terra::rast("input/covariates/income.tif")
# land use (NSW landuse 2013)
LandUse <- terra::rast("input/covariates/landuse.tif")
# population density
PopDen <- terra::rast("input/covariates/pop_den.tif")
# precipitation
Precip <- terra::rast("input/covariates/prec.tif")
# remoteness structure (ABS)
Remote <- terra::rast("input/covariates/remoteness.tif")
# slope
Slope <- terra::rast("input/covariates/slope.tif")
# soil fertility
SoilFert <- terra::rast("input/covariates/soil_fert.tif")
# soil nitrogen
SoilNit <- terra::rast("input/covariates/soil_nitrogen.tif")
# soil type
SoilType <- terra::rast("input/covariates/soil_type.tif")
# temperature
Temp <- terra::rast("input/covariates/temp.tif")

# reclassify discrete covariates as appropriate

# FEI TO WRITE CODE TO DO THIS - E.G., SEE LINES 35-37

# create raster stack of continuous covariates
StackCovsC <- terra::rast(list(Elev = Elev, Income = Income, PopDen = PopDen, Precip = Precip, Slope = Slope, SoilNit = SoilNit, Temp = Temp))

# save raster stack
saveRDS(StackCovsC, file = "output/raster_stacks/cont_covs.rds")

# create raster stack of discrete covariates
StackCovsD <- terra::rast(list(Drought = Drought, Fire = Fire, ForCode = ForCode, ForTen = ForTen, ForTenType = ForTenType, LandUse = LandUse, Remote = Remote, SoilFert = SoilFert, SoilType = SoilType))

# save raster stack
saveRDS(StackCovsD, file = "output/raster_stacks/disc_covs.rds")

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
  SUs[[i]] <- SUs[[i]] %>% mutate(Area = Shape_Area / 10000)
}

# save perocesses spatial units
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

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)

# save data
saveRDS(ZStats_Woody, file = "output/data/ZStats_Woody.rds")

# crop continuous covariate rasters by spatial units for each KMR
CropRastCovsC <- map(.x = SUs, .f = get_crop, Raster = StackCovsC)

# crop discrete covariate rasters by spatial units for each KMR
CropRastCovsD <- map(.x = SUs, .f = get_crop, Raster = StackCovsD)

# get continuous covariate values
ZStats_CovsC <- map2(.x = SUs, .y = CropRastCovsC, .f = get_zonal, Stat = "mean")

# remove "mean" label from continuous covariates names
for (i in names(ZStats_CovsC)) {
  names(ZStats_CovsC[[i]]) <- names(ZStats_CovsC[[i]]) %>% str_remove("mean.")
}

# save data
saveRDS(ZStats_CovsC, file = "output/data/ZStats_CovsC.rds")

# get discrete covariate values
ZStats_CovsD <- map2(.x = SUs, .y = CropRastCovsD, .f = get_zonal, Stat = "mode") # note here could use Stat = "frac" to get the fraction of each discrete type in each property

# remove "mode" label from discrete covariates names and covert to factors
for (i in names(ZStats_CovsD)) {
  names(ZStats_CovsD[[i]]) <- names(ZStats_CovsD[[i]]) %>% str_remove("mode.")
  ZStats_CovsD[[i]] <- ZStats_CovsD[[i]] %>% mutate_all(~as.factor(.))
}

# save data
saveRDS(ZStats_CovsD, file = "output/data/ZStats_CovsD.rds")
