# THIS CODE READS IN THE REQUIRED DATA AND ORGANISES THE DATA FOR MODEL FITTING

# Clear memory space and hard drive space before preprocessing
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=FALSE)

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
do.call(rbind, SUs) %>% filter(Shape_Area < 0)

# save processes spatial units
saveRDS(SUs, file = "output/spatial_units/sus.rds")
qsave(SUs, file = "output/spatial_units/sus.qs", preset = "fast")

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

# Check for any NA values in SA1s
do.call(rbind, SA1s) %>% filter(if_all(everything(), is.na))

# save processed SA1s
saveRDS(SA1s, file = "output/spatial_units/sa1s.rds")
qsave(SA1s, file = "output/spatial_units/sa1s.qs", preset = "fast")

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

# read in Koala habitat suitability continuous covariate for extraction together with Woody cover stack
Khab <- rast("input/covariates/KoalaHabitatSuitability.tif")

# create raster stack of woody extent and loss
Stack <- terra::rast(list(aloss = Ag_Loss, iloss = In_Loss, floss = Fo_Loss, woody = Woody, Khab = Khab))

# save raster stack
saveRDS(Stack, file = "output/raster_stacks/woodyextloss.rds")

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal2, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)

# save data
saveRDS(ZStats_Woody, file = "output/data/ZStats_Woody.rds")
qsave(ZStats_Woody, file = "output/data/ZStats_Woody.qs", preset = "fast")

ZStats_Woody_all <- do.call(rbind, ZStats_Woody)
ZStats_Woody_all %>% filter(sum.woody < sum.Khab)


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

## Continuous covariate ---- 
SUs <- qread("output/spatial_units/sus.qs")

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

# Clear memory space and hard drive space before extracting continuous covariates
rm(list = setdiff(ls(all.names = TRUE), c(ls(all.names = TRUE)[sapply(ls(all.names = TRUE), function(x) is.function(get(x)))], "SUs", "CropRastCovsC")))
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=FALSE)
gc()

# get continuous covariate values
ZStats_CovsC <- map2(.x = SUs, .y = CropRastCovsC, .f = get_zonal2, Stat = "mean")

# remove "mean" label from continuous covariates names
for (i in names(ZStats_CovsC)) {
  names(ZStats_CovsC[[i]]) <- names(ZStats_CovsC[[i]]) %>% str_remove("mean.")
}

# save data
saveRDS(ZStats_CovsC, file = "output/data/ZStats_CovsC.rds")
qsave(ZStats_CovsC, file = "output/data/ZStats_CovsC.qs", preset = "fast")

# Clearing stroage and memory space before processing Discrete variables.
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
# tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
gc()

## Discrete covariate ---- 
# load functions
source("functions.R")
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

# Manually read in  raster discrete covariates (Use this method to read in layers if the for-loop does not work or to read in additional layers)
Remoteness <- rast("input/covariates/remote2016.tif")

# create raster stack of discrete covariates
StackCovsD <- terra::rast(list(PolPref , LandTen, NSW_forten18_ForTen, NSW_forten18_TenType, NSW_forten18_ForType, NatVegReg = NatVegReg , PlanZone, LandUse, drought , Fire  , Remoteness ))
names(StackCovsD)
# save raster stack
saveRDS(StackCovsD, file = "output/raster_stacks/disc_covs.rds")

# crop discrete covariate rasters by spatial units for each KMR
CropRastCovsD <- map(.x = SUs, .f = get_crop, Raster = StackCovsD)

# get discrete covariate values
ZStats_CovsD <- map2(.x = SUs, .y = CropRastCovsD, .f = get_zonal, Stat = "mode") # note here could use Stat = "frac" to get the fraction of each discrete type in each property

# remove "mode" label from discrete covariates names and covert to factors
for (i in names(ZStats_CovsD)) {
  names(ZStats_CovsD[[i]]) <- names(ZStats_CovsD[[i]]) %>% str_remove("mode.")
  ZStats_CovsD[[i]] <- ZStats_CovsD[[i]] %>% mutate_all(~as.factor(.))
}

# save data
saveRDS(ZStats_CovsD, file = "output/data/ZStats_CovsD.rds")
qsave(ZStats_CovsD, file = "output/data/ZStats_CovsD.qs", preset = "fast")

