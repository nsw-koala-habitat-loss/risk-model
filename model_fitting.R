# THIS CODE READS IN THE REQUIRED DATA AND FITS FOREST LOSS RISK MODELS FOR EACH
# KOALA MODELLING REGION (KMR) IN NEW SOUTH WALES

# load libraries
library(sf)
library(sp)
library(raster)
#if (!require("arcgisbinding")) install.packages("arcgisbinding", repos="https://r.esri.com", type="win.binary")
#library(arcgisbinding)
library(exactextractr)
library(terra)
library(tidyverse)

# load functions
source("functions.R")

# load woody extent and woody loss rasters
Ag_Loss <- raster(paste(getwd(), "/input/clearing_data/agr_loss_1119.tif",sep=""))
In_Loss <- raster(paste(getwd(), "/input/clearing_data/inf_loss_1119.tif",sep=""))
Fo_Loss <- raster(paste(getwd(), "/input/clearing_data/for_loss_1119.tif",sep=""))
Woody <- raster(paste(getwd(), "/input/clearing_data/woody_11.tif",sep=""))

# load woody extent and woody loss rasters
Ag_Loss <- terra::rast(paste(getwd(), "/input/clearing_data/agr_loss_1119.tif",sep=""))
In_Loss <- terra::rast(paste(getwd(), "/input/clearing_data/inf_loss_1119.tif",sep=""))
Fo_Loss <- terra::rast(paste(getwd(), "/input/clearing_data/for_loss_1119.tif",sep=""))
Woody <- terra::rast(paste(getwd(), "/input/clearing_data/woody_11.tif",sep=""))

# create raster stack of woody extent and loss
Stack <- terra::rast(list(aloss=Ag_Loss, iloss=In_Loss, floss=Fo_Loss, woody=Woody))

# load spatial units for each KMR
SUs <- list(CC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Coast"),
          CST = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Southern_Tablelands"),
          DRP = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Darling_Riverine_Plains"),
          FW = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Far_West"),
          NC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "North_Coast"),
          NT = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northern_Tablelands"),
          NS = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northwest_Slopes"),
          R = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Riverina"),
          SC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "South_Coast"))

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)
