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
library(spdep)
#if (!require("spDataLarge")) install.packages('spDataLarge', repos='https://nowosad.github.io/drat/', type='source')
#library(spDataLarge)
if (!require("INLA")) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
#library(bigDM)
if (!require("diseasmapping")) install.packages("diseasemapping", repos="http://R-Forge.R-project.org")

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

# create adjacency matrices - COMMENTED OUT FOR NOW AS AGACENCY MATRIX TOO LARGE
#Adjacency <- map2(.x = SUs, .y = names(SUs), .f = get_adjacency, FileLocation = "output/neighbours/")

# TEST RUN OF AGRICULTURAL CLEARING MODEL FOR CENTRAL COAST

# get attribute table of spatial units
Covs <- SUs$CC %>% st_drop_geometry() %>% as_tibble() %>% mutate(Area = Shape_Area / 10000) %>% select(-KMR, - X, -Y, -Shape_Length, -Shape_Area)

# add area ID
Covs$areaID <- 1:nrow(Covs)

# get response data
Response <- ZStats_Woody$CC %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
              mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>%select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody)

# set up data for INLA model
R <- Response %>% select(YAg) %>% as.matrix()
NT <- Response %>% select(N) %>% as.matrix()
C <- Covs %>% select(Area) %>% mutate(Area = as.numeric(scale(Area)))
DataAg <- get_zib_format(R, NT, C)

# fit clearing versus no clearing model on its own
formula <- P ~ 0 + IntP + AreaP
ResultP <- inla(formula, data = DataAg, family = "binomial", Ntrials = NtrialsP, control.inla=list(control.vb = list(enable = FALSE)), verbose = TRUE)

# fit amount of clearing | clearing model on its own
formula <- N ~ 0 + IntN + AreaN
ResultN <- inla(formula, data = DataAg, family = "zeroinflatedbinomial0", Ntrials = NtrialsN, verbose = TRUE, control.family=(list(hyper = list(prob = list(initial = -20, fixed = TRUE)))), control.inla = list(control.vb = list(enable = FALSE)))

# fit combined ZIB model - NOTE THIS DOESN'T WORK AT THE MOMENT SO NEED TO FIX
formula <- cbind(P, N) ~ 0 + IntP + IntN + AreaP + AreaN
ResultZIB <- inla(formula, data = DataAg, family = c("binomial", "zeroinflatedbinomial0"), Ntrials = Ntrials, verbose = TRUE, control.family=list(list(), list(hyper = list(prob = list(initial = -20,fixed = TRUE)))), control.inla = list(control.vb = list(enable = FALSE)))
