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
#if (!require("diseasemapping")) install.packages("diseasemapping", repos="http://R-Forge.R-project.org")
#library(diseasemapping)

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

# proposed covariates based on workshops

#Agriculture:
#Land use
#Combined Drought Indicator
#Property size
#Property value
#Distance to nearest SUA
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

# create raster stack of continuous covariates
StackCovsC <- terra::rast(list(LVal = LVal, Elev = Elev, Nit = Nit))

# create raster stack of discrete covariates
StackCovsD <- terra::rast(list(LUse = LUse))

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

# load spatial units for SA1s
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

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)
#save data
saveRDS(ZStats_Woody, file = "output/data/ZStats_Woody.rds")

# crop continuous covariate rasters by spatial units for each KMR
CropRastCovsC <- map(.x = SUs, .f = get_crop, Raster = StackCovsC)

# crop discrete covariate rasters by spatial units for each KMR
CropRastCovsD <- map(.x = SUs, .f = get_crop, Raster = StackCovsD)

# get continuous covariate values
ZStats_CovsC <- map2(.x = SUs, .y = CropRastCovsC, .f = get_zonal, Stat = "mean")
#save data
saveRDS(ZStats_CovsC, file = "output/data/ZStats_CovsC.rds")

# get discrete covariate values
ZStats_CovsD <- map2(.x = SUs, .y = CropRastCovsD, .f = get_zonal, Stat = "mode") # note here could use Stat = "frac" to get the fraction of each discrete type in each property
#save data
saveRDS(ZStats_CovsD, file = "output/data/ZStats_CovsD.rds")

# TEST RUN OF AGRICULTURAL CLEARING MODEL FOR CENTRAL COAST

# for introducing spatially correlated random effects see: https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html

# remove "mean" label from continuous covariates names
names(ZStats_CovsC$CC) <- names(ZStats_CovsC$CC) %>% str_remove("mean.")

# remove "mode" label from discrete covariates names
names(ZStats_CovsD$CC) <- names(ZStats_CovsD$CC) %>% str_remove("mode.")

# check for multi-collinearity in continuous covariates
Corr_Cont <- cor(ZStats_CovsC$CC, use = "complete.obs")
write.csv(Corr_Cont, file="output/collinearity/cor_cont_CC.csv")

# get attribute table of spatial units and join covariates
Covs <- SUs$CC %>% st_drop_geometry() %>% as_tibble() %>% mutate(Area = Shape_Area / 10000) %>% select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(ZStats_CovsC$CC) %>% bind_cols(ZStats_CovsD$CC)  %>% mutate(SUID = 1:n())

# get response data
Response <- ZStats_Woody$CC %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody, N_now = sum.woody_now) %>%
              mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody, -sum.woody_now)

# set up data for INLA models

# response - how many cells cleared over time period
R <- Response %>% select(YAg) %>% as.matrix()
# number of trials - how many cells woody at start of time period
NT <- Response %>% select(N) %>% as.matrix()
# number of woody cells in 2023
W <- Response %>% select(N_now) %>% as.matrix()

# covariates - select relevant covariates here
C <- Covs %>% select(Area, LVal, Elev, Nit, SA1) %>% mutate(Area = as.numeric(scale(Area)), LVal = as.numeric(scale(LVal)), Elev = as.numeric(scale(Elev)), Nit = as.numeric(scale(Nit)))

# probability of clearing model

# format data for model fitting
RP <- R[which(NT > 0)] # only fit to data for properties with woody vegetation in 2011
RP <- as.vector(ifelse(RP > 0, 1, 0)) # recode to binary cleared/not cleared
NTP <- rep(1, length(RP)) # set number of trials to 1 for all properties in 2011
ResponseP <- as_tibble(cbind(RP, NTP))
names(ResponseP) <- c("P", "Ntrials")
CP <- C[which(NT > 0), ] # only fit to data for properties with woody vegetation
CP <- CP %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
DataP <- bind_cols(ResponseP, CP)

# get adjacency matrix for SA1s containing properties with forest cover
SA1IDs <- CP %>% select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
SA1Polys <- SA1s_All %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
Adj <- SA1Polys %>% get_adjacency("ModelP_CC_Adj_SA1s", "output/neighbours/")

# format data for model predictions
# this time we only include properties that have woody vegetation in 2023
RPPred <- R[which(W > 0)] # only make predictions for properties with woody vegetation in 2023
RPPred <- as.vector(ifelse(RPPred > 0, 1, 0)) # recode to binary cleared/not cleared
NTPPred <- rep(1, length(RPPred)) # set number of trials to 1 for all properties
ResponsePPred <- as_tibble(cbind(RPPred, NTPPred))
names(ResponsePPred) <- c("P", "Ntrials")
CPPred <- C %>% left_join(CP %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
CPPred <- CPPred[which(W > 0), ] # only make predictions for properties with woody vegetation in 2023
DataPPred <- bind_cols(ResponsePPred, CPPred)
DataPPred <- DataPPred %>% mutate(P = NA) # set whether cleared or not to NA so as to make predictions

# combine fitting and prediction data
DataP <- bind_rows(DataPPred, DataP)

# fit clearing versus no clearing model
formula <- P ~ Area + LVal + Elev + Nit + f(SA1ID, model = "bym", graph = Adj, scale.model = TRUE)
ResultP <- inla(formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
saveRDS(ResultP, file = "output/models/ModelP_CC_Ag.rds")

# inspect model
summary(ResultP)
ResultP$summary.fixed
head(ResultP$summary.random$SA1ID)
head(ResultP$summary.fitted.values)

# proportion cleared|clearing model

# format data for model fitting
RN <- R[which(R > 0)] # only fit to data from properties with some clearing
NTN <- NT[which(R > 0)] # only fit to data from properties with some clearing
ResponseN <- as_tibble(cbind(as.matrix(RN), as.matrix(NTN))) %>% mutate(Prop = RN / NTN) %>% mutate(Prop = ifelse(Prop < 1, Prop, 0.999)) # calculate the proportion cleared (when proportion = 1 set to 0.999)
names(ResponseN) <- c("N", "Ntrials", "Prop")
CN <- C[which(R > 0), ] # only fit to data from properties with some clearing
CN <- CN %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
DataN <- bind_cols(ResponseN, CN)

# get adjacency matrix for SA1s containing properties with forest cover
SA1IDs <- CN %>% select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
SA1Polys <- SA1s_All %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
Adj <- SA1Polys %>% get_adjacency("ModelN_CC_Adj_SA1s", "output/neighbours/")

# format data for model predictions
# this time we only include properties that have woody vegetation in 2023
RNPred <- R[which(W > 0)] # only make predictions for properties with woody vegetation in 2023
NTNPred <- NT[which(W > 0)] # only make predictions for properties with woody vegetation in 2023
ResponseNPred <- as_tibble(cbind(as.matrix(RNPred), as.matrix(NTNPred))) %>% mutate(Prop = NA) # set Prop to NA so we get predictions for these properties
names(ResponseNPred) <- c("N", "Ntrials", "Prop")
CNPred <- C %>% left_join(CN %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
CNPred <- CNPred[which(W > 0), ] # only make predictions for properties with woody vegetation in 2023
DataNPred <- bind_cols(ResponseNPred, CNPred)

# combine fitting and prediction data
DataN <- bind_rows(DataNPred, DataN)

# fit proportion cleared|clearing model
formula <- Prop ~ Area + LVal + Elev + Nit + f(SA1ID, model = "bym", graph = Adj, scale.model = TRUE) # remove random-effect for now
ResultN <- inla(formula, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
saveRDS(ResultN, file = "output/models/ModelN_CC_Ag.rds")

# look at model fit
summary(ResultN)
ResultN$summary.fixed
head(ResultN$summary.random$SA1ID)
head(ResultN$summary.fitted.values)

# spatialise predictions
Layer <- SUs$CC %>% bind_cols(R = as.vector(R), NT = as.vector(NT), W = as.vector(W), SUID = Covs$SUID) %>% dplyr::filter(W > 0) %>% mutate(ActualProp =  ifelse(NT > 0, R / NT, NA)) # only for properties with woody cover
PredictionsP <- as_tibble(ResultP$summary.fitted.values$mean[1:nrow(DataPPred)]) %>% bind_cols(SUID = as.vector(Covs$SUID[which(W > 0)]))
names(PredictionsP) <- c("PredP", "SUID")
PredictionsN <- as_tibble(ResultN$summary.fitted.values$mean[1:nrow(DataNPred)]) %>% bind_cols(SUID = as.vector(Covs$SUID[which(W > 0)]))
names(PredictionsN) <- c("PredN", "SUID")
PredictionsCombined <- PredictionsP %>% left_join(PredictionsN, by = join_by(SUID == SUID)) %>% mutate(PredAll = PredP * PredN)
Layer <- Layer %>% left_join(PredictionsCombined, by = join_by(SUID == SUID)) %>% dplyr::select(-Shape_Length, -Shape_Area)
st_write(Layer, "output/predictions/Pred_CC_Ag.shp", delete_layer = TRUE)





# OLD STUFF HERE

# set up data for combined likelihood INLA model - NOTE THIS DOESN'T WORK AT THE MOMENT SO NEED TO FIX - NOT SURE WHY

# format data for combined likelihood model
#DataAg <- get_zib_format(R, NT, C)

# fit combined ZIB model (hurdle model)
#formula <- cbind(P, N) ~ 0 + IntP + IntN + AreaP + AreaN + f(SA1IDP, model = "iid")
#ResultZIB <- inla(formula, data = DataAg, family = c("binomial", "zeroinflatedbinomial0"), Ntrials = Ntrials, verbose = TRUE, control.family=list(list(), list(hyper = list(prob = list(initial = -20,fixed = TRUE)))), control.inla = list(control.vb = list(enable = FALSE)))
