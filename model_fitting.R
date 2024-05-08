# THIS CODE READS IN THE REQUIRED DATA AND FITS FOREST LOSS RISK MODELS FOR EACH
# KOALA MODELLING REGION (KMR) IN NEW SOUTH WALES

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

# load processed data
ZStats_Woody <- readRDS(file = "output/data/ZStats_Woody.rds")
ZStats_CovsC <- readRDS(file = "output/data/ZStats_CovsC.rds")
ZStats_CovsD <- readRDS(file = "output/data/ZStats_CovsD.rds")
SUs <- readRDS(file = "output/spatial_units/sus.rds")
SA1s <- readRDS(file = "output/spatial_units/sa1s.rds")

# add area to the continuous covariates abd then rescale to have mean of zero and standard deviation of one
# and check for multi-collinearity
for (i in names(ZStats_CovsC)) {
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]] %>% mutate(Area = SUs[[i]]$Area)
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]] %>% mutate_all(~as.numeric(scale(.)))
  Corr_Cont <- cor(ZStats_CovsC[[i]], use = "complete.obs")
  write.csv(Corr_Cont, file = paste0("output/collinearity/cor_cont_", i, ".csv"))
}

# combine continuous and discrete covariates into one tibble
ZStats_Covs <- list()
for (i in 1:length(names(ZStats_CovsC))) {
    ZStats_Covs[[i]] <- bind_cols(ZStats_CovsC[[i]], ZStats_CovsD[[i]])
}
names(ZStats_Covs) <- names(ZStats_CovsC)

# DO SOMETHING HERE TO REMOVE VARIABLES THAT ARE COLLINEAR - FEI TO WORK ON THIS
# SUGGEST REMOVE VARIABLES WITH CORRELATION > 0.7 OR < -0.7 (or 0.6?)
# REMOVE RELEVANT COLUMNS FROM ZStats_CovsC FOR EACH KMR
# PROBABLY CAN'T COMPLETELY AUTOMATE THIS STEP
# EXAMPLE SUCH AS THIS TO SELECT SUBSET OF COVARIATES
for (i in names(ZStats_CovsC)) {
    ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% select(Elev)
}

# NEED TO CHANGE COVARIATES FOR EACH CLEARING TYPE HERE - FEI TO WORK ON THIS


# MODEL FITTING STEPS - STILL WORKING ON THIS

# fit models for each KMR

Test <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s)

# do model selection



# make predictions














# OLD STUFF - PLEASE KEEP THIS FOR NOW

# set up data for INLA models

# response - how many cells cleared over time period
R <- Response %>% select(YAg) %>% as.matrix()
# number of trials - how many cells woody at start of time period
NT <- Response %>% select(N) %>% as.matrix()

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
# this time we only include properties that have woody vegetation in 2011
RPPred <- R[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
RPPred <- as.vector(ifelse(RPPred > 0, 1, 0)) # recode to binary cleared/not cleared
NTPPred <- rep(1, length(RPPred)) # set number of trials to 1 for all properties
ResponsePPred <- as_tibble(cbind(RPPred, NTPPred))
names(ResponsePPred) <- c("P", "Ntrials")
CPPred <- C %>% left_join(CP %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
CPPred <- CPPred[which(NT > 0), ] # only make predictions for properties with woody vegetation in 2011
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
# this time we only include properties that have woody vegetation in 2011
RNPred <- R[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
NTNPred <- NT[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
ResponseNPred <- as_tibble(cbind(as.matrix(RNPred), as.matrix(NTNPred))) %>% mutate(Prop = NA) # set Prop to NA so we get predictions for these properties
names(ResponseNPred) <- c("N", "Ntrials", "Prop")
CNPred <- C %>% left_join(CN %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
CNPred <- CNPred[which(NT > 0), ] # only make predictions for properties with woody vegetation in 2011
DataNPred <- bind_cols(ResponseNPred, CNPred)

# combine fitting and prediction data
DataN <- bind_rows(DataNPred, DataN)

# fit proportion cleared|clearing model
formula <- Prop ~ Area + LVal + Elev + Nit + f(SA1ID, model = "bym", graph = Adj, scale.model = TRUE)
ResultN <- inla(formula, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
saveRDS(ResultN, file = "output/models/ModelN_CC_Ag.rds")

# inspect model
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

# OLD OLD STUFF HERE

# set up data for combined likelihood INLA model - NOTE THIS DOESN'T WORK AT THE MOMENT SO NEED TO FIX - NOT SURE WHY

# format data for combined likelihood model
#DataAg <- get_zib_format(R, NT, C)

# fit combined ZIB model (hurdle model)
#formula <- cbind(P, N) ~ 0 + IntP + IntN + AreaP + AreaN + f(SA1IDP, model = "iid")
#ResultZIB <- inla(formula, data = DataAg, family = c("binomial", "zeroinflatedbinomial0"), Ntrials = Ntrials, verbose = TRUE, control.family=list(list(), list(hyper = list(prob = list(initial = -20,fixed = TRUE)))), control.inla = list(control.vb = list(enable = FALSE)))
