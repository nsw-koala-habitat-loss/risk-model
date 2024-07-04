# THIS CODE READS IN THE REQUIRED DATA AND FITS FOREST LOSS RISK MODELS FOR EACH
# KOALA MODELLING REGION (KMR) IN NEW SOUTH WALES

# load libraries
library(sf)
library(sp)
library(terra)
library(tidyverse)
library(spdep)
if (!require("INLA")) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA) #INLA version 24.06.27
library(GGally)
library(ggpubr)
library(confintr)
library(readxl)
library(tictoc)
if(!require("qs")) install.packages("qs")
library(qs)
library(furrr)
library(ltm)

# load functions
source("functions.R")

# load processed data
ZStats_Woody <- qread("output/data/ZStats_Woody.qs")
ZStats_CovsC <- readRDS(file = "output/data/ZStats_CovsC.rds")
ZStats_CovsD <- readRDS(file = "output/data/ZStats_CovsD.rds")
SUs <- qread("output/spatial_units/sus.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")

# # Load proposed covariates based on workshops from lookup xlsx
# CovLookup <- readxl::read_xlsx("Input/covariates/covariate_description.xlsx", sheet = "AllLyr")


# add area to the continuous covariates abd then rescale to have mean of zero and standard deviation of one
for (i in names(ZStats_CovsC)) {
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]] %>% mutate(Area = SUs[[i]]$Area)
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]] %>% mutate(across(-matches("PC"), ~as.numeric(scale(.))))
}

# combine continuous and discrete covariates into one tibble
ZStats_Covs <- list()
for (i in 1:length(names(ZStats_CovsC))) {
  ZStats_Covs[[i]] <- bind_cols(ZStats_CovsC[[i]], ZStats_CovsD[[i]])
}
names(ZStats_Covs) <- names(ZStats_CovsC)

SUs_all <- do.call(rbind, SUs)
ZStats_Covs_all <- do.call(rbind, ZStats_Covs)
ZStats_Woody_all <- do.call(rbind, ZStats_Woody)
dim(SUs_all)
dim(ZStats_Covs_all)
dim(ZStats_Woody_all)

#'####################################
# CHECKING FOR MULTI-COLLINEARITY ####
#'####################################
ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC)
Corr_Cont <- cor(ZStats_CovsC_all, use = "complete.obs")
Corr_Cont_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
ggsave(Corr_Cont_all_plot, file = "output/collinearity/Corr_Cont_all_plot.png", width = 2000, height = 2000, units = "px")

ZStats_CovsD_all <- do.call(rbind, ZStats_CovsD)
str(ZStats_CovsD_all)

cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_all), ncol = ncol(ZStats_CovsD_all), 
                           dimnames = list(names(ZStats_CovsD_all), names(ZStats_CovsD_all)))

for (i in 1:ncol(ZStats_CovsD_all)) {
  for(j in 1:ncol(ZStats_CovsD_all)) {
    cramersv_val <- cramersv(as.data.frame(ZStats_CovsD_all[c(i,j)]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

Corr_Categ_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
ggsave(Corr_Categ_all_plot, file = "output/collinearity/Corr_Categ_all_plot.png", width = 2000, height = 2000, units = "px")

ZStats_Covs_all <- do.call(rbind, ZStats_Covs)
biserial.cor(ZStats_Covs_all$Temp, ZStats_Covs_all$Drought, use = "complete.obs", level = 1)
library(ltm)
biserial.cor(y, x, use = c("all.obs"), level = 2)
0.5847499


#'###############################################################################

# Remove covariates that are collinear ----
## SUGGEST REMOVE VARIABLES WITH CORRELATION > 0.7 OR < -0.7 (or 0.6?) -DONE
for (i in names(ZStats_Covs)) {
    ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% dplyr::select(-c(Elev, ForType))
}

# Select SUs and covariates for each ClearType ----
### This is done based on info from DAGs constructed throught expert elicitation
### DO ONCE AND SAVE RESULTS

# Filter All SUs to remove "Water" in Land Use
# SUs <- qread("output/spatial_units/sus.qs")
for(i in names(SUs)){
  SUs[[i]] <- SUs[[i]] %>% dplyr::select(-Area) %>% cbind(ZStats_Covs[[i]])  %>% filter(LandUse !="5")
  SUs[[i]] <- SUs[[i]] %>% dplyr::select(-names(ZStats_Covs[[i]]))
  ZStats_Woody[[i]] <- cbind(ZStats_Woody[[i]], ZStats_Covs[[i]]) %>% filter(LandUse !="5")
  ZStats_Woody[[i]] <- ZStats_Woody[[i]] %>% dplyr::select(-names(ZStats_Covs[[i]]))
  ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% filter(LandUse !="5")
}

## Agricultural Clearing ----
ZStats_Covs_Ag <- ZStats_Covs
for (i in names(ZStats_Covs_Ag)) {
    ZStats_Covs_Ag[[i]] <- ZStats_Covs_Ag[[i]] %>%
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, LandTen, NatVegReg, LandUse, Fire)
}

SUs_Ag <- SUs
ZStats_Woody_Ag <- ZStats_Woody
for(i in names(SUs_Ag)){
  SUs_Ag[[i]] <- cbind(SUs_Ag[[i]],ZStats_Covs_Ag[[i]])  %>% filter(NatVegReg != "0")
  SUs_Ag[[i]] <- SUs_Ag[[i]] %>% dplyr::select(-names(ZStats_Covs_Ag[[i]]))
  ZStats_Woody_Ag[[i]] <- cbind(ZStats_Woody_Ag[[i]],ZStats_Covs_Ag[[i]]) %>% filter(NatVegReg != "0")
  ZStats_Woody_Ag[[i]] <- ZStats_Woody_Ag[[i]] %>% dplyr::select(-names(ZStats_Covs_Ag[[i]]))
  ZStats_Covs_Ag[[i]] <- ZStats_Covs_Ag[[i]] %>% filter(NatVegReg != "0")
}

qsave(SUs_Ag, file = "output/spatial_units/sus_ag.qs", preset = "fast")
qsave(ZStats_Woody_Ag, file = "output/data/ZStats_Woody_Ag.qs", preset = "fast")
qsave(ZStats_Covs_Ag, file = "output/data/ZStats_Covs_Ag.qs", preset = "fast")


## Industrial Clearing ----
ZStats_Covs_In <- ZStats_Covs
for (i in names(ZStats_Covs_In)) {
  ZStats_Covs_In[[i]] <- ZStats_Covs_In[[i]] %>%
    dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, PolPref, LandTen,  PlanZone, LandUse, Fire)
}
qsave(ZStats_Covs_In, file = "output/data/ZStats_Covs_In.qs", preset = "fast")
ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")

## Forestry Clearing ----
ZStats_Covs_Fo <- ZStats_Covs
for (i in names(ZStats_Covs_Fo)) {
  ZStats_Covs_Fo[[i]] <- ZStats_Covs_Fo[[i]] %>%
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, PolPref, LandTen, ForTen, NatVegReg, LandUse, Fire)
}
qsave(ZStats_Covs_Fo, file = "output/data/ZStats_Covs_Fo.qs", preset = "fast")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")




# MODEL FITTING STEPS ----
## Load data
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
# ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")
# ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")


# Fit models for each KMR
# Test <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s)
TEST2 <- fit_model2(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)

# Test DIC variation with different CONTROL INLA Strategy  (Using deconstructed fit_model function) ----
## Store Control.inla default 
Control.INLA_Default <- control.inla()
## Model input
KMR = "CC"; ClearType = 1 ; SpatUnits = SUs_Ag ; RespData = ZStats_Woody_Ag ; CovsCD = ZStats_Covs_Ag ; SA1sPoly = SA1s; Explanatory = "All"; Verbose = FALSE; N_retry=3; Initial_Tlimit = 1000; OutputDir = NULL 

# get attribute table of spatial units and join covariates
Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())

# get response data and ensure clearing is < woody vegetation
Response <- RespData[[KMR]] %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
  mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody)

# set up data for INLA models

# response - how many cells cleared over time period
if (ClearType == 1) {
  R <- Response %>% dplyr::select(YAg) %>% as.matrix()
} else if (ClearType == 2) {
  R <- Response %>% dplyr::select(YIn) %>% as.matrix()
} else if (ClearType == 3) {
  R <- Response %>% dplyr::select(YFo) %>% as.matrix()
}

# number of trials - how many cells woody at start of time period
NT <- Response %>% dplyr::select(N) %>% as.matrix()


## Format INLA Parameters ----
### Probability of clearing model ----
RP <- R[which(NT > 0)] # only fit to data for properties with woody vegetation in 2011
RP <- as.vector(ifelse(RP > 0, 1, 0)) # recode to binary cleared/not cleared
NTP <- rep(1, length(RP)) # set number of trials to 1 for all properties in 2011
ResponseP <- as_tibble(cbind(RP, NTP))
names(ResponseP) <- c("P", "Ntrials")
CP <- Covs[which(NT > 0), ] # only fit to data for properties with woody vegetation
CP <- CP %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
DataP <- bind_cols(ResponseP, CP)
ExplV <- if(Explanatory == "All") {paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + ")} else {paste(Explanatory, , collapse=" + ")}

# get adjacency matrix for SA1s containing properties with forest cover
SA1IDs <- CP %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
SA1sPolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
AdjP <- SA1sPolyKMR %>% get_adjacency(paste0("modelP_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")

### Proportion cleared|clearing model ----
RN <- R[which(R > 0)] # only fit to data from properties with some clearing
NTN <- NT[which(R > 0)] # only fit to data from properties with some clearing
ResponseN <- as_tibble(cbind(as.matrix(RN), as.matrix(NTN))) %>% mutate(Prop = RN / NTN) %>% mutate(Prop = ifelse(Prop < 1, Prop, 0.999)) # calculate the proportion cleared (when proportion = 1 set to 0.999)
names(ResponseN) <- c("N", "Ntrials", "Prop")
CN <- Covs[which(R > 0), ] # only fit to data from properties with some clearing
CN <- CN %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
DataN <- bind_cols(ResponseN, CN)
# get adjacency matrix for SA1s containing properties with some clearing
SA1IDs <- CN %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
SA1PolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
AdjN <- SA1sPolyKMR %>% get_adjacency(paste0("modelN_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")


ExplV <- if(Explanatory == "All") {paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + ")} else {paste(Explanatory, , collapse=" + ")}

# fit clearing versus no clearing model
formula <- as.formula(paste0(paste("P", ExplV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)

Control.INLA_ORI <- control.inla(strategy = "auto", int.strategy = "auto", control.vb = INLA::control.vb(enable = FALSE))
tic("INLA_ORI")
for(i in 1:20){
  ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = Control.INLA_ORI, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultP$dic$dic))
}
toc(log = TRUE)

# fit proportion cleared|clearing model
formula <- as.formula(paste0(paste("Prop", ExplV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))

Control.INLA_ORI <- control.inla(strategy = "auto", int.strategy = "auto", control.vb = INLA::control.vb(enable = FALSE))
tic("INLA_ORI")
for(i in 1:20){
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_ORI, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultN$dic$dic))
}
toc(log = TRUE)
summary(ResultN)


Control.INLA_FASTER <- control.inla(control.vb = INLA::control.vb(enable = FALSE), strategy = "simplified.laplace" , int.strategy = "eb")
tic("INLA_FASTER")
for(i in 1:20){
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_FASTER, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultN$dic$dic))
}
toc(log = TRUE)

tic("INLA_FAST")
Control.INLA_FAST <- control.inla(control.vb = INLA::control.vb(enable = FALSE), strategy = "laplace" , int.strategy = "eb")
tic("INLA_FAST")
for(i in 1:20){
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_FAST, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultN$dic$dic))
}
toc(log = TRUE)

Control.INLA_SLOW <- control.inla(control.vb = INLA::control.vb(enable = FALSE),  strategy = "simplified.laplace" , int.strategy = "ccd")
tic("INLA_SLOW")
for(i in 1:20){
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_SLOW, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultN$dic$dic))
}
toc(log = TRUE)

Control.INLA_SLOWER <- control.inla(control.vb = INLA::control.vb(enable = FALSE), strategy = "laplace" , int.strategy = "ccd", restart = 3L)
tic("INLA_SLOWER")
for(i in 1:20){
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_SLOWER, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  print(paste0("DIC= " , ResultN$dic$dic))
}
toc(log = TRUE)

ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_SLOW, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)

## Test DIC variation with different initial value input  - This part does not work!!! yet? 
## Setting custom prior based on Full model - This part does not work!!! yet?
PriorN <- paste("table:", 
                paste(seq(0, 1, len = length(ResultN$summary.fitted.values$mean)), collapse = " "), 
                paste(as.vector(ResultN$summary.fitted.values$mean), collapse = " "), 
                collapse = " ")
ResultN2 <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = Control.INLA_SLOW, control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose, control.family = list(hyper = list(prec = list(prior = PriorN))))


ptm <- proc.time()
Test_2 <- fit_model2(KMR="CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, OutputDir = "output/models/")
proc.time() - ptm
Test_2 <- qread("output/models/Model_CC_Ag_.qs")
summary(Test_2$PModel)

formula_1 <- as.formula(paste0(paste("P", paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))


summary(Test_2$PModel)
Test_2$PModel$cpu.used

ResultP$cpu[1]

# do model selection  ## Fei to try writing this part ##

## if 
## 1 DIC of model with full set of covariates
## 2 DIC of model with reduced set of covariates
## if 2 < 1 then keep reduced set of covariates
## if 2 > 1 then keep full set of covariates
## Continue until no further reduction in DIC

ptm <- proc.time()
Test_3 <- Select_model(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, Selection =  "F", Verbose = FALSE, inla_retry = TRUE)
proc.time() - ptm
# BEst Forward selection model
# $`PopDen + LandUse + EcolCond + NatVegReg + Area + LandTen + DistCity + slope + Soil_PC1 + ScEc_PC4 + PropVal + Soil_PC3 + ScEc_PC3 + Temp`
# [1] 13773.93

# P ~  ScEc_PC2  + ScEc_PC5 + DistRoad + Soil_PC1 + Soil_PC2 + slope +      Precip  + EcolCond + Area + LandTen + NatVegReg  +      Fire + f(SA1ID, model = "bym", graph = AdjP, scale.model = TRUE)

F_test2 <- INLA_with_TimeLimit(TimeLimit = 1000, ForN_H1, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
F_test1 <- INLA_with_Retry(N_retry = 3, Initial_Tlimit = 1000, ForN_H1, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)

ptm <- proc.time()
Test_5 <- Select_model(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, Direction = "Backward", Verbose = FALSE)
proc.time() - ptm

ptm <- proc.time()
Test_6 <- Select_model(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, Direction = "Complete Forward", Verbose = FALSE, OutputDir = "output/models/")
proc.time() - ptm

ptm <- proc.time()
Test_7 <- Select_model(KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "Forward Complete", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
proc.time()-ptm

Test_7$Best_DIC_ls

ptm <- proc.time()
Test_8 <- Select_model(KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "Backward Complete", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
proc.time()-ptm

Test_8$Best_DIC_ls

tic()
for(i in 1:20){
  Mod <- fit_model2 (KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
  print(paste0("DIC= " , Mod$PModel$dic$dic+Mod$NModel$dic$dic))
}
toc()

tic()
for(i in 1:20){
  Mod <- fit_model2 (KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = c("PopDen", "DistCity", "PropVal", "Precip", "EcolCond",  "Area" , "LandUse"), Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
  print(paste0("DIC= " , Mod$PModel$dic$dic+Mod$NModel$dic$dic))
}
toc()

summary(ResultN)


Test_3A <- fit_model2(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "PopDen + ScEc_PC1 + ScEc_PC2 + ScEc_PC3 + ScEc_PC4 + ScEc_PC5 +      DistRoad + DistCity + PropVal + Soil_PC1 + Soil_PC2 + Soil_PC3 +      slope + Precip + Temp + EcolCond + Area + LandTen + NatVegReg +      LandUse + Fire", Verbose = FALSE)
P ~ PopDen + ScEc_PC1 + ScEc_PC2 + ScEc_PC3 + ScEc_PC4 + ScEc_PC5 + DistRoad + DistCity + PropVal + AgProf + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Precip + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire + f(SA1ID, model = "bym", graph = AdjP,scale.model = TRUE)

# make predictions

# JRR EXAMPLE FOR MAKING PREDICTIONS

# fit model
Fit <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s)

# get predictions from model
Preds <- predict_model(Fit)

# visualise predictions
plot(Preds$Layer)

# save predictions
st_write(Preds$Layer, "output/predictions/test.shp", delete_layer = TRUE)

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
