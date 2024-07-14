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
qsave(ZStats_Covs, file = "output/data/ZStats_Covs.qs", preset = "fast")

# Export NA SUs shapefile for checking
SUs_Covs_All <- cbind(do.call(rbind, SUs), do.call(rbind, ZStats_Covs)) 
SUs_Covs_All_Ag <- SUs_Covs_All %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, LandTen, NatVegReg, LandUse, Fire) %>% 
  filter_all(any_vars(is.na(.)))
st_write(SUs_Covs_All_Ag, "output/diagnosis/SUs_Covs_All_Ag_NA.shp", delete = TRUE)

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
ZStats_Woody <- qread("output/data/ZStats_Woody.qs")
ZStats_Covs <- qread("output/data/ZStats_Covs.qs")
for (i in names(ZStats_Covs)) {
    ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% dplyr::select(-c(Elev, ForType))
}

# Filter SUs and covariates for each ClearType ----

### This is done based on info from DAGs constructed throught expert elicitation
### Filter All SUs to remove "Water" in Land Use
### Filter All SUs to remove NAs in SA1

## Agricultural Clearing ----
## DAG Covariate
# PopDen16
# PolPref
# SocioEcon16_PC
# DistRoad
# DistCity
# LandTen
# LandUse
# prop_value
# AgProf
# TSoilPC
# elev (rm)
# slope
# prec
# temp
# drought
# fire
# EcolCond

SUs <- qread("output/spatial_units/sus.qs")
SUs_Ag <- SUs_ZStats_Ag <- SUs 
ZStats_Woody_Ag <- ZStats_Woody
ZStats_Covs_Ag <- ZStats_Covs
for(i in names(SUs_ZStats_Ag)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_Ag[[i]] <- SUs_Ag[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Agricultural Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  ZStats_Covs_Ag[[i]] <- ZStats_Covs_Ag[[i]] %>% 
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, LandTen, NatVegReg, LandUse, Drought, Fire)
  
  # Merge SUs, Covs and Woody %>% remove water land use and NAs in SA1
  SUs_ZStats_Ag[[i]] <- bind_cols(SUs_Ag[[i]], ZStats_Covs_Ag[[i]], ZStats_Woody_Ag[[i]]) %>% 
    filter(LandUse != "5", NatVegReg != "0") %>% 
    tidyr::drop_na(SA1) %>% 
    tidyr::drop_na(.)
  
  # Split SUs, Covs and Woody into individual dataframes
  SUs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% dplyr::select(names(SUs_Ag[[i]]))
  ZStats_Covs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Ag[[i]]))
  ZStats_Woody_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Ag[[i]]))
}

### Check for collinearity in Agricultural Clearing----
ZStats_CovsC_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)
Corr_Cont_Ag_all <- cor(ZStats_CovsC_Ag_all, use = "complete.obs")
Corr_Cont_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Ag_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_Cont_Ag_all_plot
# ggsave(Corr_Cont_Ag_all_plot, file = "output/collinearity/Corr_Cont_Ag_all_plot.png", width = 2000, height = 2000, units = "px")

ZStats_CovsD_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(LandTen, NatVegReg, LandUse, Drought, Fire)

cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_Ag_all), ncol = ncol(ZStats_CovsD_Ag_all), 
                           dimnames = list(names(ZStats_CovsD_Ag_all), names(ZStats_CovsD_Ag_all)))

for (i in 1:ncol(ZStats_CovsD_Ag_all)) {
  for(j in 1:ncol(ZStats_CovsD_Ag_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_Ag_all[[i]], ZStats_CovsD_Ag_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

Corr_CovD_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_CovD_Ag_all_plot
# ggsave(Corr_CovD_Ag_all_plot, file = "output/collinearity/Corr_CovD_Ag_all_plot.png", width = 2000, height = 2000, units = "px")


qsave(SUs_Ag, file = "output/spatial_units/sus_ag.qs", preset = "fast")
qsave(ZStats_Woody_Ag, file = "output/data/ZStats_Woody_Ag.qs", preset = "fast")
qsave(ZStats_Covs_Ag, file = "output/data/ZStats_Covs_Ag.qs", preset = "fast")

# SUs_Woody_Covs_all <- cbind(do.call(rbind, SUs_Ag), do.call(rbind, ZStats_Woody_Ag), do.call(rbind, ZStats_Covs_Ag))
# SUs_Woody_Covs_all_NA <- SUs_Woody_Covs_all %>% filter_all(any_vars(is.na(.)))
# st_write(SUs_Woody_Covs_all_NA, "output/diagnosis/SUs_Woody_Covs_all_NA.shp", delete = TRUE)
# SUs_Ag_all_v <- vect(SUs_Ag_all)
# ggplot()+tidyterra::geom_spatvector(data = SUs_Ag_all_v, aes(fill = KMR)) + theme_minimal()

NA_Val <- t(map_dfr(ZStats_Covs_Ag, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(ZStats_Covs_Ag, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 3))))
colnames(NA_Val) <- names(ZStats_Covs_Ag)
colnames(NA_Pct) <- names(ZStats_Covs_Ag)
NA_Val
NA_Pct


NA_Val <- t(map_dfr(SUs_Ag, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(SUs_Ag, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 5))))
colnames(NA_Val) <- names(SUs_Ag)
colnames(NA_Pct) <- names(SUs_Ag)
NA_Val
NA_Pct

## Forestry Clearing ----
## DAG_Covariate
# PopDen16
# PolPref
# SocioEcon16_PC
# DistRoad
# DistCity
# LandTen
# NSW_forten18_ForTen
# NSW_forten18_ForType
# NatVegReg
# LandUse
# prop_value
# AgProf
# TSoilPC
# elev
# slope
# prec
# temp
# drought
# fire
# EcolCond
# Area

SUs <- qread("output/spatial_units/sus.qs")
SUs_Fo <- SUs_ZStats_Fo <- SUs 
ZStats_Woody_Fo <- ZStats_Woody
ZStats_Covs_Fo <- ZStats_Covs

for(i in names(SUs_ZStats_Fo)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_Fo[[i]] <- SUs_Fo[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Forestry Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  # Drop LandTen as the infor captured in TenType and Drop ForTen with limited SUs in each category 
  ZStats_Covs_Fo[[i]] <- ZStats_Covs_Fo[[i]] %>% 
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, NatVegReg, LandUse, Drought, Fire, TenType)
  
  # Merge SUs, Covs and Woody 
  SUs_ZStats_Fo[[i]] <- bind_cols(SUs_Fo[[i]], ZStats_Covs_Fo[[i]], ZStats_Woody_Fo[[i]]) %>% 
    
    # Remove water land use type and NAs in SA1
    filter(LandUse != "5") %>% 
    tidyr::drop_na(SA1) %>% 
    
    # Remove SUs in LLS excluded area (non-rural) 
    # Remove Multiple-use public forest(2) and Nature conservation reserve (3) Forest Tenure Type
    filter(NatVegReg != "0", TenType == "1") %>% 
    
    # Remove NAs
    tidyr::drop_na(.)
  
  # Split SUs, Covs and Woody into individual dataframes
  SUs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% dplyr::select(names(SUs_Fo[[i]]))
  ZStats_Covs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Fo[[i]])) %>% dplyr::select(-TenType)
  ZStats_Woody_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Fo[[i]]))
}

### Check for collinearity in Forestry Clearing----
ZStats_CovsC_Fo_all <- do.call(rbind, ZStats_Covs_Fo) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)
Corr_Cont_Fo_all <- cor(ZStats_CovsC_Fo_all, use = "complete.obs")
Corr_Cont_Fo_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Fo_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_Cont_Fo_all_plot
ggsave(Corr_Cont_Fo_all_plot, file = "output/collinearity/Corr_Cont_Fo_all_plot.png", width = 2000, height = 2000, units = "px")

ZStats_CovsD_Fo_all <- do.call(rbind, ZStats_Covs_Fo) %>% 
  dplyr::select(NatVegReg, LandUse, Fire, Drought )

cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_Fo_all), ncol = ncol(ZStats_CovsD_Fo_all), 
                           dimnames = list(names(ZStats_CovsD_Fo_all), names(ZStats_CovsD_Fo_all)))

for (i in 1:ncol(ZStats_CovsD_Fo_all)) {
  for(j in 1:ncol(ZStats_CovsD_Fo_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_Fo_all[[i]], ZStats_CovsD_Fo_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

Corr_CovD_Fo_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_CovD_Fo_all_plot
ggsave(Corr_CovD_Fo_all_plot, file = "output/collinearity/Corr_CovD_Fo_all_plot.png", width = 2000, height = 2000, units = "px")

qsave(SUs_Fo, file = "output/spatial_units/sus_fo.qs", preset = "fast")
qsave(ZStats_Woody_Fo, file = "output/data/ZStats_Woody_Fo.qs", preset = "fast")
qsave(ZStats_Covs_Fo, file = "output/data/ZStats_Covs_Fo.qs", preset = "fast")

SUs_Fo <- qread("output/spatial_units/sus_fo.qs")
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
names(ZStats_Covs_Fo)

map(ZStats_Covs_Fo, ~summary(.))

NA_Val <- t(map_dfr(ZStats_Covs_Fo, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(ZStats_Covs_Fo, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 3))))
colnames(NA_Val) <- names(ZStats_Covs_Fo)
colnames(NA_Pct) <- names(ZStats_Covs_Fo)
NA_Val
NA_Pct


NA_Val <- t(map_dfr(SUs_Fo, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(SUs_Fo, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 5))))
colnames(NA_Val) <- names(SUs_Fo)
colnames(NA_Pct) <- names(SUs_Fo)
NA_Val
NA_Pct

## Infrastructure Clearing ----
## DAG_Covariate
# PopDen16
# PopGro16
# PolPref
# SocioEcon16_PC
# DistRoad
# DistCity
# LandTen
# PlanZone
# LandUse
# prop_value
# TSoilPC
# elev
# slope
# prec
# temp
# drought
# fire
# EcolCond
# Area
# NatVegReg

SUs <- qread("output/spatial_units/sus.qs")
SUs_In <- SUs_ZStats_In <- SUs
ZStats_Woody_In <- ZStats_Woody
ZStats_Covs_In <- ZStats_Covs
colnames(ZStats_Covs_In$CST)

for(i in names(SUs_ZStats_In)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_In[[i]] <- SUs_In[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Forestry Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  # Drop NatVegReg due as it's correlated to LandUse and PlanZone
  # Drop fire, and Soil due to uncertain direct effect on Infrastructure Clearing
  ZStats_Covs_In[[i]] <- ZStats_Covs_In[[i]] %>% 
    dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, slope, Precip, Temp, EcolCond, Area, LandTen, PlanZone, LandUse, Drought, Remote)
  
  # Merge SUs, Covs and Woody 
  SUs_ZStats_In[[i]] <- bind_cols(SUs_In[[i]], ZStats_Covs_In[[i]], ZStats_Woody_In[[i]]) %>% 
    
    # Remove "5" (water) and "1" (Conservation and natural environments) land use type and NAs in SA1
    filter(LandUse != "5") %>% 
    tidyr::drop_na(SA1) %>% 
    
    # Remove SUs in "0" (Major Cities), "1" (Inner Regional), "2" (Outer Regional)
    filter(Remote %in% c("0", "1", "2")) %>% 
    
    # Remove NAs
    tidyr::drop_na(.)
    
  
  # Split SUs, Covs and Woody into individual dataframes
  SUs_In[[i]] <- SUs_ZStats_In[[i]] %>% dplyr::select(names(SUs_In[[i]]))
  ZStats_Covs_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_In[[i]]))%>% dplyr::select(-Remote)
  ZStats_Woody_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_In[[i]])) 
}

### Check for collinearity in Forestry Clearing----
ZStats_CovsC_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, slope, Precip, Temp, EcolCond, Area)
Corr_Cont_In_all <- cor(ZStats_CovsC_In_all, use = "complete.obs")
Corr_Cont_In_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_In_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_Cont_In_all_plot
ggsave(Corr_Cont_In_all_plot, file = "output/collinearity/Corr_Cont_In_all_plot.png", width = 2000, height = 2000, units = "px")

ZStats_CovsD_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(LandTen, PlanZone, LandUse, Drought)

cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_In_all), ncol = ncol(ZStats_CovsD_In_all), 
                           dimnames = list(names(ZStats_CovsD_In_all), names(ZStats_CovsD_In_all)))

for (i in 1:ncol(ZStats_CovsD_In_all)) {
  for(j in 1:ncol(ZStats_CovsD_In_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_In_all[[i]], ZStats_CovsD_In_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

Corr_CovD_In_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_CovD_In_all_plot
ggsave(Corr_CovD_In_all_plot, file = "output/collinearity/Corr_CovD_In_all_plot.png", width = 2000, height = 2000, units = "px")

qsave(SUs_In, file = "output/spatial_units/sus_In.qs", preset = "fast")
qsave(ZStats_Woody_In, file = "output/data/ZStats_Woody_In.qs", preset = "fast")
qsave(ZStats_Covs_In, file = "output/data/ZStats_Covs_In.qs", preset = "fast")

map(ZStats_Covs_In, ~summary(.))
map(ZStats_Covs_In, ~nrow(.))
# 
# PlanZone_sum <- matrix(NA, nrow = length(as.factor(as.integer(0:6))), ncol = length(names(ZStats_Covs_In)),
#                        dimnames = list(as.factor(as.integer(0:6)), names(ZStats_Covs_In)))
# 
# PlanZone_pct <- PlanZone_sum
# sort(unique(ZStats_Covs_In$CC$PlanZone))
# for(i in 1:ncol(PlanZone_sum)){
#   for(j in sort(unique(ZStats_Covs_In$CC$PlanZone))){
#     PlanZone_sum[j,i] <- sum(ZStats_Covs_In[[i]]$PlanZone == j)
#     PlanZone_pct[j,i] <- round(sum(ZStats_Covs_In[[i]]$PlanZone == j)/sum(!is.na(ZStats_Covs_In[[i]]$PlanZone)),3)
#   }
# }
# 
# PlanZone_sum
# PlanZone_pct
# SUs_In_all <- do.call(rbind, SUs_In)
# Zstats_Covs_In_all_Ru <- do.call(rbind, ZStats_Covs_In) %>% 
#   filter(PlanZone == "5") %>% 
#   dplyr::select(PlanZone)
# 
# ZStats_Covs_In_Pz <- ZStats_Covs_In
# for(i in names(ZStats_Covs_In_Pz)){
#   ZStats_Covs_In_Pz[[i]] <- ZStats_Covs_In_Pz[[i]] %>% filter(PlanZone == "5") %>% dplyr::select(PlanZone)
# }
# map(ZStats_Covs_In_Pz, ~summary(.))
# 
# SUs_ZStats_PlanZone <- bind_cols(do.call(rbind, SUs_In), do.call(rbind, ZStats_Covs_In)) %>% 
#   dplyr::select(names(SUs_In$CC), LandTen, PlanZone)
# 
# st_write(SUs_ZStats_PlanZone, "output/diagnosis/SUs_ZStats_PlanZone.shp", delete = TRUE)

NA_Val <- t(map_dfr(ZStats_Covs_In, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(ZStats_Covs_In, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 3))))
colnames(NA_Val) <- names(ZStats_Covs_In)
colnames(NA_Pct) <- names(ZStats_Covs_In)
NA_Val
NA_Pct


NA_Val <- t(map_dfr(SUs_In, ~map(., ~sum(is.na(.)))))
NA_Pct <- t(map_dfr(SUs_In, 
                    ~map(., 
                         ~round((sum(is.na(.)))/length(.), 5))))
colnames(NA_Val) <- names(SUs_In)
colnames(NA_Pct) <- names(SUs_In)
NA_Val
NA_Pct



# MODEL FITTING STEPS ----

## Load data
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
# Test <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s)
TEST2 <- fit_model2(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)


# Test fit models for each KMR and ClearType ----
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Ag)

## Check DIC variations ----
# tic()
# for(i in 1:5){
#   Mod <- fit_model2 (KMR = "NC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "PopDen + ScEc_PC1 + AgProf + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
#   print(paste0("DIC= " , Mod$PModel$dic$dic+Mod$NModel$dic$dic))
# }
# toc()
inla.setOption(num.threads="2:1")
### Agricultural clearing in each KMR ----
ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:3){
    Model_Ag <- fit_model2(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_Ag$PModel$dic$dic + Model_Ag$NModel$dic$dic, 2)))
  }
  toc(log = TRUE)
  cat("\n\n")
}
proc.time() - ptm

### Forestry clearing in each KMR ----
SUs_Fo <- qread("output/spatial_units/sus_fo.qs")
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Fo)
summary(ZStats_Woody_Fo$FW$sum.woody)
ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:3){
    Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_Fo$PModel$dic$dic + Model_Fo$NModel$dic$dic, 3)))
  }
  toc(log = TRUE)
  cat("\n\n")
}
proc.time() - ptm

### Infrastructure clearing in each KMR ----
SUs_In <- qread("output/spatial_units/sus_In.qs")
ZStats_Woody_In <- qread("output/data/ZStats_Woody_In.qs")
ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_In)

ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:3){
    Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_In$PModel$dic$dic + Model_In$NModel$dic$dic, 3)))
  }
  toc(log = TRUE)
  cat("\n\n")
}
proc.time() - ptm


## Test DIC variation with different initial value input  - This part does not work!!! yet? ----
## Setting custom prior based on Full model - This part does not work!!! yet? ----
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


# Model selections----
# (This part is copied from different scripts for running model selection for each clearing type concurrently)

## Agricultural clearing in each KMR ----
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Ag)

inla.setOption(num.threads="8:1")

ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm

ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm

## Forestry clearing in each KMR ----
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
SUs_Fo <- qread("output/spatial_units/sus_fo.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Fo)


# Check model selection step Compare DIC values for Forward complete and backward complete model selection
## Check BC files for each KMR, 
## Extract lowest dDIC valuw based on DIC list
## Compare that with PModel and Nmodel DIC values
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
KMRs <- names(SUs_Ag)
rm(SUs_Ag)

SelModed_BC_Fl <- list.files("output/models/", pattern = glob2rx("SelModel_*_Ag_BC.qs"), full.names = TRUE)

for(kmr in KMRs){
  SelModel <- qread(paste0("output/models/SelModel_", kmr, "_Ag_BC.qs"))
  MinDIC <- SelModel$DIC_ls[which.min(unlist(SelModel$DIC_ls) + sapply(strsplit(names(SelModel$DIC_ls), "\\+"), length)*2)]
  SelModelDIC <- SelModel$PModel$dic$dic + SelModel$NModel$dic$dic
  print(kmr)
  print(MinDIC)
  print(SelModelDIC)
  cat("\n\n")
}

## Check FC files for each KMR

for(kmr in KMRs){
  SelModel <- qread(paste0("output/models/SelModel_", kmr, "_Ag_FC.qs"))
  MinDIC <- SelModel$DIC_ls[which.min(unlist(SelModel$DIC_ls) + sapply(strsplit(names(SelModel$DIC_ls), "\\+"), length)*2)]
  SelExplV <- names(SelModel$CovsCD[[kmr]])
  SelModelDIC <- SelModel$PModel$dic$dic + SelModel$NModel$dic$dic
  print(kmr)
  print(MinDIC)
  print(SelModelDIC)
  print(SelExplV)
  cat("\n\n")
}

# Compare DIC values for Forward complete and backward complete model selection

for(kmr in KMRs){
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_FC.qs"))
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls) + sapply(strsplit(names(SelModel_BC$DIC_ls), "\\+"), length)*2)]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls) + sapply(strsplit(names(SelModel_FC$DIC_ls), "\\+"), length)*2)]
  MinDIC_BC_expl <- unlist(strsplit(names(MinDIC_BC), "\\+")) %>% map( ~str_trim(.)) %>% unlist() %>% sort() %>% paste(collapse = " + ")
  MinDIC_FC_expl <- unlist(strsplit(names(MinDIC_FC), "\\+")) %>% map( ~str_trim(.)) %>% unlist() %>% sort() %>% paste(collapse = " + ")
  SelModelDIC_BC <- SelModel_BC$PModel$dic$dic + SelModel_BC$NModel$dic$dic
  SelModelDIC_FC <- SelModel_FC$PModel$dic$dic + SelModel_FC$NModel$dic$dic
  cat(kmr, ":\n")
  print(MinDIC_BC)
  print(MinDIC_FC)
  print(MinDIC_BC_expl)
  print(MinDIC_FC_expl)
  cat("\n\n")
}

# CC :
#   $`PopDen + ScEc_PC2 + ScEc_PC4 + DistCity + PropVal + Soil_PC1 + slope + Precip + Temp + EcolCond + Area + NatVegReg + LandUse + Fire`
# [1] 6592.657
# 
# $`LandUse + Area + EcolCond + PopDen + Precip + PropVal + Temp + NatVegReg + Soil_PC1 + slope + ScEc_PC2 + ScEc_PC4 + DistCity + Fire`
# [1] 6592.656
# 
# [1] "Area + DistCity + EcolCond + Fire + LandUse + NatVegReg + PopDen + Precip + PropVal + ScEc_PC2 + ScEc_PC4 + slope + Soil_PC1 + Temp"
# [1] "Area + DistCity + EcolCond + Fire + LandUse + NatVegReg + PopDen + Precip + PropVal + ScEc_PC2 + ScEc_PC4 + slope + Soil_PC1 + Temp"
# 
# 
# CST :
#   $`PopDen + DistRoad + AgProf + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Precip + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire`
# [1] 21126.37
# 
# $`Area + EcolCond + Fire + NatVegReg + PopDen + LandTen + Soil_PC1 + LandUse + Precip + Soil_PC2 + Soil_PC3 + slope + DistRoad + AgProf`
# [1] 21126.36
# 
# [1] "AgProf + Area + DistRoad + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + slope + Soil_PC1 + Soil_PC2 + Soil_PC3"
# [1] "AgProf + Area + DistRoad + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + slope + Soil_PC1 + Soil_PC2 + Soil_PC3"
# 
# 
# DRP :
#   $`PopDen + ScEc_PC4 + DistRoad + DistCity + Soil_PC1 + Soil_PC2 + Soil_PC3 + Precip + EcolCond + Area + LandTen + LandUse`
# [1] 5999.661
# 
# $`Area + LandUse + EcolCond + Soil_PC1 + LandTen + PopDen + Soil_PC3 + DistRoad + DistCity + ScEc_PC4 + ScEc_PC2 + Precip`
# [1] 5999.574
# 
# [1] "Area + DistCity + DistRoad + EcolCond + LandTen + LandUse + PopDen + Precip + ScEc_PC4 + Soil_PC1 + Soil_PC2 + Soil_PC3"
# [1] "Area + DistCity + DistRoad + EcolCond + LandTen + LandUse + PopDen + Precip + ScEc_PC2 + ScEc_PC4 + Soil_PC1 + Soil_PC3"
# 
# 
# FW :
#   $`PopDen + ScEc_PC1 + DistRoad + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + EcolCond + Area + LandTen + LandUse`
# [1] 15213.06
# 
# $`Area + LandUse + Soil_PC3 + LandTen + EcolCond + PopDen + Soil_PC1 + slope + DistRoad + ScEc_PC1 + Soil_PC2`
# [1] 15213.06
# 
# [1] "Area + DistRoad + EcolCond + LandTen + LandUse + PopDen + ScEc_PC1 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3"
# [1] "Area + DistRoad + EcolCond + LandTen + LandUse + PopDen + ScEc_PC1 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3"
# 
# 
# NC :
#   $`PopDen + ScEc_PC1 + AgProf + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire`
# [1] 18995.14
# 
# $`Area + LandUse + EcolCond + NatVegReg + LandTen + Soil_PC1 + PopDen + Soil_PC3 + Soil_PC2 + slope + AgProf + Temp + Fire + ScEc_PC4 + ScEc_PC2`
# [1] 18995.77
# 
# [1] "AgProf + Area + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + ScEc_PC1 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3 + Temp"
# [1] "AgProf + Area + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + ScEc_PC2 + ScEc_PC4 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3 + Temp"
# 
# 
# NT :
#   $`PopDen + ScEc_PC5 + DistRoad + DistCity + Soil_PC1 + Soil_PC2 + slope + Precip + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire`
# [1] 9147.844
# 
# $`Area + EcolCond + NatVegReg + PopDen + LandTen + DistCity + Fire + Soil_PC1 + Precip + LandUse + slope + Soil_PC2 + DistRoad + Temp + ScEc_PC5`
# [1] 9147.845
# 
# [1] "Area + DistCity + DistRoad + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + ScEc_PC5 + slope + Soil_PC1 + Soil_PC2 + Temp"
# [1] "Area + DistCity + DistRoad + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + ScEc_PC5 + slope + Soil_PC1 + Soil_PC2 + Temp"
# 
# 
# NS :
#   $`PopDen + ScEc_PC3 + DistRoad + DistCity + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse`
# [1] 20845.19
# 
# $`Area + EcolCond + Soil_PC1 + LandTen + Soil_PC2 + LandUse + NatVegReg + DistRoad + Soil_PC3 + DistCity + PopDen + ScEc_PC3 + slope + Precip + Temp`
# [1] 20843.28
# 
# [1] "Area + DistCity + DistRoad + EcolCond + LandTen + LandUse + NatVegReg + PopDen + ScEc_PC3 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3 + Temp"
# [1] "Area + DistCity + DistRoad + EcolCond + LandTen + LandUse + NatVegReg + PopDen + Precip + ScEc_PC3 + slope + Soil_PC1 + Soil_PC2 + Soil_PC3 + Temp"
# 
# 
# R :
#   $`PopDen + ScEc_PC3 + DistRoad + PropVal + Soil_PC3 + slope + Precip + EcolCond + Area + LandTen + LandUse`
# [1] 8604.176
# 
# $`LandUse + Area + EcolCond + Soil_PC3 + DistRoad + PopDen + slope + Precip + LandTen + DistCity + ScEc_PC5 + PropVal`
# [1] 8602.907
# 
# [1] "Area + DistRoad + EcolCond + LandTen + LandUse + PopDen + Precip + PropVal + ScEc_PC3 + slope + Soil_PC3"
# [1] "Area + DistCity + DistRoad + EcolCond + LandTen + LandUse + PopDen + Precip + PropVal + ScEc_PC5 + slope + Soil_PC3"
# 
# 
# SC :
#   $`PopDen + ScEc_PC3 + PropVal + Precip + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire`
# [1] 3374.913
# 
# $`LandUse + Area + NatVegReg + EcolCond + Fire + LandTen + Precip + PopDen + PropVal + ScEc_PC3 + Soil_PC3`
# [1] 3372.506
# 
# [1] "Area + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + PropVal + ScEc_PC3"
# [1] "Area + EcolCond + Fire + LandTen + LandUse + NatVegReg + PopDen + Precip + PropVal + ScEc_PC3 + Soil_PC3"

# # BEst Forward selection model
# # $`PopDen + LandUse + EcolCond + NatVegReg + Area + LandTen + DistCity + slope + Soil_PC1 + ScEc_PC4 + PropVal + Soil_PC3 + ScEc_PC3 + Temp`
# # [1] 13773.93
# 
# # P ~  ScEc_PC2  + ScEc_PC5 + DistRoad + Soil_PC1 + Soil_PC2 + slope +      Precip  + EcolCond + Area + LandTen + NatVegReg  +      Fire + f(SA1ID, model = "bym", graph = AdjP, scale.model = TRUE)
# 
# F_test2 <- INLA_with_TimeLimit(TimeLimit = 1000, ForN_H1, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
# F_test1 <- INLA_with_Retry(N_retry = 3, Initial_Tlimit = 1000, ForN_H1, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
# 
# ptm <- proc.time()
# Test_5 <- Select_model(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, Direction = "Backward", Verbose = FALSE)
# proc.time() - ptm
# 
# ptm <- proc.time()
# Test_6 <- Select_model(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, Direction = "Complete Forward", Verbose = FALSE, OutputDir = "output/models/")
# proc.time() - ptm
# 
# ptm <- proc.time()
# Test_7 <- Select_model(KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "Forward Complete", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# proc.time()-ptm
# 
# Test_7$Best_DIC_ls
# 
# ptm <- proc.time()
# Test_8 <- Select_model(KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "Backward Complete", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# proc.time()-ptm
# 
# Test_8$Best_DIC_ls
# 
# tic()
# for(i in 1:20){
#   Mod <- fit_model2 (KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
#   print(paste0("DIC= " , Mod$PModel$dic$dic+Mod$NModel$dic$dic))
# }
# toc()
# 

# 
# summary(ResultN)
# 
# 
# Test_3A <- fit_model2(KMR = "CC", ClearType = 1, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "PopDen + ScEc_PC1 + ScEc_PC2 + ScEc_PC3 + ScEc_PC4 + ScEc_PC5 +      DistRoad + DistCity + PropVal + Soil_PC1 + Soil_PC2 + Soil_PC3 +      slope + Precip + Temp + EcolCond + Area + LandTen + NatVegReg +      LandUse + Fire", Verbose = FALSE)
# P ~ PopDen + ScEc_PC1 + ScEc_PC2 + ScEc_PC3 + ScEc_PC4 + ScEc_PC5 + DistRoad + DistCity + PropVal + AgProf + Soil_PC1 + Soil_PC2 + Soil_PC3 + slope + Precip + Temp + EcolCond + Area + LandTen + NatVegReg + LandUse + Fire + f(SA1ID, model = "bym", graph = AdjP,scale.model = TRUE)

# Model predictions ----

# JRR EXAMPLE FOR MAKING PREDICTIONS

# fit model
Fit <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s)

# get predictions from model
Preds <- predict_model(Fit)

# visualise predictions
plot(Preds$Layer)

# save predictions
st_write(Preds$Layer, "output/predictions/test.shp", delete_layer = TRUE)

# Model prediction across all KMRs ----
## Agri----
SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
KMRs <- names(SUs_Ag)
rm(SUs_Ag)

for (kmr in KMRs){
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_FC.qs"))
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls) + sapply(strsplit(names(SelModel_BC$DIC_ls), "\\+"), length)*2)]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls) + sapply(strsplit(names(SelModel_FC$DIC_ls), "\\+"), length)*2)]
  MinDIC_BC_dDIC <- unlist(MinDIC_BC) + sapply(strsplit(names(MinDIC_BC), "\\+"), length)*2
  MinDIC_FC_dDIC <- unlist(MinDIC_FC) + sapply(strsplit(names(MinDIC_FC), "\\+"), length)*2
  Best_Mod <- if_else(MinDIC_BC_dDIC < MinDIC_FC_dDIC, "BC", "FC")
  cat("Best Model for  ", kmr, ": ", Best_Mod, "\n")
  Preds <- predict_model(qread(paste0("output/models/SelModel_", kmr, "_Ag_", Best_Mod, ".qs")))
  qsave(Preds, file = file.path(paste0("output/predictions/Preds_", kmr, "_Ag",".qs")), preset = "fast")
}

SelModel_FW_Ag_FC <- qread("output/models/SelModel_FW_Ag_FC.qs")
SelModel_FW_Ag_FC_SU <- SelModel_FW_Ag_FC$SpatUnits$FW
st_write(SelModel_FW_Ag_FC_SU, "output/diagnosis/SelModel_FW_Ag_FC_SU.shp", delete = TRUE)
st_write(Layer, "output/diagnosis/Pred_FW_Ag.shp", delete = TRUE)
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Ag.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  Preds <- qread(Pred_Flist[i])
  st_write(Preds$Layer, file.path(paste0("output/predictions/Preds_", KMRs[i], "_Ag",".shp")), delete_layer = TRUE)
  assign(paste0("Preds_", KMRs[i], "Ag"), qread(Pred_Flist[i]))
}

Pred_shp_Files <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Ag.shp"), full.names = TRUE)
for(i in 1:length(Pred_shp_Files)){
  assign(paste0("Preds_", KMRs[i], "_Ag"), st_read(Pred_shp_Files[i]))
}
gc()

Pred_Ag <- rbind(Preds_CC_Ag, Preds_CST_Ag, Preds_DRP_Ag, Preds_FW_Ag, Preds_NC_Ag, Preds_NT_Ag, Preds_NS_Ag, Preds_R_Ag, Preds_SC_Ag)
st_write(Pred_Ag, "output/predictions/Pred_Ag.shp", delete_layer = TRUE)

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
