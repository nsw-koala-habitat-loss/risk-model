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

# load Pre-processed data----
ZStats_Woody <- qread("output/data/ZStats_Woody.qs")
ZStats_CovsC <- qread("output/data/ZStats_CovsC.qs")
ZStats_CovsD <- qread("output/data/ZStats_CovsD.qs")
SUs <- qread("output/spatial_units/sus.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")

# # Load proposed covariates based on workshops from lookup xlsx
# CovLookup <- readxl::read_xlsx("Input/covariates/covariate_description.xlsx", sheet = "AllLyr")

# add area to the continuous covariates and then rescale all cont. covariates to have mean of zero and standard deviation of one except for PCA covariates and Koala Habitat suitability
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

# # Export NA SUs shapefile for checking
# SUs_Covs_All <- cbind(do.call(rbind, SUs), do.call(rbind, ZStats_Covs)) 
# SUs_Covs_All_Ag <- SUs_Covs_All %>% 
#   dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, LandTen, NatVegReg, LandUse, Fire) %>% 
#   filter_all(any_vars(is.na(.)))
# st_write(SUs_Covs_All_Ag, "output/diagnosis/SUs_Covs_All_Ag_NA.shp", delete = TRUE)


# CHECKING FOR MULTI-COLLINEARITY ---- 
## across all clearing types and all KMRs

## Continuous Covariates ----
ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC)
Corr_Cont <- cor(ZStats_CovsC_all, use = "complete.obs")
Corr_Cont_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+  ## highlight variables with correlations  > 0.5 OR < -0.5
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
ggsave(Corr_Cont_all_plot, file = "output/collinearity/Corr_Cont_all_plot.png", width = 2000, height = 2000, units = "px")

## Discrete Covariates ----
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
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+   ## highlight variables with correlations  > 0.5 OR < -0.5
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
ggsave(Corr_Categ_all_plot, file = "output/collinearity/Corr_Categ_all_plot.png", width = 2000, height = 2000, units = "px")

## Remove covariates that are collinear ----
### Remove variables with correlations  > 0.6 OR < -0.6
ZStats_Woody <- qread("output/data/ZStats_Woody.qs")
ZStats_Covs <- qread("output/data/ZStats_Covs.qs")
for (i in names(ZStats_Covs)) {
    ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% dplyr::select(-c(Elev, ForType))
}

# Select all covariates for each clearing type----
## Agricultural Clearing ----

SUs <- qread("output/spatial_units/sus.qs")
SUs_Ag <- SUs_ZStats_Ag <- SUs 
ZStats_Woody_Ag <- ZStats_Woody
ZStats_Covs_Ag <- ZStats_Covs
ZStats_Khab_Ag <- list()

for(i in names(SUs_ZStats_Ag)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_Ag[[i]] <- SUs_Ag[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Agricultural Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  ZStats_Covs_Ag[[i]] <- ZStats_Covs_Ag[[i]] %>%
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, LandTen, NatVegReg, LandUse, Drought, Fire)
  
  # Merge SUs, Covs and Woody
  SUs_ZStats_Ag[[i]] <- bind_cols(SUs_Ag[[i]], ZStats_Covs_Ag[[i]], ZStats_Woody_Ag[[i]]) %>% 
    
    # Remove water land use type and NAs in SA1
    filter(LandUse != "5", NatVegReg != "0") %>%
    
    # Remove NAs
    tidyr::drop_na(SA1) %>% tidyr::drop_na(.) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse),
                  NatVegReg = droplevels(NatVegReg))
  
  ## Split SUs, Covs, Woody and KoalaHabitat into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% dplyr::select(names(SUs_Ag[[i]]))
  ## Select covariates for Agricultural Clearing 
  ZStats_Covs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Ag[[i]]))
  ## Select woody extent and loss columns
  ZStats_Woody_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Ag[[i]]))
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select("sum.woody","sum.aloss", "sum.Khab" )
}

### Check for collinearity in Agricultural Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates 
ZStats_CovsC_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)

#### Calculate correlation matrix
Corr_Cont_Ag_all <- cor(ZStats_CovsC_Ag_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Ag_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
# Corr_Cont_Ag_all_plot
ggsave(Corr_Cont_Ag_all_plot, file = "output/collinearity/Corr_Cont_Ag_all_plot.png", width = 2000, height = 2000, units = "px")

#### Discrete Covariates ----
#### Combine all Covariates for all KMRs then select discrete covariates
ZStats_CovsD_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(LandTen, NatVegReg, LandUse, Drought, Fire)

#### Calculate Cramers V matrix
cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_Ag_all), ncol = ncol(ZStats_CovsD_Ag_all), 
                           dimnames = list(names(ZStats_CovsD_Ag_all), names(ZStats_CovsD_Ag_all)))

#### Calculate Cramers V for each pair of discrete covariates
for (i in 1:ncol(ZStats_CovsD_Ag_all)) {
  for(j in 1:ncol(ZStats_CovsD_Ag_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_Ag_all[[i]], ZStats_CovsD_Ag_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

#### Plot Cramers V matrix
Corr_CovD_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
# Corr_CovD_Ag_all_plot
ggsave(Corr_CovD_Ag_all_plot, file = "output/collinearity/Corr_CovD_Ag_all_plot.png", width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Ag, file = "output/spatial_units/SUs_Ag.qs", preset = "fast")
qsave(ZStats_Woody_Ag, file = "output/data/ZStats_Woody_Ag.qs", preset = "fast")
qsave(ZStats_Covs_Ag, file = "output/data/ZStats_Covs_Ag.qs", preset = "fast")
qsave(ZStats_Khab_Ag, file = "output/data/ZStats_Khab_Ag.qs", preset = "fast")

### Check for NAs in the data ----
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

SUs <- qread("output/spatial_units/sus.qs")
SUs_Fo <- SUs_ZStats_Fo <- SUs 
ZStats_Woody_Fo <- ZStats_Woody
ZStats_Covs_Fo <- ZStats_Covs
ZStats_Khab_Fo <- list()

for(i in names(SUs_ZStats_Fo)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_Fo[[i]] <- SUs_Fo[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Forestry Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  # Drop LandTen as the infor captured in TenType and Drop ForTen with limited SUs in each category 
  ZStats_Covs_Fo[[i]] <- ZStats_Covs_Fo[[i]] %>% 
    dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area, NatVegReg, LandUse, Drought, Fire, TenType)
  
  # Merge SUs, Covs and Woody 
  SUs_ZStats_Fo[[i]] <- bind_cols(SUs_Fo[[i]], ZStats_Covs_Fo[[i]], ZStats_Woody_Fo[[i]], ZStats_Khab_Fo[[i]]) %>% 
    
    # Remove water land use type and NAs in SA1
    # Remove SUs in LLS excluded area (non-rural) # Remove Multiple-use public forest(2) and Nature conservation reserve (3) Forest Tenure Type
    filter(LandUse != "5", NatVegReg != "0", TenType == "1") %>% 
    
    # Remove NAs
    tidyr::drop_na(.) %>% tidyr::drop_na(SA1) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse),
                  NatVegReg = droplevels(NatVegReg))
  
  ## Split SUs, Covs, Woody and KoalaHabitat into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% dplyr::select(names(SUs_Fo[[i]]))
  ## Select covariates for Forestry Clearing
  ZStats_Covs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Fo[[i]])) %>% dplyr::select(-TenType)
  ## Select woody extent and loss columns
  ZStats_Woody_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Fo[[i]]))
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Khab_Fo[[i]]))
}

### Check for collinearity in Forestry Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates 
ZStats_CovsC_Fo_all <- do.call(rbind, ZStats_Covs_Fo) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)

#### Calculate correlation matrix
Corr_Cont_Fo_all <- cor(ZStats_CovsC_Fo_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_Fo_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Fo_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
# Corr_Cont_Fo_all_plot
ggsave(Corr_Cont_Fo_all_plot, file = "output/collinearity/Corr_Cont_Fo_all_plot.png", width = 2000, height = 2000, units = "px")

#### Discrete Covariates ----
#### Combine all Covariates for all KMRs then select discrete covariates
ZStats_CovsD_Fo_all <- do.call(rbind, ZStats_Covs_Fo) %>% 
  dplyr::select(NatVegReg, LandUse, Fire, Drought )

#### Calculate Cramers V matrix
cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_Fo_all), ncol = ncol(ZStats_CovsD_Fo_all), 
                           dimnames = list(names(ZStats_CovsD_Fo_all), names(ZStats_CovsD_Fo_all)))

#### Calculate Cramers V for each pair of discrete covariates
for (i in 1:ncol(ZStats_CovsD_Fo_all)) {
  for(j in 1:ncol(ZStats_CovsD_Fo_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_Fo_all[[i]], ZStats_CovsD_Fo_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}

#### Plot Cramers V matrix
Corr_CovD_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
# Corr_CovD_Ag_all_plot
ggsave(Corr_CovD_Ag_all_plot, file = "output/collinearity/Corr_CovD_Ag_all_plot.png", width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Fo, file = "output/spatial_units/sus_fo.qs", preset = "fast")
qsave(ZStats_Woody_Fo, file = "output/data/ZStats_Woody_Fo.qs", preset = "fast")
qsave(ZStats_Covs_Fo, file = "output/data/ZStats_Covs_Fo.qs", preset = "fast")
qsave(ZStats_Khab_Fo, file = "output/data/ZStats_Khab_Fo.qs", preset = "fast")

### Check for NAs in the data ----
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
SUs <- qread("output/spatial_units/sus.qs")
SUs_In <- SUs_ZStats_In <- SUs
ZStats_Woody_In <- ZStats_Woody
ZStats_Covs_In <- ZStats_Covs
ZStats_Khab_In <- list()

for(i in names(SUs_ZStats_In)){
  
  # Remove Area column from SUs before merging with ZStats to prevent duplication
  SUs_In[[i]] <- SUs_In[[i]] %>% dplyr::select(-Area)
  
  # Select covariates for Forestry Clearing (Based on DAGs)
  # Drop PolPref: Many KMR has only 1 Category and causing model fitting issues
  # Drop NatVegReg due as it's correlated to LandUse and PlanZone
  # Drop fire, and Soil due to uncertain direct effect on Infrastructure Clearing
  ZStats_Covs_In[[i]] <- ZStats_Covs_In[[i]] %>% 
    dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, slope, Precip, Temp, EcolCond, Area, LandTen, PlanZone, LandUse, Drought, Remoteness)
  
  # Merge SUs, Covs and Woody 
  SUs_ZStats_In[[i]] <- bind_cols(SUs_In[[i]], ZStats_Covs_In[[i]], ZStats_Woody_In[[i]], ZStats_Khab_In[[i]]) %>% 
    
    # Remove SUs in "0" (Major Cities), "1" (Inner Regional), "2" (Outer Regional)
    filter(Remoteness %in% c("0", "1", "2")) %>% 
    
    # Remove "5" (water) and "1" (Conservation and natural environments) land use type and NAs in SA1
    filter(LandUse != "5") %>% 
    
    # Remove NAs
    tidyr::drop_na(SA1) %>% tidyr::drop_na(.) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse))
    
  
  ## Split SUs, Covs and Woody into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_In[[i]] <- SUs_ZStats_In[[i]] %>% dplyr::select(names(SUs_In[[i]]))
  ## Select covariates for Infrastructure Clearing
  ZStats_Covs_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_In[[i]]))%>% dplyr::select(-Remoteness)
  ## Select woody extent and loss columns
  ZStats_Woody_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_In[[i]])) 
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Khab_In[[i]]))
}

### Check for collinearity in Forestry Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates
ZStats_CovsC_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, slope, Precip, Temp, EcolCond, Area)

#### Calculate correlation matrix
Corr_Cont_In_all <- cor(ZStats_CovsC_In_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_In_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_In_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
# Corr_Cont_In_all_plot
ggsave(Corr_Cont_In_all_plot, file = "output/collinearity/Corr_Cont_In_all_plot.png", width = 2000, height = 2000, units = "px")

#### Discrete Covariates ----
#### Combine all Covariates for all KMRs then select discrete covariates
ZStats_CovsD_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(LandTen, PlanZone, LandUse, Drought)
#### Calculate Cramers V matrix
cramers_v_matrix <- matrix(NA, nrow = ncol(ZStats_CovsD_In_all), ncol = ncol(ZStats_CovsD_In_all), 
                           dimnames = list(names(ZStats_CovsD_In_all), names(ZStats_CovsD_In_all)))
#### Calculate Cramers V for each pair of discrete covariates
for (i in 1:ncol(ZStats_CovsD_In_all)) {
  for(j in 1:ncol(ZStats_CovsD_In_all)) {
    cramersv_val <- cramersv(chisq.test(ZStats_CovsD_In_all[[i]], ZStats_CovsD_In_all[[j]]))
    cramers_v_matrix[i,j] <- cramersv_val
  }
}
#### Plot Cramers V matrix
Corr_CovD_In_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = cramers_v_matrix, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
# Corr_CovD_In_all_plot
ggsave(Corr_CovD_In_all_plot, file = "output/collinearity/Corr_CovD_In_all_plot.png", width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_In, file = "output/spatial_units/sus_In.qs", preset = "fast")
qsave(ZStats_Woody_In, file = "output/data/ZStats_Woody_In.qs", preset = "fast")
qsave(ZStats_Covs_In, file = "output/data/ZStats_Covs_In.qs", preset = "fast")
qsave(ZStats_Khab_In, file = "output/data/ZStats_Khab_In.qs", preset = "fast")


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
# 
# map(ZStats_Woody_In, ~summary(.))
# sum(ZStats_Woody_In$CST$sum.woody>0)
# length(ZStats_Woody_In$CST$sum.woody)
# summary(ZStats_Covs_In$CC$PopGro)
# 
# lm_mod <- lm(ZStats_Covs_In$CC$PopGro ~ 1, na.action = na.omit)
# b <- boxcox(lm_mod, plotit = FALSE)
# ZStats_CovsC_LB[j,i] <- b$x[which.max(b$y)]
# 
# hist_PopGro <- ggplot(ZStats_Covs_In$CC, aes(x = PopGro+1)) +
#   geom_histogram(binwidth = 0.1) +
#   labs(x = "PopGro", y = "Frequency")+
#   theme_pubr()
# hist_PopGro

### Check for NAs in the data ----
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

## Agricultural ----
## Load data
# ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
# ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
# SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
# SA1s <- qread("output/spatial_units/sa1s.qs")
# # Test <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s)
# TEST2 <- fit_model2(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# 


## Check DIC variations ----
### Agricultural clearing in each KMR ----

ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Ag)

ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:5){
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
# summary(ZStats_Woody_Fo$FW$sum.woody)
ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:5){
    Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_Fo$PModel$dic$dic + Model_Fo$NModel$dic$dic, 2)))
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
inla.setOption(num.threads="8:1")
# map(ZStats_Covs_In, ~summary(.))

ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:5){
    Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_In$PModel$dic$dic + Model_In$NModel$dic$dic, 2)))
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


# MODEL SELECTIONS----
### This part can be run independently without running the model fitting steps above
### Each clearing types and selection directions was set to run independently and can be run in parallel
### Adjust inla.setOption(num.threads= "8:1") to allocate the number of threads to be used in 1 inla run. The number need to be reduce if running parallel.

## Agricultural----
## Load data
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Ag)

# inla.setOption(num.threads="8:1")

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

## Forestry----
## Load data
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
SUs_Fo <- qread("output/spatial_units/sus_Fo.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_Fo)

# inla.setOption(num.threads="8:1")

## Forward model selection
ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm

## Backward model selection
ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm

## Infrastructure----
## Load data
ZStats_Woody_In <- qread("output/data/ZStats_Woody_In.qs")
ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")
SUs_In <- qread("output/spatial_units/sus_In.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
KMRs <- names(SUs_In)

# inla.setOption(num.threads="8:1")
## Forward model selection
ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm

## Backward model selection
ptm <- proc.time()
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = "output/models/")
}
proc.time() - ptm



## Check model selection step ---- 
### Check errors in Model selections steps ---
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
KMRs <- names(SUs_Ag)
rm(SUs_Ag)
ClrTyps <- c("Ag", "Fo", "In")

for(ClrTyp in ClrTyps){
  for(kmr in KMRs){
    SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_" , ClrTyp,"_FC.qs"))
    SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_" , ClrTyp,"_FC.qs"))
    if(length(SelModel_FC$ERROR_ls) == 0){
      cat(ClrTyp, kmr, " Forward selection: NO ERROR! :-)", "\n")
    }else if(length(SelModel_FC$ERROR_ls) > 0){
      cat(ClrTyp, kmr, " Forward selection error: ","\n")
      print(SelModel_FC$ERROR_ls)
      cat("\n\n")}
    if(length(SelModel_BC$ERROR_ls) == 0){
      cat(ClrTyp, kmr, " Backward selection: NO ERROR! :-)", "\n")
    }else if(length(SelModel_BC$ERROR_ls) > 0){
      cat(ClrTyp, kmr, " Backward selection error: ", "\n")
      print(SelModel_BC$ERROR_ls)
      cat("\n\n")}
  }}


### Compare DIC values for Forward complete and backward complete model selection ----
## Also verify the minimum DIC values recorded in the model selection steps (WHILE LOOP) and the final selected model  

SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
KMRs <- names(SUs_Ag)
rm(SUs_Ag)
ClrTyps <- c("Ag", "Fo", "In")

for(ClrTyp in ClrTyps){  ## Agriculture, Forestry, Infrastructure
  for(kmr in KMRs){  ## KMRs
    
    # read model selection results for both forward and backward selection
    SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_" , ClrTyp , "_BC.qs"))
    SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_" , ClrTyp , "_FC.qs"))
    
    # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
    MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
    MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
    MinDIC_BC_expl <- unlist(strsplit(names(MinDIC_BC), "\\+")) %>% map( ~str_trim(.)) %>% unlist() %>% sort() %>% paste(collapse = " + ")
    MinDIC_FC_expl <- unlist(strsplit(names(MinDIC_FC), "\\+")) %>% map( ~str_trim(.)) %>% unlist() %>% sort() %>% paste(collapse = " + ")
    
    # Get DIC values for the selected model
    SelModelDIC_BC <- SelModel_BC$PModel$dic$dic + SelModel_BC$NModel$dic$dic
    SelModelDIC_FC <- SelModel_FC$PModel$dic$dic + SelModel_FC$NModel$dic$dic
    SelModel_BC_Covs <- names(SelModel_BC$CovsCD[[kmr]]) %>% sort() %>% paste(collapse = " + ")
    SelModel_FC_Covs <- names(SelModel_FC$CovsCD[[kmr]]) %>% sort() %>% paste(collapse = " + ")
    
    # Print results
    cat(ClrTyp, kmr, ":\n")
    
    # Backward complete model selection
    cat("BC:\n")
    cat("DIC_ls: Min DIC:", unlist(MinDIC_BC), "\n")
    cat(MinDIC_BC_expl, "\n")
    cat("Sel_Model: DIC:", SelModelDIC_BC, "\n" )
    cat(SelModel_BC_Covs, "\n\n")
    
    if(MinDIC_BC_expl != SelModel_BC_Covs){message("Selected model does not match the minimum DIC model!!! \n\n")}
    
    # Forward complete model selection
    cat("FC:\n")
    cat("DIC_ls: Min DIC:", unlist(MinDIC_FC), "\n")
    cat(MinDIC_FC_expl, "\n")
    cat("Sel_Model: DIC:", SelModelDIC_FC, "\n")
    cat(SelModel_FC_Covs, "\n\n")
    
    if(MinDIC_FC_expl != SelModel_FC_Covs){message("Selected model does not match the minimum DIC model!!! \n\n")}
    
    # Check if the selected model is the same for both forward and backward selection
    if(MinDIC_BC_expl == MinDIC_FC_expl){
      cat("Forward and backward selection produce the same result! NICE!!\n\n")
    }else{
      message("Forward and backward selection produce different result!!! \n\n")
    }
  }}

## Export model selection results and coefficient estimates ----

### Agriculture ----
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
KMRs <- names(ZStats_Covs_Ag)
kmr <- KMRs[1]
Model_Ag <- fit_model2(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_Ag <- summary(Model_Ag$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="CC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  if(MinDIC_BC_DIC < MinDIC_FC_DIC){
    Best_Mod <- SelModel_BC} else{
      Best_Mod <- SelModel_FC}
  cat("Best Model for  ", kmr, ": ", if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC"), "\n")
  
  # Get the covariate of the selected model
  Cov <- summary(Best_Mod$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
  
  # Get the coefficient of the selected model
  Cof_PModel <- summary(Best_Mod$PModel)$fixed[,1]
  Cof_PModel <- if_else(Cof_PModel>0, "+", "-")
  Cof_NModel <- summary(Best_Mod$NModel)$fixed[,1]
  Cof_NModel <- if_else(Cof_NModel>0, "+", "-")
  Cof <- paste(Cof_PModel, Cof_NModel, sep = " / ")
  Cov_cof <- as.data.frame(cbind(Cov, Cof))
  rownames(Cov_cof) <- NULL
  names(Cov_cof) <- c("Cov", kmr)
  Cov_ls_Ag <- left_join(Cov_ls_Ag, Cov_cof, by  = join_by("Covariate" == "Cov"))
}

# Prepare the table for the selected model
Cov_ls_Ag_tab <- Cov_ls_Ag %>% filter_all(any_vars(!is.na(.))) %>%
  mutate(Variable = c("(Intercept)", "Population density", "Socio-Economic PC1", "Socio-Economic PC2", "Socio-Economic PC3", "Socio-Economic PC4",  "Socio-Economic PC5", 
                      "Distance to road", "Distance to urban centre", "Property value", "Agricultural profit", "Soil PC1", "Soil PC2", "Soil PC3", 
                      "Slope", "Rainfall", "Temperature", "Ecological condition", "Property size", "Land Tenure2", "Land Tenure3", "Land Tenure4",
                      "Land use regulations", "Land use type1", "Land use type2", "Land use type3", "Land use type4", "Drought", "Fire")) %>% 
  dplyr::select(Variable, everything()) %>% 
  dplyr::select(-Covariate) %>% 
  arrange(Variable) %>% 
  filter(Variable != "Land use type5",
         !(is.na(CC) & is.na(CST) & is.na(DRP) & is.na(FW) & is.na(NC) & is.na(NT) & is.na(NS) & is.na(R) & is.na(SC) ))
# Export the table
write.csv(Cov_ls_Ag_tab, "output/models/Cov_ls_Ag.csv", row.names = FALSE, na = "")

### Forestry ----
SUs_Fo <- qread("output/spatial_units/SUs_Fo.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
KMRs <- names(ZStats_Covs_Fo)
kmr <- KMRs[1]
Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_Fo <- summary(Model_Fo$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="CC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  if(MinDIC_BC_DIC < MinDIC_FC_DIC){
    Best_Mod <- SelModel_BC} else{
      Best_Mod <- SelModel_FC}
  cat("Best Model for  ", kmr, ": ", if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC"), "\n")
  
  # Get the covariate of the selected model
  Cov <- summary(Best_Mod$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
  
  # Get the coefficient of the selected model
  Cof_PModel <- summary(Best_Mod$PModel)$fixed[,1]
  Cof_PModel <- if_else(Cof_PModel>0, "+", "-")
  Cof_NModel <- summary(Best_Mod$NModel)$fixed[,1]
  Cof_NModel <- if_else(Cof_NModel>0, "+", "-")
  Cof <- paste(Cof_PModel, Cof_NModel, sep = " / ")
  Cov_cof <- as.data.frame(cbind(Cov, Cof))
  rownames(Cov_cof) <- NULL
  names(Cov_cof) <- c("Cov", kmr)
  Cov_ls_Fo <- left_join(Cov_ls_Fo, Cov_cof, by  = join_by("Covariate" == "Cov"))
}

Cov_ls_tab <- Cov_ls_Fo %>% 
  mutate(Variable = c("(Intercept)", "Population density", "Socio-Economic PC1", "Socio-Economic PC2", "Socio-Economic PC3", "Socio-Economic PC4",  "Socio-Economic PC5", 
                      "Distance to road", "Distance to urban centre", "Property value", "Agricultural profit", "Soil PC1", "Soil PC2", "Soil PC3", 
                      "Slope", "Rainfall", "Temperature", "Ecological condition", "Property size", "Land use regulations",
                      "Land use type1", "Land use type2", "Land use type3", "Land use type4", "Drought", "Fire")) %>%
  dplyr::select(Variable, everything()) %>%
  dplyr::select(-Covariate) %>%
  arrange(Variable) %>%
  filter(Variable != "Land use type5",
         !(is.na(CC) & is.na(CST) & is.na(DRP) & is.na(FW) & is.na(NC) & is.na(NT) & is.na(NS) & is.na(R) & is.na(SC) ))
Cov_ls_tab
write.csv(Cov_ls_tab, "output/models/Cov_ls_Fo.csv", row.names = FALSE, na = "")

### Infrastructure ----
SUs_In <- qread("output/spatial_units/SUs_In.qs")
SA1s <- qread("output/spatial_units/sa1s.qs")
ZStats_Woody_In <- qread("output/data/ZStats_Woody_In.qs")
ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")
KMRs <- names(ZStats_Covs_In)
kmr <- KMRs[1]
Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_In <- summary(Model_In$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="CC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_In_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_In_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  if(MinDIC_BC_DIC < MinDIC_FC_DIC){
    Best_Mod <- SelModel_BC} else{
      Best_Mod <- SelModel_FC}
  cat("Best Model for  ", kmr, ": ", if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC"), "\n")
  
  # Get the covariate of the selected model
  Cov <- summary(Best_Mod$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
  
  # Get the coefficient of the selected model
  Cof_PModel <- summary(Best_Mod$PModel)$fixed[,1]
  Cof_PModel <- if_else(Cof_PModel>0, "+", "-")
  Cof_NModel <- summary(Best_Mod$NModel)$fixed[,1]
  Cof_NModel <- if_else(Cof_NModel>0, "+", "-")
  Cof <- paste(Cof_PModel, Cof_NModel, sep = " / ")
  Cov_cof <- as.data.frame(cbind(Cov, Cof))
  rownames(Cov_cof) <- NULL
  names(Cov_cof) <- c("Cov", kmr)
  Cov_ls_In <- left_join(Cov_ls_In, Cov_cof, by  = join_by("Covariate" == "Cov"))
}

Cov_ls_tab <- Cov_ls_In %>%
  mutate(Variable = c("(Intercept)", "Population density", "Population growth", "Socio-Economic PC1", "Socio-Economic PC2", "Socio-Economic PC3", "Socio-Economic PC4",  "Socio-Economic PC5", 
                      "Distance to road", "Distance to urban centre", "Property value", "Slope", "Rainfall", "Temperature", "Ecological condition", "Property size", "Land Tenure2", "Land Tenure3", "Land Tenure4",
                      "Planning Zone1", "Planning Zone2", "Planning Zone3", "Planning Zone4", "Land use type1", "Land use type2", "Land use type3", "Land use type4", "Drought")) %>%
  dplyr::select(Variable, everything()) %>%
  dplyr::select(-Covariate) %>%
  arrange(Variable) %>% 
  filter(Variable != "Land use type5",
                            !(is.na(CC) & is.na(CST) & is.na(DRP) & is.na(FW) & is.na(NC) & is.na(NT) & is.na(NS) & is.na(R) & is.na(SC) ))

write.csv(Cov_ls_tab, "output/models/Cov_ls_In.csv", row.names = FALSE, na = "")

# MODEL PREDICTIONS ----
# JRR EXAMPLE FOR MAKING PREDICTIONS
# 
# # fit model
# Fit <- fit_model(KMR = "CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s)
# 
# # get predictions from model
# Preds <- predict_model(Fit)
# 
# # visualise predictions
# plot(Preds$Layer)
# 
# # save predictions
# st_write(Preds$Layer, "output/predictions/test.shp", delete_layer = TRUE)

# Model prediction across all KMRs ----
KMR_shp <- st_read("input/spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")

## Agriculture ----
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
KMRs <- names(SUs_Ag)
rm(SUs_Ag)

MindDIC_ls <- list()

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Ag_FC.qs"))
 
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  Best_Mod <- if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC")
  cat("Best Model for  ", kmr, ": ", Best_Mod, "\n")
  Preds <- predict_model2(qread(paste0("output/models/SelModel_", kmr, "_Ag_", Best_Mod, ".qs")))
  qsave(Preds, file = file.path(paste0("output/predictions/Preds_", kmr, "_Ag",".qs")), preset = "fast")
}

# 
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Ag.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  Preds <- qread(Pred_Flist[i])
  st_write(Preds$Layer, file.path(paste0("output/predictions/Preds_", KMRs[i], "_Ag",".shp")), delete_layer = TRUE)
}

Pred_shp_Files <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Ag.shp"), full.names = TRUE)
for(i in 1:length(Pred_shp_Files)){
  assign(paste0("Preds_", KMRs[i], "_Ag"), st_read(Pred_shp_Files[i]))
}

Pred_Ag <- rbind(Preds_CC_Ag, Preds_CST_Ag, Preds_DRP_Ag, Preds_FW_Ag, Preds_NC_Ag, Preds_NT_Ag, Preds_NS_Ag, Preds_R_Ag, Preds_SC_Ag)

st_write(Pred_Ag, "output/predictions/Pred_Ag.shp", delete_layer = TRUE)


Pred_Ag <- st_read("output/predictions/Pred_Ag.shp")
Pred_Ag_plot <- ggplot(Pred_Ag)+
  geom_sf(aes(fill = PredAll), lwd = 0)+
  scale_fill_viridis_c(direction = -1, option = "A")+
  geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.1)+
  theme_void()
ggsave(filename = "output/figures/Pred_Ag_plot.png", plot = Pred_Ag_plot, units = "cm", width = 15, height = 15, dpi = 300)

## Forestry ----
SUs_Fo <- qread("output/spatial_units/SUs_Fo.qs")
KMRs <- names(SUs_Fo)
rm(SUs_Fo)

MindDIC_ls <- list()

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  Best_Mod <- if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC")
  cat("Best Model for  ", kmr, ": ", Best_Mod, "\n")
  Preds <- predict_model2(qread(paste0("output/models/SelModel_", kmr, "_Fo_", Best_Mod, ".qs")))
  qsave(Preds, file = file.path(paste0("output/predictions/Preds_", kmr, "_Fo",".qs")), preset = "fast")
}

# 
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Fo.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  Preds <- qread(Pred_Flist[i])
  st_write(Preds$Layer, file.path(paste0("output/predictions/Preds_", KMRs[i], "_Fo",".shp")), delete_layer = TRUE)
}

Pred_shp_Files <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Fo.shp"), full.names = TRUE)
for(i in 1:length(Pred_shp_Files)){
  assign(paste0("Preds_", KMRs[i], "_Fo"), st_read(Pred_shp_Files[i]))
}

Pred_Fo <- rbind(Preds_CC_Fo, Preds_CST_Fo, Preds_DRP_Fo, Preds_FW_Fo, Preds_NC_Fo, Preds_NT_Fo, Preds_NS_Fo, Preds_R_Fo, Preds_SC_Fo)

st_write(Pred_Fo, "output/predictions/Pred_Fo.shp", delete_layer = TRUE)


Pred_Fo <- st_read("output/predictions/Pred_Fo.shp")
Pred_Fo_plot <- ggplot(Pred_Fo)+
  geom_sf(aes(fill = PredAll), lwd = 0)+
  scale_fill_viridis_c(direction = -1, option = "A")+
  geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.1)+
  theme_void()
ggsave(filename = "output/figures/Pred_Fo_plot.png", plot = Pred_Fo_plot, units = "cm", width = 15, height = 15, dpi = 300)


## Infrastructure ----
SUs_In <- qread("output/spatial_units/SUs_In.qs")
KMRs <- names(SUs_In)
rm(SUs_In)

MindDIC_ls <- list()

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_In_BC.qs"))
  SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_In_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  Best_Mod <- if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC")
  cat("Best Model for  ", kmr, ": ", Best_Mod, "\n")
  Preds <- predict_model2(qread(paste0("output/models/SelModel_", kmr, "_In_", Best_Mod, ".qs")))
  qsave(Preds, file = file.path(paste0("output/predictions/Preds_", kmr, "_In",".qs")), preset = "fast")
}

# 
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_In.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  Preds <- qread(Pred_Flist[i])
  st_write(Preds$Layer, file.path(paste0("output/predictions/Preds_", KMRs[i], "_In",".shp")), delete_layer = TRUE)
}

Pred_shp_Files <- list.files("output/predictions/", pattern = glob2rx("Preds_*_In.shp"), full.names = TRUE)
for(i in 1:length(Pred_shp_Files)){
  assign(paste0("Preds_", KMRs[i], "_In"), st_read(Pred_shp_Files[i]))
}

Pred_In <- rbind(Preds_CC_In, Preds_CST_In, Preds_DRP_In, Preds_FW_In, Preds_NC_In, Preds_NT_In, Preds_NS_In, Preds_R_In, Preds_SC_In)

st_write(Pred_In, "output/predictions/Pred_In.shp", delete_layer = TRUE)


Pred_In <- st_read("output/predictions/Pred_In.shp")
Pred_In_plot <- ggplot(Pred_In)+
  geom_sf(aes(fill = PredAll), lwd = 0)+
  scale_fill_viridis_c(direction = -1, option = "A")+
  geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.1)+
  theme_void()
ggsave(filename = "output/figures/Pred_In_plot.png", plot = Pred_In_plot, units = "cm", width = 15, height = 15, dpi = 300)

# Koala habitat loss risk ----

## Agriculture ----
### Load SU Zonal data for Woody, Koala habitat and Predicted Risk 
ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
ZStats_Khab_Ag <- qread("output/data/ZStats_Khab_Ag.qs")
KMRs <- names(ZStats_Woody_Ag)
kmr=KMRs[1]

### load prediction data for each KMR then make lists based on KMR
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Ag.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  assign(basename(tools::file_path_sans_ext(Pred_Flist[i])), qread(Pred_Flist[i]))
}

Pred_Ag_ls <- list("CC" = Preds_CC_Ag$Layer, "CST" = Preds_CST_Ag$Layer, "DRP" = Preds_DRP_Ag$Layer, "FW" = Preds_FW_Ag$Layer, "NC" = Preds_NC_Ag$Layer, "NT" = Preds_NT_Ag$Layer, "NS" = Preds_NS_Ag$Layer, "R" = Preds_R_Ag$Layer, "SC" = Preds_SC_Ag$Layer)
qsave(Pred_Ag_ls, file = "output/predictions/Pred_Ag.qs", preset = "fast")

Khab_Woody_Ag <- ZStats_Khab_Ag
ZStats_WoodyRem_Ag <- ZStats_Woody_Ag
Khab_risk_Ag <- ZStats_Khab_Ag

for(kmr in KMRs){
  
  #### Filter data for properties with woody vegetation in 2011 and ensure clearing is < woody vegetation
  ZStats_WoodyRem_Ag[[kmr]] <- ZStats_WoodyRem_Ag[[kmr]] %>% 
    mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
    mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% 
    dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody, Khab = sum.Khab) %>%
    filter(N > 0)
  # NT <-  Response %>% dplyr::select(N) %>% as.matrix()
  
  # Warning message if woody vegetation is greater than Koala habitat
  if(nrow(ZStats_WoodyRem_Ag[[kmr]][ZStats_WoodyRem_Ag[[kmr]]$N < ZStats_WoodyRem_Ag[[kmr]]$Khab,])>0){message("KMR=", kmr, ": Woody vegetation greater than Koala habitat!!! \n\n")}
  
  ### Koala habitat risk data (Woody, Risk, Khab/Woody, Khab*Risk)
  Khab_risk_Ag[[kmr]] <- bind_cols(Pred_Ag_ls[[kmr]], ZStats_WoodyRem_Ag[[kmr]]) %>%
    mutate(Khab = ifelse(Khab < N, Khab, N),
      Khab_P = Khab/N, KhabRisk = PredAll*Khab_P) %>% 
    dplyr::select(Woody = N, Risk = PredAll, Khab_P, KhabRisk)
}

qsave(Khab_risk_Ag, file = "output/predictions/Khab_risk_Ag.qs", preset = "fast")
Khab_risk_Ag_all <- do.call(rbind, Khab_risk_Ag)
st_write(Khab_risk_Ag_all, "output/predictions/Khab_risk_Ag.shp", delete_layer = TRUE)

## Forestry ----
### Load SU Zonal data for Woody, Koala habitat and Predicted Risk 
ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
ZStats_Khab_Fo <- qread("output/data/ZStats_Khab_Fo.qs")
KMRs <- names(ZStats_Woody_Fo)
kmr=KMRs[1]

### load prediction data for each KMR then make lists based on KMR
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_Fo.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  assign(basename(tools::file_path_sans_ext(Pred_Flist[i])), qread(Pred_Flist[i]))
}

Pred_Fo_ls <- list("CC" = Preds_CC_Fo$Layer, "CST" = Preds_CST_Fo$Layer, "DRP" = Preds_DRP_Fo$Layer, "FW" = Preds_FW_Fo$Layer, "NC" = Preds_NC_Fo$Layer, "NT" = Preds_NT_Fo$Layer, "NS" = Preds_NS_Fo$Layer, "R" = Preds_R_Fo$Layer, "SC" = Preds_SC_Fo$Layer)
qsave(Pred_Fo_ls, file = "output/predictions/Pred_Fo.qs", preset = "fast")

Khab_Woody_Fo <- ZStats_Khab_Fo
ZStats_WoodyRem_Fo <- ZStats_Woody_Fo
Khab_risk_Fo <- ZStats_Khab_Fo

for(kmr in KMRs){
  
  #### Filter data for properties with woody vegetation in 2011 and ensure clearing is < woody vegetation
  ZStats_WoodyRem_Fo[[kmr]] <- ZStats_WoodyRem_Fo[[kmr]] %>% 
    mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
    mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% 
    dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody, Khab = sum.Khab) %>%
    filter(N > 0)
  # NT <-  Response %>% dplyr::select(N) %>% as.matrix()
  
  # Warning message if woody vegetation is greater than Koala habitat
  if(nrow(ZStats_WoodyRem_Fo[[kmr]][ZStats_WoodyRem_Fo[[kmr]]$N < ZStats_WoodyRem_Fo[[kmr]]$Khab,])>0){message("KMR=", kmr, ": Woody vegetation greater than Koala habitat!!! \n\n")}
  
  ### Koala habitat risk data (Woody, Risk, Khab/Woody, Khab*Risk)
  Khab_risk_Fo[[kmr]] <- bind_cols(Pred_Fo_ls[[kmr]], ZStats_WoodyRem_Fo[[kmr]]) %>%
    mutate(Khab = ifelse(Khab < N, Khab, N),
           Khab_P = Khab/N, KhabRisk = PredAll*Khab_P) %>% 
    dplyr::select(Woody = N, Risk = PredAll, Khab_P, KhabRisk)
}

qsave(Khab_risk_Fo, file = "output/predictions/Khab_risk_Fo.qs", preset = "fast")
Khab_risk_Fo_all <- do.call(rbind, Khab_risk_Fo)
st_write(Khab_risk_Fo_all, "output/predictions/Khab_risk_Fo.shp", delete_layer = TRUE)

## Infrastructure ----
### Load SU Zonal data for Woody, Koala habitat and Predicted Risk 
ZStats_Woody_In <- qread("output/data/ZStats_Woody_In.qs")
ZStats_Khab_In <- qread("output/data/ZStats_Khab_In.qs")
KMRs <- names(ZStats_Woody_In)
kmr=KMRs[1]

### load prediction data for each KMR then make lists based on KMR
Pred_Flist <- list.files("output/predictions/", pattern = glob2rx("Preds_*_In.qs"), full.names = TRUE)

for (i in 1:length(Pred_Flist)){
  assign(basename(tools::file_path_sans_ext(Pred_Flist[i])), qread(Pred_Flist[i]))
}

Pred_In_ls <- list("CC" = Preds_CC_In$Layer, "CST" = Preds_CST_In$Layer, "DRP" = Preds_DRP_In$Layer, "FW" = Preds_FW_In$Layer, "NC" = Preds_NC_In$Layer, "NT" = Preds_NT_In$Layer, "NS" = Preds_NS_In$Layer, "R" = Preds_R_In$Layer, "SC" = Preds_SC_In$Layer)
qsave(Pred_In_ls, file = "output/predictions/Pred_In.qs", preset = "fast")

Khab_Woody_In <- ZStats_Khab_In
ZStats_WoodyRem_In <- ZStats_Woody_In
Khab_risk_In <- ZStats_Khab_In

for(kmr in KMRs){
  
  #### Filter data for properties with woody vegetation in 2011 and ensure clearing is < woody vegetation
  ZStats_WoodyRem_In[[kmr]] <- ZStats_WoodyRem_In[[kmr]] %>% 
    mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
    mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% 
    dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody, Khab = sum.Khab) %>%
    filter(N > 0)
  # NT <-  Response %>% dplyr::select(N) %>% as.matrix()
  
  # Warning message if woody vegetation is greater than Koala habitat
  if(nrow(ZStats_WoodyRem_In[[kmr]][ZStats_WoodyRem_In[[kmr]]$N < ZStats_WoodyRem_In[[kmr]]$Khab,])>0){message("KMR=", kmr, ": Woody vegetation greater than Koala habitat!!! \n\n")}
  
  ### Koala habitat risk data (Woody, Risk, Khab/Woody, Khab*Risk)
  Khab_risk_In[[kmr]] <- bind_cols(Pred_In_ls[[kmr]], ZStats_WoodyRem_In[[kmr]]) %>%
    mutate(Khab = ifelse(Khab < N, Khab, N),
           Khab_P = Khab/N, KhabRisk = PredAll*Khab_P) %>% 
    dplyr::select(Woody = N, Risk = PredAll, Khab_P, KhabRisk)
}

qsave(Khab_risk_In, file = "output/predictions/Khab_risk_In.qs", preset = "fast")
Khab_risk_In_all <- do.call(rbind, Khab_risk_In)
st_write(Khab_risk_In_all, "output/predictions/Khab_risk_In.shp", delete_layer = TRUE)


# # Plot observed vegetation loss ----
# ## Agriculture ----
# SUs_ZStats_Woody_Ag <- bind_cols(do.call(rbind, qread("output/spatial_units/SUs_Ag.qs")), 
#                                  do.call(rbind, qread("output/data/ZStats_Woody_Ag.qs")))
# SUs_ZStats_Woody_Ag <- SUs_ZStats_Woody_Ag %>% mutate(sum.aloss_lg = log(sum.aloss+1))
# WoodLoss_Ag <- ggplot(SUs_ZStats_Woody_Ag)+
#   geom_sf(aes(fill = sum.aloss_lg), lwd = 0)+
#   scale_fill_viridis_c(direction = -1, option = "magma")+
#   geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.5)+
#   theme_void()
# ggsave(filename = "output/predictions/WoodLoss_Ag.png", plot = WoodLoss_Ag, units = "cm", width = 15, height = 15, dpi = 300)
# 
# ## Forestry ----
# SUs_ZStats_Woody_Fo <- bind_cols(do.call(rbind, qread("output/spatial_units/sus_Fo.qs")), 
#                                  do.call(rbind, qread("output/data/ZStats_Woody_Fo.qs")))
# SUs_ZStats_Woody_Fo <- SUs_ZStats_Woody_Fo %>% mutate(sum.floss_lg = log(sum.floss+1))
# WoodLoss_Fo <- ggplot(SUs_ZStats_Woody_Fo)+
#   geom_sf(aes(fill = sum.floss_lg), lwd = 0)+
#   scale_fill_viridis_c(direction = -1)+
#   geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.5)+
#   theme_void()
# ggsave(filename = "output/predictions/WoodLoss_Fo.png", plot = WoodLoss_Fo, units = "cm", width = 15, height = 15, dpi = 300)
# 
# ## Infrastructure ----
# SUs_ZStats_Woody_In <- bind_cols(do.call(rbind, qread("output/spatial_units/sus_In.qs")), 
#                                  do.call(rbind, qread("output/data/ZStats_Woody_In.qs")))
# SUs_ZStats_Woody_In <- SUs_ZStats_Woody_In %>% mutate(sum.iloss_lg = log(sum.iloss+1))
# WoodLoss_In <- ggplot(SUs_ZStats_Woody_In)+
#   geom_sf(aes(fill = sum.iloss_lg), lwd = 0)+
#   scale_fill_viridis_c(direction = -1, option = "cividis")+
#   geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.5)+
#   theme_void()
# ggsave(filename = "output/predictions/WoodLoss_In.png", plot = WoodLoss_In, units = "cm", width = 15, height = 15, dpi = 300)


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
