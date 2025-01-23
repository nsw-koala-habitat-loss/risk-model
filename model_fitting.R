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

# file directories
OUTPUT_DIR <- "D:/Data/NSW_Deforestation/risk-model/output/"

# load Pre-processed data----
ZStats_Woody <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody.qs"))
ZStats_CovsC <- qread(file.path(OUTPUT_DIR, "data/ZStats_CovsC.qs"))
ZStats_CovsD <- qread(file.path(OUTPUT_DIR, "data/ZStats_CovsD.qs"))
SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/sus.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))

# # Load proposed covariates based on workshops from lookup xlsx
# CovLookup <- readxl::read_xlsx("Input/covariates/covariate_description.xlsx", sheet = "AllLyr")

# add area to the continuous covariates
# then rescale all cont. covariates to have mean of zero and SD of one 
# except for PCA covariates and Koala Habitat suitability
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
qsave(ZStats_Covs, file = file.path(OUTPUT_DIR, "data/ZStats_Covs.qs"), preset = "fast")


# CHECKING FOR MULTI-COLLINEARITY ---- 
## across all clearing types and all KMRs

## Continuous Covariates ----
ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC)
Corr_Cont <- cor(ZStats_CovsC_all, use = "complete.obs")
Corr_Cont_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+  ## highlight variables with correlations  > 0.5 OR < -0.5
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
# Export correlation plot
ggsave(Corr_Cont_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_all_plot.png"), width = 2000, height = 2000, units = "px")

## Discrete Covariates ----
ZStats_CovsD_all <- do.call(rbind, ZStats_CovsD)

# Make Cramers V matrix
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
# Export correlation plot
ggsave(Corr_Categ_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Categ_all_plot.png"), width = 2000, height = 2000, units = "px")

## Remove covariates that are collinear ----
### Remove variables with correlations  > 0.6 OR < -0.6
ZStats_Woody <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody.qs"))
ZStats_Covs <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs.qs"))
for (i in names(ZStats_Covs)) {
    ZStats_Covs[[i]] <- ZStats_Covs[[i]] %>% dplyr::select(-c(Elev, ForType))
}

# Select covariates for each clearing type----
## Agricultural Clearing ----

SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/sus.qs"))
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
Corr_Cont_Ag_all_plot
# ggsave(Corr_Cont_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

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
Corr_CovD_Ag_all_plot
# ggsave(Corr_CovD_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Ag, file = file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"), preset = "fast")
qsave(ZStats_Woody_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"), preset = "fast")
qsave(ZStats_Covs_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"), preset = "fast")
qsave(ZStats_Khab_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_Ag.qs"), preset = "fast")

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

SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/sus.qs"))
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
ggsave(Corr_Cont_Fo_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_Fo_all_plot.png"), width = 2000, height = 2000, units = "px")

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
ggsave(Corr_CovD_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Fo, file = file.path(OUTPUT_DIR, "spatial_units/sus_fo.qs"), preset = "fast")
qsave(ZStats_Woody_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"), preset = "fast")
qsave(ZStats_Covs_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"), preset = "fast")
qsave(ZStats_Khab_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_Fo.qs"), preset = "fast")

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
SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/sus.qs"))
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
ggsave(Corr_Cont_In_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_In_all_plot.png"), width = 2000, height = 2000, units = "px")

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
ggsave(Corr_CovD_In_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_In_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_In, file = file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"), preset = "fast")
qsave(ZStats_Woody_In, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"), preset = "fast")
qsave(ZStats_Covs_In, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"), preset = "fast")
qsave(ZStats_Khab_In, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_In.qs"), preset = "fast")

### Check for NAs in the data ----
map(ZStats_Covs_In, ~summary(.))
map(ZStats_Covs_In, ~nrow(.))

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

## Check DIC variations ----
#### If DIC varies significantly in repeated runs, then we need to set custom priors
### Agricultural clearing in each KMR ----

ZStats_Woody_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"))
ZStats_Covs_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"))
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
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
SUs_Fo <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_fo.qs"))
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))
ZStats_Covs_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
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
SUs_In <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"))
ZStats_Woody_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"))
ZStats_Covs_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
KMRs <- names(SUs_In)
inla.setOption(num.threads="8:1")
# map(ZStats_Covs_In, ~summary(.))

ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:5){
    Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "Area + DistCity + DistRoad + Drought + EcolCond + LandTen + LandUse + PlanZone + PopDen + PopGro + Precip + PropVal + ScEc_PC1 + ScEc_PC3 + ScEc_PC4 + ScEc_PC5 + Temp", Verbose = FALSE, N_retry=5, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_In$PModel$dic$dic + Model_In$NModel$dic$dic, 2)))
  }
  toc(log = TRUE)
  cat("\n\n")
}
proc.time() - ptm


# ptm <- proc.time()
# Test_2 <- fit_model2(KMR="CC", ClearType = 1, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, OutputDir = file.path(OUTPUT_DIR, "models/")
# proc.time() - ptm
# Test_2 <- qread(file.path(OUTPUT_DIR, "models/Model_CC_Ag_.qs"))
# summary(Test_2$PModel)
# 
# formula_1 <- as.formula(paste0(paste("P", paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
# 
# 
# summary(Test_2$PModel)


# MODEL SELECTIONS----
### Each clearing types and selection directions was set to run independently and can be run in parallel

## Agricultural----
## Load data
ZStats_Woody_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"))
ZStats_Covs_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"))
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
KMRs <- names(SUs_Ag)

### Adjust inla.setOption(num.threads= "8:1") to allocate the number of threads to be used in 1 inla run. The number need to be reduce if running parallel.
# inla.setOption(num.threads="8:1")

## Forward model selection
tic("Select_model: Ag forward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, 
               RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, 
               SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)

## Forward model selection
tic("Select_model: Ag backward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, 
               RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, 
               SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)

## Forestry----
## Load data
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))
ZStats_Covs_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"))
SUs_Fo <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_Fo.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
KMRs <- names(SUs_Fo)

### Adjust inla.setOption(num.threads= "8:1") to allocate the number of threads to be used in 1 inla run. The number need to be reduce if running parallel.
# inla.setOption(num.threads="8:1")

## Forward model selection
tic("Select_model: Fo forward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, 
               RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, 
               SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)

## Backward model selection
tic("Select_model: Fo backward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, 
               RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, 
               SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)

## Infrastructure----
## Load data
ZStats_Woody_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"))
ZStats_Covs_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"))
SUs_In <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
KMRs <- names(SUs_In)

# inla.setOption(num.threads="8:1")

## Forward model selection
toc("Select_model: In forward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, 
               RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, 
               SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)


## Backward model selection
tic("Select_model: In backward selection")
for(kmr in KMRs){
  Select_model(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, 
               RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, 
               SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, 
               N_retry=3, Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
}
toc(log = TRUE)


## Check model selection step ---- 
### Check errors in Model selections steps ---
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
KMRs <- names(SUs_Ag)
rm(SUs_Ag)
ClrTyps <- c("Ag", "Fo", "In")

for(ClrTyp in ClrTyps){
  for(kmr in KMRs){
    SelModel_FC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_" , ClrTyp,"_FC.qs")))
    SelModel_BC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_" , ClrTyp,"_FC.qs")))
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

SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
KMRs <- names(SUs_Ag)
rm(SUs_Ag)
ClrTyps <- c("Ag", "Fo", "In")

for(ClrTyp in ClrTyps){  ## Agriculture, Forestry, Infrastructure
  for(kmr in KMRs){  ## KMRs
    
    # read model selection results for both forward and backward selection
    SelModel_BC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_" , ClrTyp , "_BC.qs")))
    SelModel_FC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_" , ClrTyp , "_FC.qs")))
    
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
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
ZStats_Woody_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"))
ZStats_Covs_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"))
KMRs <- names(ZStats_Covs_Ag)
kmr <- KMRs[1]
Model_Ag <- fit_model2(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_Ag <- summary(Model_Ag$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="CC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_Ag_BC.qs")))
  SelModel_FC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_Ag_FC.qs")))
  
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
write.csv(Cov_ls_Ag_tab, file.path(OUTPUT_DIR, "models/Cov_ls_Ag.csv"), row.names = FALSE, na = "")

### Forestry ----
SUs_Fo <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Fo.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))
ZStats_Covs_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"))
KMRs <- names(ZStats_Covs_Fo)
kmr <- KMRs[1]
Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_Fo <- summary(Model_Fo$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="CC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_Fo_BC.qs")))
  SelModel_FC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_Fo_FC.qs")))
  
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
write.csv(Cov_ls_tab, file.path(OUTPUT_DIR, "models/Cov_ls_Fo.csv"), row.names = FALSE, na = "")

### Infrastructure ----
SUs_In <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_In.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
ZStats_Woody_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"))
ZStats_Covs_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"))
KMRs <- names(ZStats_Covs_In)
kmr <- KMRs[1]
Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_In <- summary(Model_In$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)
kmr="NC"

for (kmr in KMRs){
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_In_BC.qs")))
  SelModel_FC <- qread(file.path(OUTPUT_DIR, paste0("models/SelModel_", kmr, "_In_FC.qs")))
  
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

write.csv(Cov_ls_tab, file.path(OUTPUT_DIR, "models/Cov_ls_In.csv"), row.names = FALSE, na = "")

# MODEL PREDICTIONS ----
## Refit the best model ----
### Identify model with lower DIC (Between FC & BC)
### Refit the model to save internal GMRF approximations for inla.posterior.sample(): control.compute = (list(config = TRUE))
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_Ag.qs"))
KMRs <- names(SUs_Ag)
ClearTypes <- c(1, 2, 3)

walk(ClearTypes, function(ClearTypes_val){
  walk(KMRs, function(kmr_val){
    print(paste0("Refit model for ", kmr_val, " with ClearType ", ClearTypes_val))
    refit_model(KMR = kmr_val, ClearType = ClearTypes_val, ModelDir = file.path(OUTPUT_DIR, "models/"))
  })
})

refit_model(KMR = "CC", ClearType = 1, ModelDir = file.path(OUTPUT_DIR, "models/"))

# Prediction by sampling ----

### Posterior sample generated by high-performance computing platform (‘Bunya’ AMD epyc3 Milan compute cores with memory 28GB to 1.5TB) administrated by The University of Queensland Research Computing Centre (2024).
### Reference: The University of Queensland Research Computing Centre. 2024. Bunya supercomputer. Brisbane, Queensland, Australia. https://dx.doi.org/10.48610/wf6c-qy55
### Note:See /SLURM/ folder and /R/ folder for the code to run the prediction on HPC

## Example for the prediction by sampling
# Load directory
MODEL <- qread(file.path(OUTPUT_DIR, "models/Model_CC_Ag.qs"))
CT <- if(MODEL$ClearType == "1"){"Ag"} else if(MODEL$ClearType == "2"){"In"} else if(MODEL$ClearType == "3"){"Fo"}
N <- 50000 # Number of samples
cat("\n\nRun prediction by sampling " , N , "times for\nKMR: ",  MODEL$KMR, "\nClear Type: ", CT, "\n\n")

# Predictions
Pred <- predict_model3(model = MODEL, N = N, RandEff = "SA1ID")

# Save Predictions
output_name <- file.path(file.path(OUTPUT_DIR, "predictions/Pred_CC_Ag.qs"))
cat("\n\nSave Predictions to: ", output_name, "\n")
qsave(Pred, output_name)

# Generate shapefile output----
walk(1:3, ~Combine_Predictions(ClearType = .x, Prediction_DIR = file.path(OUTPUT_DIR, "predictions/"), WRITE_SHP = TRUE, WRITE_DATA = TRUE))

## Save predictions to File Geodatabase (GDB)
## If the Geodatabase export is done in Combine_Predictions function then this part can be skipped
# gdb_path <- file.path(OUTPUT_DIR, "predictions.gdb"
# Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
# st_write(Pred_Ag, dsn = gdb_path, layer = "Pred_Ag", driver = "OpenFileGDB")
# Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))
# st_write(Pred_In, dsn = gdb_path, layer = "Pred_In", driver = "OpenFileGDB")
# Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))
# st_write(Pred_Fo, dsn = gdb_path, layer = "Pred_Fo", driver = "OpenFileGDB")

# Koala habitat loss risk ----

## Load predictions output data and Woody vegetation data (containing Koala habitat loss column)
Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))
Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))

ZStats_Woody_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"))
ZStats_Woody_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"))
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))

## Lookup dataframe for KMR acronyms and KMR Full names
KMR_DF <- Pred_Ag %>% st_drop_geometry() %>%  dplyr::select(KMR) %>% unique() %>% mutate(KMR_a = str_extract_all(KMR, "\\b[A-Za-z]") %>% sapply(paste, collapse = ""))

Khab_data_Ag <- Prep_Khab(ZStats_Woody = ZStats_Woody_Ag, KMR_name_df = KMR_DF)
Khab_data_In <- Prep_Khab(ZStats_Woody = ZStats_Woody_In, KMR_name_df = KMR_DF)
Khab_data_Fo <- Prep_Khab(ZStats_Woody = ZStats_Woody_Fo, KMR_name_df = KMR_DF)

## Calculate Koala habitat loss risk
Khab_risk_Ag <- Get_Khab_loss_risk(Pred_data = Pred_Ag, Khab_data = Khab_data_Ag)
Khab_risk_In <- Get_Khab_loss_risk(Pred_data = Pred_In, Khab_data = Khab_data_In)
Khab_risk_Fo <- Get_Khab_loss_risk(Pred_data = Pred_Fo, Khab_data = Khab_data_Fo)

## Save the output
qsave(Khab_risk_Ag, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"), preset = "fast")
qsave(Khab_risk_In, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"), preset = "fast")
qsave(Khab_risk_Fo, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"), preset = "fast")

Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))
Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))

st_write(obj = Khab_risk_Ag, dsn = gdb_path, layer = "Khab_risk_Ag", driver = "OpenFileGDB")
st_write(obj = Khab_risk_In, dsn = gdb_path, layer = "Khab_risk_In", driver = "OpenFileGDB")
st_write(obj = Khab_risk_Fo, dsn = gdb_path, layer = "Khab_risk_Fo", driver = "OpenFileGDB")
st_layers(gdb_path)
st_write(Khab_risk_Ag, file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.shp"), delete_layer = TRUE, append = FALSE)
st_write(Khab_risk_In, file.path(OUTPUT_DIR, "predictions/Khab_risk_In.shp"), delete_layer = TRUE, append = FALSE)
st_write(Khab_risk_Fo, file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.shp"), delete_layer = TRUE, append = FALSE)
