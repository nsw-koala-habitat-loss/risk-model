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
options(max.print = 999)

# load functions
source("functions.R")

## Note: change the path to the input and output directory

INPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/input")
OUTPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/output")

# load Pre-processed data----
ZStats_Woody <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody.qs"))
ZStats_CovsC <- qread(file.path(OUTPUT_DIR, "data/ZStats_CovsC.qs"))
ZStats_CovsD <- qread(file.path(OUTPUT_DIR, "data/ZStats_CovsD.qs"))
SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/sus.qs"))
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))

## Check all data has same number of rows
walk(names(SUs), ~{
  if (nrow(SUs[[.x]]) != nrow(ZStats_Woody[[.x]])) {
    stop(paste0("Number of rows in SUs and ZStats_Woody do not match for ", .x))
  } else {
    print(paste0("Number of rows in SUs and ZStats_Woody match for ", .x, ": ", nrow(SUs[[.x]])))
  }
  if (nrow(SUs[[.x]]) != nrow(ZStats_CovsC[[.x]])) {
    stop(paste0("Number of rows in SUs and ZStats_CovsC do not match for ", .x))
  } else {
    print(paste0("Number of rows in SUs and ZStats_CovsC match for ", .x, ": ", nrow(SUs[[.x]])))
  }
  if (nrow(SUs[[.x]]) != nrow(ZStats_CovsD[[.x]])) {
    stop(paste0("Number of rows in SUs and ZStats_CovsD do not match for ", .x))
  } else {
    print(paste0("Number of rows in SUs and ZStats_CovsD match for ", .x, ": ", nrow(SUs[[.x]])))
  }
})

# Rescale all cont. covariates 
## mean of zero and SD of one except for PCA covariates and Koala Habitat suitability
for (i in names(ZStats_CovsC)) {
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]]
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
ZStats_Woody_Ag <- ZStats_Woody #%>% map(function(x)nrow(x))
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

    # Only include SUs in Freehold and leasehold land tenure type (removing 25658 SUs)
    filter(LandTen %in% c("1", "2")) %>%
    
    # Remove NAs
    tidyr::drop_na(SA1) %>% tidyr::drop_na(.) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse),
                  NatVegReg = droplevels(NatVegReg),
                  LandTen = droplevels(LandTen))
  
  ## Split SUs, Covs, Woody and KoalaHabitat into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% dplyr::select(names(SUs_Ag[[i]]))
  print(paste0("Number of SUs in ", i, ": ", nrow(SUs_Ag[[i]])))
  ## Select covariates for Agricultural Clearing 
  ZStats_Covs_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Ag[[i]]))%>% dplyr::select(-LandTen)
  print(paste0("Number of covariates in ", i, ": ", nrow(ZStats_Covs_Ag[[i]])))
  ## Select woody extent and loss columns
  ZStats_Woody_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Ag[[i]]))
  print(paste0("Number of woody extent and loss columns in ", i, ": ", nrow(ZStats_Woody_Ag[[i]])))
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_Ag[[i]] <- SUs_ZStats_Ag[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select("sum.woody","sum.aloss", "sum.Khab" )
  print(paste0("Number of Koala habitat suitability columns in ", i, ": ", nrow(ZStats_Khab_Ag[[i]])))
}

### Check for collinearity in Agricultural Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates 
ZStats_CovsC_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf, Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)
nrow(ZStats_CovsC_Ag_all)
#### Calculate correlation matrix
Corr_Cont_Ag_all <- cor(ZStats_CovsC_Ag_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_Ag_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Ag_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_Cont_Ag_all_plot
ggsave(Corr_Cont_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Discrete Covariates ----
#### Combine all Covariates for all KMRs then select discrete covariates
ZStats_CovsD_Ag_all <- do.call(rbind, ZStats_Covs_Ag) %>% 
  dplyr::select(NatVegReg, LandUse, Drought, Fire)

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
ggsave(Corr_CovD_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Ag, file = file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"), preset = "fast")
nrow(SUs_Ag %>% bind_rows())
qsave(ZStats_Woody_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"), preset = "fast")
nrow(ZStats_Woody_Ag %>% bind_rows())
qsave(ZStats_Covs_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"), preset = "fast")
nrow(ZStats_Covs_Ag %>% bind_rows())
qsave(ZStats_Khab_Ag, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_Ag.qs"), preset = "fast")
nrow(ZStats_Khab_Ag %>% bind_rows())

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
    # Remove SUs in LLS excluded area (non-rural) 
    # Remove Multiple-use public forest(2) and Nature conservation reserve (3) Forest Tenure Type
    filter(LandUse != "5", NatVegReg != "0", TenType == "1") %>% 
    
    # Remove NAs
    tidyr::drop_na(.) %>% tidyr::drop_na(SA1) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse),
                  NatVegReg = droplevels(NatVegReg),
                  TenType = droplevels(TenType))
  
  ## Split SUs, Covs, Woody and Koala Habitat into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% dplyr::select(names(SUs_Fo[[i]]))
  print(paste0("Number of SUs in ", i, ": ", nrow(SUs_Fo[[i]])))
  ## Select covariates for Forestry Clearing
  ZStats_Covs_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_Fo[[i]])) %>% dplyr::select(-TenType)
  print(paste0("Number of covariates in ", i, ": ", nrow(ZStats_Covs_Fo[[i]])))
  ## Select woody extent and loss columns
  ZStats_Woody_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_Fo[[i]]))
  print(paste0("Number of woody extent and loss columns in ", i, ": ", nrow(ZStats_Woody_Fo[[i]])))
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_Fo[[i]] <- SUs_ZStats_Fo[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Khab_Fo[[i]]))
  print(paste0("Number of Koala habitat suitability columns in ", i, ": ", nrow(ZStats_Khab_Fo[[i]])))
}

### Check for collinearity in Forestry Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates 
ZStats_CovsC_Fo_all <- do.call(rbind, ZStats_Covs_Fo) %>% 
  dplyr::select(PopDen, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, AgProf,Soil_PC1, Soil_PC2, Soil_PC3, slope, Precip, Temp, EcolCond, Area)
nrow(ZStats_CovsC_Fo_all)
#### Calculate correlation matrix
Corr_Cont_Fo_all <- cor(ZStats_CovsC_Fo_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_Fo_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_Fo_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = FALSE, alpha = FALSE)
Corr_Cont_Fo_all_plot
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
Corr_CovD_Ag_all_plot
ggsave(Corr_CovD_Ag_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_Ag_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_Fo, file = file.path(OUTPUT_DIR, "spatial_units/sus_fo.qs"), preset = "fast")
nrow(SUs_Fo %>% bind_rows())
qsave(ZStats_Woody_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"), preset = "fast")
nrow(ZStats_Woody_Fo %>% bind_rows())
qsave(ZStats_Covs_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"), preset = "fast")
nrow(ZStats_Covs_Fo %>% bind_rows())
qsave(ZStats_Khab_Fo, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_Fo.qs"), preset = "fast")
nrow(ZStats_Khab_Fo %>% bind_rows())


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
    
    # Keep only SUs in Freehold and leasehold land tenure type (removing 133215 SUs)
    filter(LandTen %in% c("1", "2")) %>%

    # Remove NAs
    tidyr::drop_na(SA1) %>% tidyr::drop_na(.) %>% 
    
    # Remove drop unused levels
    dplyr::mutate(LandUse = droplevels(LandUse))
    
  
  ## Split SUs, Covs and Woody into individual dataframes
  ## Select columns corresponding to Spatial Units
  SUs_In[[i]] <- SUs_ZStats_In[[i]] %>% dplyr::select(names(SUs_In[[i]]))
  print(paste0("Number of SUs in ", i, ": ", nrow(SUs_In[[i]])))
  ## Select covariates for Infrastructure Clearing
  ZStats_Covs_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Covs_In[[i]]))%>% dplyr::select(-Remoteness, -LandTen)
  print(paste0("Number of covariates in ", i, ": ", nrow(ZStats_Covs_In[[i]])))
  ## Select woody extent and loss columns
  ZStats_Woody_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Woody_In[[i]])) 
  ## Select Woody extent and loss columns and koala habitat suitability columns
  ZStats_Khab_In[[i]] <- SUs_ZStats_In[[i]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(names(ZStats_Khab_In[[i]]))
}

### Check for collinearity in Infrastructure Clearing----
#### Continuous Covariates ----
#### Combine all Covariates for all KMRs then select continuous covariates
ZStats_CovsC_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(PopDen, PopGro, ScEc_PC1, ScEc_PC2, ScEc_PC3, ScEc_PC4, ScEc_PC5, DistRoad, DistCity, PropVal, slope, Precip, Temp, EcolCond, Area)
nrow(ZStats_CovsC_In_all)
#### Calculate correlation matrix
Corr_Cont_In_all <- cor(ZStats_CovsC_In_all, use = "complete.obs")

#### Plot correlation matrix
Corr_Cont_In_all_plot <- ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont_In_all, label = TRUE, hjust = 1, layout.exp = 2)+ 
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+ 
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) + 
  guides(color = "none", alpha = "none")
Corr_Cont_In_all_plot
ggsave(Corr_Cont_In_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_Cont_In_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Discrete Covariates ----
#### Combine all Covariates for all KMRs then select discrete covariates
ZStats_CovsD_In_all <- do.call(rbind, ZStats_Covs_In) %>% 
  dplyr::select(PlanZone, LandUse, Drought)
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
Corr_CovD_In_all_plot
ggsave(Corr_CovD_In_all_plot, file = file.path(OUTPUT_DIR, "collinearity/Corr_CovD_In_all_plot.png"), width = 2000, height = 2000, units = "px")

#### Export data
qsave(SUs_In, file = file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"), preset = "fast")
nrow(SUs_In %>% bind_rows())
qsave(ZStats_Woody_In, file = file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"), preset = "fast")
nrow(ZStats_Woody_In %>% bind_rows())
qsave(ZStats_Covs_In, file = file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"), preset = "fast")
nrow(ZStats_Covs_In %>% bind_rows())
qsave(ZStats_Khab_In, file = file.path(OUTPUT_DIR, "data/ZStats_Khab_In.qs"), preset = "fast")
nrow(ZStats_Khab_In %>% bind_rows())

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
nrow(ZStats_Woody_Ag %>% bind_rows())
ZStats_Covs_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Ag.qs"))
nrow(ZStats_Covs_Ag %>% bind_rows())
SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs_Ag.qs"))
nrow(SUs_Ag %>% bind_rows())
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
# inla.setOption(num.threads="8:1")
# map(ZStats_Covs_In, ~summary(.))

ptm <- proc.time()
for(kmr in KMRs){
  tic(paste0(kmr))
  cat(paste0("Running full model for: ", kmr, "\n"))
  for(i in 1:5){
    Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Verbose = FALSE, N_rerun = 1, N_retry=5, Initial_Tlimit = 1000, OutputDir = NULL)
    print(paste0(kmr, ":  DIC= " , round(Model_In$PModel$dic$dic + Model_In$NModel$dic$dic, 2)))
  }
  toc(log = TRUE)
  cat("\n\n")
}
proc.time() - ptm


## Try informative priors for Infrastructure clearing in SC KMR ----
### Model convergence issues in SC KMR for Infrastructure clearing
### Try informative priors based on preliminary model runs
### Then use higher precision (lower INLA tolerance) for model fitting
kmr <- "SC"

for(i in 1:20){
  Model_In <- fit_model3(KMR = kmr , ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, 
                         Verbose = FALSE, N_rerun=2 , N_retry=4, NMod_TOLN = 1e-14, Initial_Tlimit = 1000, OutputDir = NULL)
  cat(paste0(kmr, ":  DIC= " , round(Model_In$PModel$dic$dic + Model_In$NModel$dic$dic, 2), "\n"))
  cat(paste0("PModel dMLIK = ", abs(Model_In$PModel$mlik[1,1] - Model_In$PModel$mlik[2,1]), "\n"))
  cat(paste0("NModel dMLIK = ", abs(Model_In$NModel$mlik[1,1] - Model_In$NModel$mlik[2,1]), "\n"))
}

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
nrow(ZStats_Woody_Fo %>% bind_rows())
ZStats_Covs_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"))
nrow(ZStats_Covs_Fo %>% bind_rows())
SUs_Fo <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_Fo.qs"))
nrow(SUs_Fo %>% bind_rows())
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
nrow(SA1s %>% bind_rows())
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
nrow(ZStats_Woody_In %>% bind_rows())
ZStats_Covs_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_In.qs"))
nrow(ZStats_Covs_In %>% bind_rows())
SUs_In <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"))
nrow(SUs_In %>% bind_rows())
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
nrow(SA1s %>% bind_rows())
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


### Rerun model selection for KMR==SC because of the DIC convergence issue ----
### Use select_model2 function that use informative prior and use higher precision (lower INLA tolerance) for model fitting
kmr <- "SC"
inla.setOption(num.threads="16:1")

#### Forward model selection
toc("Select_model: In forward selection")
Select_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, 
               RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, 
               SA1sPoly = SA1s, Direction = "FC", Verbose = FALSE, 
               NMod_TOLN = 1e-13, N_rerun = 3, N_retry=3, 
               Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
toc(log = TRUE)

#### Backward model selection
tic("Select_model: In backward selection")
Select_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, 
             RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, 
             SA1sPoly = SA1s, Direction = "BC", Verbose = FALSE, 
             NMod_TOLN = 1e-13, N_rerun = 3, N_retry=3, 
             Initial_Tlimit = 1000, OutputDir = file.path(OUTPUT_DIR, "models/"))
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
    
    # Get the difference in mlik integration and Gaussian for both models
    SelModel_BC_PModel_dMLIK <- abs(SelModel_BC$PModel$mlik[1,1] - SelModel_BC$PModel$mlik[2,1])
    SelModel_BC_NModel_dMLIK <- abs(SelModel_BC$NModel$mlik[1,1] - SelModel_BC$NModel$mlik[2,1])
    SelModel_FC_PModel_dMLIK <- abs(SelModel_FC$PModel$mlik[1,1] - SelModel_FC$PModel$mlik[2,1])
    SelModel_FC_NModel_dMLIK <- abs(SelModel_FC$NModel$mlik[1,1] - SelModel_FC$NModel$mlik[2,1])

    # Print results
    cat(ClrTyp, kmr, ":\n")
    
    # Backward complete model selection
    cat("BC:\n")
    cat("DIC_ls: Min DIC:", unlist(MinDIC_BC), "\n")
    cat(MinDIC_BC_expl, "\n")
    cat("Sel_Model: DIC:", SelModelDIC_BC, "\n" )
    cat(SelModel_BC_Covs, "\n\n")
    cat("PModel dMLIK:", SelModel_BC_PModel_dMLIK, "\n")
    cat("NModel dMLIK:", SelModel_BC_NModel_dMLIK, "\n")

    if(MinDIC_BC_expl != SelModel_BC_Covs){message("Selected model does not match the minimum DIC model!!! \n")}
    if(SelModel_BC_PModel_dMLIK > 2 | SelModel_BC_NModel_dMLIK > 2){message("Selected model has dMLIK values greater than 2!!! \n\n")}
    
    # Forward complete model selection
    cat("FC:\n")
    cat("DIC_ls: Min DIC:", unlist(MinDIC_FC), "\n")
    cat(MinDIC_FC_expl, "\n")
    cat("Sel_Model: DIC:", SelModelDIC_FC, "\n")
    cat(SelModel_FC_Covs, "\n\n")
    cat("PModel dMLIK:", SelModel_FC_PModel_dMLIK, "\n")
    cat("NModel dMLIK:", SelModel_FC_NModel_dMLIK, "\n")
    if(MinDIC_FC_expl != SelModel_FC_Covs){message("Selected model does not match the minimum DIC model!!! \n\n")}
    if(SelModel_FC_PModel_dMLIK > 2 | SelModel_FC_NModel_dMLIK > 2){message("Selected model has dMLIK values greater than 2!!! \n\n")}
    
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
                      "Slope", "Rainfall", "Temperature", "Ecological condition", "Property size",
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
nrow(SUs_Fo %>% bind_rows())
SA1s <- qread(file.path(OUTPUT_DIR, "spatial_units/sa1s.qs"))
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))
ZStats_Covs_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs_Fo.qs"))
KMRs <- names(ZStats_Covs_Fo)
kmr <- KMRs[1]
Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_Fo <- summary(Model_Fo$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)

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
kmr <- KMRs[9]
Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, N_rerun = 0, NMod_TOLN = 1e-11, Initial_Tlimit = 1000, OutputDir = NULL)
Cov_ls_In <- summary(Model_In$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate)

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
                      "Distance to road", "Distance to urban centre", "Property value", "Slope", "Rainfall", "Temperature", "Ecological condition", "Property size",
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
ClrTyps <- c("Ag", "In", "Fo")
ClearTypes <- c(1, 2, 3)

walk(ClearTypes, function(ClearTypes_val){
  walk(KMRs, function(kmr_val){
    print(paste0("Refit model for ", kmr_val, " with ClearType ", ClearTypes_val))
    refit_model(KMR = kmr_val, ClearType = ClearTypes_val, N_rerun = if(ClearTypes_val==2){6}else{2}, ModelDir = file.path(OUTPUT_DIR, "models/"))
  })
})

## Refit KMR SC with ClearType In again with refit_model and lower tolerance level to force convergence ----
kmr <- "SC"
ClearType <- 2
ClrTyp <- "In"

refit_model(KMR = kmr, ClearType = ClearType, N_rerun = 4, NMod_TOLN = 1e-13, ModelDir = file.path(OUTPUT_DIR, "models/"))

## Check MLIK integration and Gaussian approximation for the refitted model
walk(ClrTyps, function(ClrTyp){
  walk(KMRs, function(kmr){
    MODEL <- qread(file.path(OUTPUT_DIR, paste0("models/Model_", kmr, "_", ClrTyp, ".qs")))
    dMLIK_PModel <- abs(MODEL$PModel$mlik[1,1] - MODEL$PModel$mlik[2,1])
    dMLIK_NModel <- abs(MODEL$NModel$mlik[1,1] - MODEL$NModel$mlik[2,1])
    cat("Model: ", kmr, " ClearType: ", ClrTyp, "\n")
    cat("DIC: ", MODEL$PModel$dic$dic + MODEL$NModel$dic$dic, "\n")
    cat("PModel dMLIK: ", dMLIK_PModel, "\n")
    cat("NModel dMLIK: ", dMLIK_NModel, "\n")
  })
})

# Export diagnostic numbers and plots for the refitted model: -----
## Export csv dMLIK values----
dMLIK_ClrTyps <- map(ClrTyps, function(ClrTyp){
   dMLIK_KMR<- map(KMRs, function(kmr_val){
    MODEL <- qread(file.path(OUTPUT_DIR, paste0("models/Model_", kmr_val, "_", ClrTyp, ".qs")))
    dMLIK_PModel <- abs(MODEL$PModel$mlik[1,1] - MODEL$PModel$mlik[2,1])
    dMLIK_NModel <- abs(MODEL$NModel$mlik[1,1] - MODEL$NModel$mlik[2,1])
    return(data.frame(KMR = kmr_val, ClearType = ClrTyp, dMLIK_PModel = dMLIK_PModel, dMLIK_NModel = dMLIK_NModel))
  })
})
dMLIK_ClrTyps_BD <- bind_rows(dMLIK_ClrTyps)
row.names(dMLIK_ClrTyps_BD) <- NULL
write.csv(dMLIK_ClrTyps_BD, file = file.path(OUTPUT_DIR, "models/dMLIK.csv"), row.names = FALSE, na = "")


# Prediction by sampling ----

### Posterior sample generated by high-performance computing platform (‘Bunya’ AMD epyc3 Milan compute cores with memory 28GB to 1.5TB) administrated by The University of Queensland Research Computing Centre (2024).
### Reference: The University of Queensland Research Computing Centre. 2024. Bunya supercomputer. Brisbane, Queensland, Australia. https://dx.doi.org/10.48610/wf6c-qy55
### See "Generate_R_script.R" "Generate_SLURM_Script.R" for generating customised the script to run the prediction on HPC

## Example prediction by sampling from the posterior distribution ----
# Load directory
MODEL <- qread(file.path(OUTPUT_DIR, "models/Model_CC_Ag.qs"))

CT <- if(MODEL$ClearType == "1"){"Ag"} else if(MODEL$ClearType == "2"){"In"} else if(MODEL$ClearType == "3"){"Fo"}
N <- 20 # Number of samples
cat("\n\nRun prediction by sampling " , N , "times for\nKMR: ",  MODEL$KMR, "\nClear Type: ", CT, "\n\n")

# Predictions
Pred <- predict_model3(model = MODEL, N = N, RandEff = "SA1ID")
qread(file.path(OUTPUT_DIR, "spatial_units/sus_Ag.qs")) %>% map(nrow)
qread(file.path(OUTPUT_DIR, "spatial_units/sus_In.qs")) %>% map(nrow)
qread(file.path(OUTPUT_DIR, "spatial_units/sus_Fo.qs")) %>% map(nrow)

# Save Predictions
output_name <- file.path(file.path(OUTPUT_DIR, "predictions/Pred_CC_Ag.qs"))
cat("\n\nSave Predictions to: ", output_name, "\n")
qsave(Pred, output_name)

# Generate vector output----
## Check if geodatabase already exists, if not create one
# gdb_path <- file.path(OUTPUT_DIR, "predictions/Predictions.gdb")
# if(!file.exists(gdb_path)){
#   dir.create(gdb_path)
# }
walk(1:3, ~Combine_Predictions(ClearType = .x, Prediction_DIR = file.path(OUTPUT_DIR, "predictions/"), WRITE_SHP = TRUE, WRITE_DATA = TRUE))

## Combine_Predictions(ClearType = 2, Prediction_DIR = file.path(OUTPUT_DIR, "predictions/"), WRITE_SHP = TRUE, WRITE_DATA = TRUE)

# Extract Beta distribution precision (Phi) for each KMR ----
KMRs <- c("CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC")
NModel_Phi <- map(KMRs, ~qread(
  file.path(OUTPUT_DIR, "predictions", paste0("Pred_", .x, "_In.qs")))$NModel_Phi
  )
NModel_Phi <- unlist(NModel_Phi) 
names(NModel_Phi) <- KMRs

qsave(NModel_Phi, file = file.path(OUTPUT_DIR, "predictions/NModel_Phi_In.qs"))

# Koala habitat loss risk ----
gdb_path <- file.path(OUTPUT_DIR, "predictions/Predictions.gdb")

## Load predictions output data and Woody vegetation data (containing Koala habitat loss column)
Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))
Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))

SUs_Ag <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_Ag.qs"))
SUs_In <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_In.qs"))
SUs_Fo <- qread(file.path(OUTPUT_DIR, "spatial_units/sus_Fo.qs"))

ZStats_Woody_Ag <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Ag.qs"))
ZStats_Woody_In <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_In.qs"))
ZStats_Woody_Fo <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody_Fo.qs"))

## Add SUID and KMR to ZStats_Woody data
for(i in names(ZStats_Woody_Ag)){
  ZStats_Woody_Ag[[i]]$SUID <- SUs_Ag[[i]]$SUID
  ZStats_Woody_In[[i]]$SUID <- SUs_In[[i]]$SUID
  ZStats_Woody_Fo[[i]]$SUID <- SUs_Fo[[i]]$SUID
  ZStats_Woody_Ag[[i]]$KMR <- SUs_Ag[[i]]$KMR
  ZStats_Woody_In[[i]]$KMR <- SUs_In[[i]]$KMR
  ZStats_Woody_Fo[[i]]$KMR <- SUs_Fo[[i]]$KMR
}

## Calculate Koala habitat loss risk
Khab_risk_Ag <- Get_Khab_loss_risk(Pred_data = Pred_Ag, Khab_data = ZStats_Woody_Ag %>% do.call(rbind,.))
Khab_risk_In <- Get_Khab_loss_risk(Pred_data = Pred_In, Khab_data = ZStats_Woody_In %>% do.call(rbind,.))
Khab_risk_Fo <- Get_Khab_loss_risk(Pred_data = Pred_Fo, Khab_data = ZStats_Woody_Fo %>% do.call(rbind,.))


## Save the output
### Save to a geodatabase
st_layers(gdb_path)

st_write(obj = Khab_risk_Ag, dsn = gdb_path, layer = "Khab_risk_Ag", driver = "OpenFileGDB", append = FALSE)
st_write(obj = Khab_risk_In, dsn = gdb_path, layer = "Khab_risk_In", driver = "OpenFileGDB", append = FALSE)
st_write(obj = Khab_risk_Fo, dsn = gdb_path, layer = "Khab_risk_Fo", driver = "OpenFileGDB", append = FALSE)
st_layers(gdb_path)

qsave(Khab_risk_Ag, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"), preset = "fast")
qsave(Khab_risk_In, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"), preset = "fast")
qsave(Khab_risk_Fo, file = file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"), preset = "fast")

# Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
# Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))
# Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))

# Combining all Predictions and Koala habitat loss risk ----
## Load prediction data
Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))
Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))

## Load Koala habitat loss risk data
Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))
Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))

## Load spatial units data
SUs <- qread(file.path(OUTPUT_DIR, "spatial_units/SUs.qs")) %>% bind_rows()

## Load woody vegetation data
ZStats_Woody <- qread(file.path(OUTPUT_DIR, "data/ZStats_Woody.qs")) %>% bind_rows()

## Load covariate data
ZStats_Covs <- qread(file.path(OUTPUT_DIR, "data/ZStats_Covs.qs")) %>% bind_rows()

## Combine all data
SUs_ZStats <- bind_cols(SUs %>% dplyr::select(-Area), 
                        ZStats_Woody, 
                        ZStats_Covs) %>% 
              mutate(WoodClr = sum.aloss+sum.iloss+sum.floss,
                     Woody = ifelse(sum.woody < WoodClr, WoodClr, sum.woody)) %>%
              filter(LandTen %in% c("1", "2")) %>%
              dplyr::select(KMR, SA1, SUID, Woody, WoodClr, Area) 
## Create a vector of KMRs
KMRs_Long <- unique(SUs_ZStats$KMR)

# Create a named vector of KMRs with their full names
KMRs_Short <- sapply(KMRs_Long, function(name) {
  # Remove non-alphabetic characters
  clean_name <- gsub("[^a-zA-Z ]", "", name)
  # Extract the first letter of each word and concatenate them
  abbreviation <- paste0(toupper(substr(unlist(strsplit(clean_name, " ")), 1, 1)), collapse = "")
  return(abbreviation)
})

## Prepare data frames for joining
Khab_risk_Ag_DF <- st_drop_geometry(Khab_risk_Ag) %>% dplyr::select(KMR, SUID, KhabRisk_Ag = KhabRisk)
Khab_risk_In_DF <- st_drop_geometry(Khab_risk_In) %>% dplyr::select(KMR, SUID, KhabRisk_In = KhabRisk)
Khab_risk_Fo_DF <- st_drop_geometry(Khab_risk_Fo) %>% dplyr::select(KMR, SUID, KhabRisk_Fo = KhabRisk)

Pred_Ag_DF <- st_drop_geometry(Pred_Ag) %>% left_join(Khab_risk_Ag_DF, by = c("KMR", "SUID")) %>% dplyr::select(KMR, SUID, PredP_Ag =PredP, PredN_Ag = PredN, PredAll_Ag = PredAll, KhabRisk_Ag)
Pred_In_DF <- st_drop_geometry(Pred_In) %>% left_join(Khab_risk_In_DF, by = c("KMR", "SUID")) %>% dplyr::select(KMR, SUID, PredP_In =PredP, PredN_In = PredN, PredAll_In = PredAll, KhabRisk_In)
Pred_Fo_DF <- st_drop_geometry(Pred_Fo) %>% left_join(Khab_risk_Fo_DF, by = c("KMR", "SUID")) %>% dplyr::select(KMR, SUID, PredP_Fo =PredP, PredN_Fo = PredN, PredAll_Fo = PredAll, KhabRisk_Fo)
## Check if the number of rows match
nrow(Pred_In_DF) == nrow(Khab_risk_In_DF)
nrow(Pred_Ag_DF) == nrow(Khab_risk_Ag_DF)
nrow(Pred_Fo_DF) == nrow(Khab_risk_Fo_DF)

## Join all predictions and Koala habitat loss risk
Pred_ALL_DF <- full_join(Pred_In_DF, Pred_Ag_DF, by = c("KMR", "SUID")) %>%
    full_join(Pred_Fo_DF, by = c("KMR", "SUID")) %>%
    ## Calculate overall prediction and Koala habitat loss risk across all three clearing types
    mutate(PredAll_Ag_n = ifelse(is.na(PredAll_Ag), 1, 1 - PredAll_Ag),
           PredAll_In_n = ifelse(is.na(PredAll_In), 1, 1 - PredAll_In),
           PredAll_Fo_n = ifelse(is.na(PredAll_Fo), 1, 1 - PredAll_Fo),
           KhabRisk_Ag_n = ifelse(is.na(KhabRisk_Ag), 1, 1 - KhabRisk_Ag),
           KhabRisk_In_n = ifelse(is.na(KhabRisk_In), 1, 1 - KhabRisk_In),
           KhabRisk_Fo_n = ifelse(is.na(KhabRisk_Fo), 1, 1 - KhabRisk_Fo),
           Pred_All = 1 - (PredAll_In_n * PredAll_Ag_n * PredAll_Fo_n),
           KhabRisk_All = 1 - (KhabRisk_Ag_n * KhabRisk_In_n * KhabRisk_Fo_n)) %>% 
    dplyr::select(KMR, SUID, Pred_In = PredAll_In, Pred_Ag = PredAll_Ag, Pred_Fo = PredAll_Fo, 
                  KhabRisk_Ag, KhabRisk_In, KhabRisk_Fo, KhabRisk_All, Pred_All)

## Join with spatial units data
SUs_Pred_SF <- left_join(SUs_ZStats %>% dplyr::select(KMR, SUID, SA1), Pred_ALL_DF, by = c("KMR", "SUID"))

## Save the output to a geopackage
gpkg_path <- file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg")
if(!dir.exists(gpkg_path)){dir.create(gpkg_path)}

st_write(SUs_Pred_SF, dsn = gpkg_path, layer = "SUs_Predictions", append = FALSE)