## Plotting model selection results

#Load library & functions ----
library(tidyverse)
library(qs)
library(sf)
library(spdep)
library(INLA)
library(ggpubr)
library(colorspace)
library(ggspatial)
library(ggnewscale)
library(ggrepel)
library(cowplot)
library(patchwork)
library(shadowtext)
theme_set(theme_pubr())
source("functions.R")

# Prepare data ----

# Example usage for '_Ag', '_Fo', and '_In'
Cov_ls_Ag_long <- Get_cov_coeff_long(1)
Cov_ls_In_long <- Get_cov_coeff_long(2)
Cov_ls_Fo_long <- Get_cov_coeff_long(3)


## Agriculture ----
# 
# SUs_Ag <- qread("output/spatial_units/sus_ag.qs")
# SA1s <- qread("output/spatial_units/sa1s.qs")
# ZStats_Woody_Ag <- qread("output/data/ZStats_Woody_Ag.qs")
# ZStats_Covs_Ag <- qread("output/data/ZStats_Covs_Ag.qs")
# KMRs <- names(ZStats_Covs_Ag)
# kmr <- KMRs[1]
# Model_Ag <- fit_model2(KMR = kmr, ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# Cov_ls <- summary(Model_Ag$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% arrange(Covariate)
# 
# Cov_ls_Ag_long <- data.frame()
# for (kmr in KMRs){
#   
#   # read model selection results
#   MODEL <- qread(paste0("output/models/Model_", kmr, "_Ag.qs"))
# 
#   # Get Covariates and coefficients
#   Cov <- summary(MODEL$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
#   
#   # Separate the coefficient into positive or negative
#   Cof_PModel <- summary(MODEL$PModel)$fixed[,1]
#   Cof_PModel <- if_else(Cof_PModel>0, 1, -1)
#   Cof_NModel <- summary(MODEL$NModel)$fixed[,1]
#   Cof_NModel <- if_else(Cof_NModel>0, 1, -1)
#   
#   # Combine Covariate and coefficients
#   Cov_cof <- as.data.frame(cbind(Cov, Cof_PModel, Cof_NModel))
#   rownames(Cov_cof) <- NULL
#   Cov_ls_Ag <- left_join(Cov_ls, Cov_cof, by  = join_by("Covariate" == "Cov")) %>%  mutate(kmr = kmr)
#   Cov_ls_Ag_long <- rbind(Cov_ls_Ag_long, Cov_ls_Ag)
# }
# 
# # Change all NA to 0
# Cov_ls_Ag_long <- Cov_ls_Ag_long %>% 
#   mutate(Cof_PModel_Ag = ifelse(is.na(Cof_PModel), 0, Cof_PModel),
#          Cof_NModel_Ag = ifelse(is.na(Cof_NModel), 0, Cof_NModel))
# 
# # Split Pmodel and Nmodel into 2 dataframes
# Cov_ls_Ag_long_PMod <- Cov_ls_Ag_long %>%
#   dplyr::select(Covariate, Cof_PModel, kmr) %>% 
#   dplyr::mutate(Cof_PModel = ifelse(is.na(Cof_PModel), 0, Cof_PModel))
# 
# Cov_ls_Ag_long_NMod <- Cov_ls_Ag_long %>%
#   dplyr::select(Covariate, Cof_NModel, kmr) %>% 
#   
#   # duplicate the each rows 3 times and assign x and y coordinates for plotting triangles 
#   uncount(weight = 3) %>% 
#   mutate(x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c,
#          Cof_NModel = ifelse(is.na(Cof_NModel), 0, Cof_NModel))
# 
# # diverge_hcl(5, palette = "Blue_Red3")[4:5] 
# # diverge_hcl(5, palette = "Red_Green")[4:5]
# # diverge_hcl(5, palette = "Purple_Brown")[4:5]
# # Test plot
# # ggplot() +
# #   geom_tile(data = Cov_ls_Ag_long_PMod, aes(x = kmr, y = Covariate, fill = Cof_PModel))+
# #   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
# #   geom_polygon(data = Cov_ls_Ag_long_NMod, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel), color = "grey80")+
# #   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey80") ,
# #                     breaks = c(-1, 1, 0),
# #                     labels = c("negative", "positive", "Not Selected"))
# 
# 
# ## Forestry ----
# SUs_Fo <- qread("output/spatial_units/sus_Fo.qs")
# SA1s <- qread("output/spatial_units/sa1s.qs")
# ZStats_Woody_Fo <- qread("output/data/ZStats_Woody_Fo.qs")
# ZStats_Covs_Fo <- qread("output/data/ZStats_Covs_Fo.qs")
# KMRs <- names(ZStats_Covs_Fo)
# kmr <- KMRs[1]
# Model_Fo <- fit_model2(KMR = kmr, ClearType = 3, SpatUnits = SUs_Fo, RespData = ZStats_Woody_Fo, CovsCD = ZStats_Covs_Fo, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# Cov_ls <- summary(Model_Fo$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% arrange(Covariate)
# # kmr="CC"
# 
# Cov_ls_Fo_long <- data.frame()
# 
# for (kmr in KMRs){
#   
#   # read model selection results
#   SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_BC.qs"))
#   SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_Fo_FC.qs"))
#   
#   # Get Covariates based on DIC
#   MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
#   MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
#   MinDIC_BC_DIC <- unlist(MinDIC_BC)
#   MinDIC_FC_DIC <- unlist(MinDIC_FC)
#   
#   # Assign best model based on DIC
#   if(MinDIC_BC_DIC < MinDIC_FC_DIC){
#     Best_Mod <- SelModel_BC} else{
#       Best_Mod <- SelModel_FC}
#   cat("Best Model for  ", kmr, ": ", if_else(MinDIC_BC_dDIC < MinDIC_FC_dDIC, "BC", "FC"), "\n")
#   
#   # Get Covariates and coefficients
#   Cov <- summary(Best_Mod$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
#   
#   # Separate the coefficient into positive or negative
#   Cof_PModel <- summary(Best_Mod$PModel)$fixed[,1]
#   Cof_PModel <- if_else(Cof_PModel>0, 1, -1)
#   Cof_NModel <- summary(Best_Mod$NModel)$fixed[,1]
#   Cof_NModel <- if_else(Cof_NModel>0, 1, -1)
#   
#   # Combine Covariate and coefficients
#   Cov_cof <- as.data.frame(cbind(Cov, Cof_PModel, Cof_NModel))
#   rownames(Cov_cof) <- NULL
#   Cov_ls_Fo <- left_join(Cov_ls, Cov_cof, by  = join_by("Covariate" == "Cov")) %>%  mutate(kmr = kmr)
#   Cov_ls_Fo_long <- rbind(Cov_ls_Fo_long, Cov_ls_Fo)
# }
# 
# # Change all NA to 0
# Cov_ls_Fo_long <- Cov_ls_Fo_long %>% 
#   mutate(Cof_PModel_Fo = ifelse(is.na(Cof_PModel), 0, Cof_PModel),
#          Cof_NModel_Fo = ifelse(is.na(Cof_NModel), 0, Cof_NModel))
# 
# # Cov_ls_Fo_long_PMod <- Cov_ls_Fo_long %>%
# #   dplyr::select(Covariate, Cof_PModel, kmr) %>% 
# #   dplyr::mutate(Cof_PModel = ifelse(is.na(Cof_PModel), 0, Cof_PModel))
# # 
# # Cov_ls_Fo_long_NMod <- Cov_ls_Fo_long %>%
# #   dplyr::select(Covariate, Cof_NModel, kmr) %>% 
# #   uncount(weight = 3) %>% 
# #   mutate(x = as.numeric(as.factor(kmr)),
# #          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
# #          y = as.numeric(as.factor(Covariate)),
# #          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
# #          x1 = x + x_c, y1 = y + y_c,
# #          Cof_NModel = ifelse(is.na(Cof_NModel), 0, Cof_NModel))
# 
# ## Infrastructure ----
# SUs_In <- qread("output/spatial_units/sus_In.qs")
# SA1s <- qread("output/spatial_units/sa1s.qs")
# ZStats_Woody_In <- qread("output/data/ZStats_Woody_In.qs")
# ZStats_Covs_In <- qread("output/data/ZStats_Covs_In.qs")
# KMRs <- names(ZStats_Covs_In)
# kmr <- KMRs[1]
# Model_In <- fit_model2(KMR = kmr, ClearType = 2, SpatUnits = SUs_In, RespData = ZStats_Woody_In, CovsCD = ZStats_Covs_In, SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL)
# Cov_ls <- summary(Model_In$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% arrange(Covariate)
# # kmr="CC"
# 
# Cov_ls_In_long <- data.frame()
# 
# for (kmr in KMRs){
#   
#   # read model selection results
#   SelModel_BC <- qread(paste0("output/models/SelModel_", kmr, "_In_BC.qs"))
#   SelModel_FC <- qread(paste0("output/models/SelModel_", kmr, "_In_FC.qs"))
#   
#   # Get Covariates based on DIC
#   MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
#   MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
#   MinDIC_BC_DIC <- unlist(MinDIC_BC)
#   MinDIC_FC_DIC <- unlist(MinDIC_FC)
#   
#   # Assign best model based on DIC
#   if(MinDIC_BC_DIC < MinDIC_FC_DIC){
#     Best_Mod <- SelModel_BC} else{
#       Best_Mod <- SelModel_FC}
#   cat("Best Model for  ", kmr, ": ", if_else(MinDIC_BC_dDIC < MinDIC_FC_dDIC, "BC", "FC"), "\n")
#   
#   # Get Covariates and coefficients
#   Cov <- summary(Best_Mod$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% dplyr::select(Covariate) %>% unlist()
#   
#   # Separate the coefficient into positive or negative
#   Cof_PModel <- summary(Best_Mod$PModel)$fixed[,1]
#   Cof_PModel <- if_else(Cof_PModel>0, 1, -1)
#   Cof_NModel <- summary(Best_Mod$NModel)$fixed[,1]
#   Cof_NModel <- if_else(Cof_NModel>0, 1, -1)
#   
#   # Combine Covariate and coefficients
#   Cov_cof <- as.data.frame(cbind(Cov, Cof_PModel, Cof_NModel))
#   rownames(Cov_cof) <- NULL
#   Cov_ls_In <- left_join(Cov_ls, Cov_cof, by  = join_by("Covariate" == "Cov")) %>%  mutate(kmr = kmr)
#   Cov_ls_In_long <- rbind(Cov_ls_In_long, Cov_ls_In)
# }
# 
# # Change all NA to 0
# Cov_ls_In_long <- Cov_ls_In_long %>% 
#   mutate(Cof_PModel_In = ifelse(is.na(Cof_PModel), 0, Cof_PModel),
#          Cof_NModel_In = ifelse(is.na(Cof_NModel), 0, Cof_NModel))


# Combine results from all clearing types ----
Cov_df <- full_join(Cov_ls_Ag_long, Cov_ls_Fo_long, by = c("Covariate", "kmr")) %>% 
  full_join(Cov_ls_In_long, by = c("Covariate", "kmr")) %>% 
  filter(!Covariate == "LandUse5") %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, Cof_PModel_Fo, Cof_NModel_Fo, Cof_PModel_In, Cof_NModel_In) %>% 
  arrange(kmr, Covariate)
qsave(Cov_df, "output/data/Cov_df.qs")
Cov_df <- qread("output/data/Cov_df.qs")

# Checking for patterns for KMR along the coast
# Cov_df <- Cov_df %>% 
#   filter(kmr %in% c("CC", "NC", "SC"))
Cov_df2 <- Cov_df %>%
  uncount(weight = 3) %>%
  mutate(x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

# Create long format for plotting
Cov_df_long <- rbind(Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                     Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                     Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In))

Cov_df2_long <- rbind(Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x1, y1) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                      Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x1, y1) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                      Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In, x1, y1) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In))

qsave(Cov_df_long, "output/data/Cov_df_long_ForPlotting.qs")
qsave(Cov_df2_long, "output/data/Cov_df2_long_ForPlotting.qs")

Cov_df_long <- qread("output/data/Cov_df_long_ForPlotting.qs")
Cov_df2_long <- qread("output/data/Cov_df2_long_ForPlotting.qs")

# unique(Cov_df_long[,1])
# write.csv(unique(Cov_df_long[,1]), "output/data/Covariate_List.csv")

Cov_plot <- ggplot() +
  geom_tile(data = Cov_df_long, aes(x = kmr, y = Covariate, fill = Cof_PModel), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_long, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel), color = "grey80")+
  facet_wrap(vars(Model))+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  scale_y_discrete(labels = c("(Intercept)", 
                              "Agricultural profit", 
                              "Property Size", 
                              "Distance to urban center", 
                              "Distance to road", 
                              "Drought", 
                              "Ecological condition", 
                              "Fire", 
                              expression("Land tenure\n(Leasehold)"), 
                              expression("Land tenure\n(Crown purposes)"), 
                              "Land tenure (Other crown land)", 
                              "Land use (Production-Natural Env)", 
                              "Land use (Production-Dryland)", 
                              "Land use (Production-Irrigated)", 
                              "Land use (Intensive Uses)", 
                              "NVR Cat1 & Cat2 regulated", 
                              "NVR Cat2 Vulnerable / Sensitive", 
                              "Planning zone (Environment)", 
                              "Planning zone (Others)", 
                              "Planning zone (Residential)", 
                              "Planning zone (Rural)", 
                              "Population density", 
                              "Population Growth", 
                              "Rainfall", 
                              "Property value", 
                              expression("Socio-Economic PC1\n(Lower income)"), 
                              expression("Socio-Economic PC2\n(% Australia born, Eng. speaking)"), 
                              expression("Socio-Economic PC3\n(% 1-parent fam. with/without children u15)"), 
                              expression("Socio-Economic PC4\n(% coup. fam. no children u15, large household)"), 
                              expression("Socio-Economic PC5\n(% coup. fam. with children under 15)"), 
                              "Slope", 
                              expression("Soil PC1\n(High bulk density, sand content)"), 
                              expression("Soil PC2\n(High organic carbon, silt content)"), 
                              expression("Soil PC3\n(High total nitrogen, avail. water capacity)"), 
                              "Temperature"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
Cov_plot

ggsave("output/figures/Cov_plot.png", Cov_plot, width = 16, height = 17.5, dpi = 300)


# Extract numbers for results  -----

## Get the covariate that is selected most times across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(.)),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo)+ abs(Cof_PModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariates that is selected most times and negatively associated across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(if_else(.=="-1", "1", "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_NModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_NModel_Fo) + abs(Cof_PModel_In) + abs(Cof_NModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariates that is selected most times and positively associated with the prob of clearing (PMod) across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(if_else(.=="1", "1", "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariates that is selected most times and negatively associated with the prob of clearing (PMod) across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(if_else(.=="-1", "1", "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariates that is selected most times and positvely associated with the prob of clearing given clearing (NMod) across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(if_else(.=="1", "1", "0"))),
         Cof_total = abs(Cof_NModel_Ag) + abs(Cof_NModel_Fo) + abs(Cof_NModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariates that is selected most times and negatively associated with the prob of clearing given clearing (NMod) across all KMR and Clear Type
Cov_df %>% 
  mutate(across(3:8, ~as.integer(if_else(.=="-1", "1", "0"))),
         Cof_total = abs(Cof_NModel_Ag) + abs(Cof_NModel_Fo) + abs(Cof_NModel_In)) %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

# Covariate selection barplot ----

## Manipulate data for plotting
Cov_df_sum <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cov_Ag = Cof_PModel_Ag, Cov_Fo = Cof_PModel_Fo, Cov_In = Cof_PModel_In) %>%
  mutate(Cov_Ag = abs(as.numeric(Cov_Ag)),
         Cov_Fo = abs(as.numeric(Cov_Fo)),
         Cov_In = abs(as.numeric(Cov_In)),
         Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Property Size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire",
                               Covariate  == "LandTen2" ~ "Land tenure",
                               Covariate  == "LandTen3" ~ "Land tenure",
                               Covariate  == "LandTen4" ~ "Land tenure",
                               Covariate  == "LandUse1" ~ "Land use",
                               Covariate  == "LandUse2" ~ "Land use",
                               Covariate  == "LandUse3" ~ "Land use",
                               Covariate  == "LandUse4" ~ "Land use",
                               Covariate  == "NatVegReg1" ~ "Land use regulations",
                               Covariate  == "NatVegReg2" ~ "Land use regulations",
                               Covariate  == "PlanZone1" ~ "Planning zone",
                               Covariate  == "PlanZone2" ~ "Planning zone",
                               Covariate  == "PlanZone3" ~ "Planning zone",
                               Covariate  == "PlanZone4" ~ "Planning zone",
                               Covariate  == "PopDen" ~ "Population density",
                               Covariate  == "PopGro" ~ "Population Growth",
                               Covariate  == "Precip" ~ "Rainfall",
                               Covariate  == "PropVal" ~ "Property value",
                               Covariate  == "ScEc_PC1" ~ "Socio-Economic PC1",
                               Covariate  == "ScEc_PC2" ~ "Socio-Economic PC2",
                               Covariate  == "ScEc_PC3" ~ "Socio-Economic PC3",
                               Covariate  == "ScEc_PC4" ~ "Socio-Economic PC4",
                               Covariate  == "ScEc_PC5" ~ "Socio-Economic PC5",
                               Covariate  == "slope" ~ "Slope",
                               Covariate  == "Soil_PC1" ~ "Soil PC1",
                               Covariate  == "Soil_PC2" ~ "Soil PC2",
                               Covariate  == "Soil_PC3" ~ "Soil PC3",
                               Covariate  == "Temp" ~ "Temperature",
                               .default = Covariate)) %>% 
  distinct()

## Manipulate Agriculture data for plotting
Cov_df_sum_Ag <- Cov_df_sum %>% 
  dplyr::select(Covariate, kmr, Cov_Ag) %>% 
  drop_na() %>% 
  group_by(Covariate) %>%
  summarise(Cov_Ag = sum(Cov_Ag)) %>% 
  arrange(desc(Cov_Ag)) %>% 
  filter(Covariate != "(Intercept)")

## Manipulate Forestry data for plotting
Cov_df_sum_Fo <- Cov_df_sum %>% 
  dplyr::select(Covariate, kmr, Cov_Fo) %>% 
  drop_na() %>% 
  group_by(Covariate) %>%
  summarise(Cov_Fo = sum(Cov_Fo)) %>% 
  arrange(desc(Cov_Fo)) %>% 
  filter(Covariate != "(Intercept)")

## Manipulate Infrastructure data for plotting
Cov_df_sum_In <- Cov_df_sum %>% 
  dplyr::select(Covariate, kmr, Cov_In) %>% 
  drop_na() %>% 
  group_by(Covariate) %>%
  summarise(Cov_In = sum(Cov_In)) %>% 
  arrange(desc(Cov_In)) %>% 
  filter(Covariate != "(Intercept)")

## Agriculture Covariate selection barplot
Cov_sum_Ag_plot <- ggplot(Cov_df_sum_Ag, aes(y = reorder(Covariate, Cov_Ag), x = Cov_Ag)) +
  geom_bar(stat = "identity") +
  labs(title = "Agriculture", x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y=element_blank(),
        plot.title = element_text(hjust=0.9 , vjust = -10))
Cov_sum_Ag_plot

## Forestry Covariate selection barplot
Cov_sum_Fo_plot <- ggplot(Cov_df_sum_Fo, aes(y = reorder(Covariate, Cov_Fo), x = Cov_Fo)) +
  geom_bar(stat = "identity") +
  labs(title = "Forestry", x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x=element_blank())
Cov_sum_Fo_plot

## Infrastructure Covariate selection barplot
Cov_sum_In_plot <- ggplot(Cov_df_sum_In, aes(y = reorder(Covariate, Cov_In), x = Cov_In)) +
  geom_bar(stat = "identity") +
  labs(title = "Infrastructure", y = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x=element_blank())
Cov_sum_In_plot

## Combine all plots
Cov_sum_plot <- ggarrange(Cov_sum_Ag_plot+rremove("xlab")+rremove("ylab"), 
                          Cov_sum_Fo_plot+rremove("xlab")+rremove("ylab"), 
                          Cov_sum_In_plot+rremove("xlab")+rremove("ylab"), 
                          ncol = 3, nrow = 1)

## Annotate plots
Cov_sum_plot <- annotate_figure(Cov_sum_plot, bottom = text_grob("Number of times covariate selected"))
Cov_sum_plot
## export barplot
ggsave("output/figures/Cov_sum_plot.png", Cov_sum_plot2, width = 11, height = 6, dpi = 300, bg = "white")

## Barplot to show positive and negatively associated covariates. (Not used)

Cov_df_sum2 <- Cov_df %>% 
  mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Property Size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire",
                               Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
                               Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
                               Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
                               Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
                               Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
                               Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
                               Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
                               Covariate  == "NatVegReg1" ~ "NVR Cat1 & Cat2 regulated",
                               Covariate  == "NatVegReg2" ~ "NVR Cat2 Vulnerable / Sensitive",
                               Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
                               Covariate  == "PlanZone2" ~ "Planning zone (Others)",
                               Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
                               Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
                               Covariate  == "PopDen" ~ "Population density",
                               Covariate  == "PopGro" ~ "Population Growth",
                               Covariate  == "Precip" ~ "Rainfall",
                               Covariate  == "PropVal" ~ "Property value",
                               Covariate  == "ScEc_PC1" ~ "Socio-Economic PC1",
                               Covariate  == "ScEc_PC2" ~ "Socio-Economic PC2",
                               Covariate  == "ScEc_PC3" ~ "Socio-Economic PC3",
                               Covariate  == "ScEc_PC4" ~ "Socio-Economic PC4",
                               Covariate  == "ScEc_PC5" ~ "Socio-Economic PC5",
                               Covariate  == "slope" ~ "Slope",
                               Covariate  == "Soil_PC1" ~ "Soil PC1",
                               Covariate  == "Soil_PC2" ~ "Soil PC2",
                               Covariate  == "Soil_PC3" ~ "Soil PC3",
                               Covariate  == "Temp" ~ "Temperature",
                               .default = Covariate)) %>% 
  distinct() %>% 
  group_by(Covariate) %>%
  summarise(Cof_PModel_Ag_pos = sum(as.numeric(Cof_PModel_Ag == 1)),
            Cof_PModel_Ag_neg = sum(as.numeric(Cof_PModel_Ag == -1))*-1,
            Cof_NModel_Ag_pos = sum(as.numeric(Cof_NModel_Ag == 1)),
            Cof_NModel_Ag_neg = sum(as.numeric(Cof_NModel_Ag == -1))*-1,
            Cof_PModel_Ag_sum = sum(as.numeric(Cof_PModel_Ag == 1)) + sum(as.numeric(Cof_PModel_Ag == -1)),
            Cof_PModel_Fo_pos = sum(as.numeric(Cof_PModel_Fo == 1)),
            Cof_PModel_Fo_neg = sum(as.numeric(Cof_PModel_Fo == -1))*-1,
            Cof_NModel_Fo_pos = sum(as.numeric(Cof_NModel_Fo == 1)),
            Cof_NModel_Fo_neg = sum(as.numeric(Cof_NModel_Fo == -1))*-1,
            Cof_PModel_Fo_sum = sum(as.numeric(Cof_PModel_Fo == 1)) + sum(as.numeric(Cof_PModel_Fo == -1)),
            Cof_PModel_In_pos = sum(as.numeric(Cof_PModel_In == 1)),
            Cof_PModel_In_neg = sum(as.numeric(Cof_PModel_In == -1))*-1,
            Cof_NModel_In_pos = sum(as.numeric(Cof_NModel_In == 1)),
            Cof_NModel_In_neg = sum(as.numeric(Cof_NModel_In == -1))*-1,
            Cof_PModel_In_sum = sum(as.numeric(Cof_PModel_In == 1)) + sum(as.numeric(Cof_PModel_In == -1)))

Cov_df_sum2 %>% arrange(desc(Cof_PModel_Fo_sum)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
Cov_df_sum2 %>% arrange(Cof_PModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
Cov_df_sum2 %>% arrange(desc(Cof_PModel_Fo_pos)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
Cov_df_sum2 %>% arrange(Cof_PModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
Cov_df_sum2 %>% arrange(desc(Cof_NModel_Fo_pos)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
Cov_df_sum2 %>% arrange(Cof_NModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)

Cov_df_sum2 %>% arrange(desc(Cof_PModel_In_sum)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
Cov_df_sum2 %>% arrange(Cof_PModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
Cov_df_sum2 %>% arrange(desc(Cof_PModel_In_pos)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
Cov_df_sum2 %>% arrange(Cof_PModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
Cov_df_sum2 %>% arrange(desc(Cof_NModel_In_pos)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
Cov_df_sum2 %>% arrange(Cof_NModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)

# Cov_df_sum_Ag2 <- Cov_df_sum2 %>% 
#   mutate(Cof = abs(Cof_PModel_Ag_pos) + abs(Cof_PModel_Ag_neg)) %>%
#   arrange(desc(Cof)) %>% 
#   dplyr::select(Covariate,  Cof_PModel_Ag_pos,  Cof_PModel_Ag_neg,  Cof_NModel_Ag_pos,  Cof_NModel_Ag_neg) %>% 
#   drop_na()

Cov_df_sum_Ag2 <- rbind(Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_pos) %>% mutate(Model = "PModel"),
                        Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_neg) %>% mutate(Model = "PModel"),
                        Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_pos) %>% mutate(Model = "NModel"),
                        Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_neg) %>% mutate(Model = "NModel"))

Cov_df_sum_Ag2 <- rbind(Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_pos, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "PModel"),
                        Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_neg, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "PModel"),
                        Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_pos, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "NModel"),
                        Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_neg, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "NModel")) %>% 
  drop_na()

Cov_sum_Ag2_plot <- ggplot(data = Cov_df_sum_Ag2, aes(x = reorder(Covariate, -Cof_Ag_sum), y = Cof_Ag, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = "Agriculture",
       y = "Number of times Covariate selected") +
  scale_y_continuous(breaks = seq(-8, 8, by = 4),labels = abs(seq(-8, 8, by = 4)), guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust=0.9 , vjust = -10))
ggsave("output/figures/Cov_sum_Ag2_plot.png", Cov_sum_Ag2_plot, width = 11, height = 11, dpi = 300, bg = "white")
  
Cov_sum_Ag_plot <- ggplot(data = Cov_df_sum_Ag %>% 
                            dplyr::select(Covariate, Cof_PModel_Ag_pos, Cof_PModel_Ag_neg, Cof_NModel_Ag_pos,Cof_PModel_Ag_neg) %>% 
                            drop_na(), 
                          aes(x = Covariate), y = Cov_Ag) +
  geom_bar(stat = "identity") +
  labs(title = "Agriculture",
       y = "Number of times Covariate selected") +
  scale_y_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust=0.9 , vjust = -10))
Cov_sum_Ag_plot  


# Combined tile and bar plot ----
## Prepare data ----
Cov_df <- qread("output/data/Cov_df.qs") %>% 
  mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Property Size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire",
                               Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
                               Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
                               Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
                               Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
                               Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
                               Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
                               Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
                               Covariate  == "NatVegReg1" ~ "NVR Cat1 & Cat2 regulated",
                               Covariate  == "NatVegReg2" ~ "NVR Cat2 Vulnerable / Sensitive",
                               Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
                               Covariate  == "PlanZone2" ~ "Planning zone (Others)",
                               Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
                               Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
                               Covariate  == "PopDen" ~ "Population density",
                               Covariate  == "PopGro" ~ "Population Growth",
                               Covariate  == "Precip" ~ "Rainfall",
                               Covariate  == "PropVal" ~ "Property value",
                               Covariate  == "ScEc_PC1" ~ "Socio-Economic PC1 (Lower income)",
                               Covariate  == "ScEc_PC2" ~ "Socio-Economic PC2 (% Australia born, Eng. speaking)",
                               Covariate  == "ScEc_PC3" ~ "Socio-Economic PC3 (% 1-parent fam. with/without children u15)",
                               Covariate  == "ScEc_PC4" ~ "Socio-Economic PC4 (% coup. fam. no children u15)",
                               Covariate  == "ScEc_PC5" ~ "Socio-Economic PC5 (% coup. fam. with children under 15)",
                               Covariate  == "slope" ~ "Slope",
                               Covariate  == "Soil_PC1" ~ "Soil PC1 (High bulk density, sand content)",
                               Covariate  == "Soil_PC2" ~ "Soil PC2 (High organic carbon, silt content)",
                               Covariate  == "Soil_PC3" ~ "Soil PC3 (High total nitrogen, avail. water capacity)",
                               Covariate  == "Temp" ~ "Temperature",
                               .default = Covariate))

## Agriculture ----
Cov_df_Ag <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df2_Ag <- Cov_df_Ag %>% 
  uncount(weight = 3) %>% 
  mutate(x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Cov_df_Ag_sum <- Cov_df_Ag %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Ag))))

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_Ag, aes(x = kmr, y = Covariate, fill = Cof_PModel_Ag), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_Ag, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Ag), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(axis.title.y=element_blank())
Tile_plot

Bar_plot <- ggplot(data = Cov_df_Ag_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
Bar_plot
TileBar_plot_Ag <- Tile_plot + Bar_plot + plot_layout(width = c(5,3.5))
TileBar_plot_Ag
ggsave("output/figures/TileBar_plot_Ag.png", TileBar_plot_Ag, width = 11, height = 8, dpi = 300, bg = "white")

## Forestry ----
Cov_df_Fo <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df2_Fo <- Cov_df_Fo %>% 
  uncount(weight = 3) %>% 
  mutate(x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Cov_df_Fo_sum <- Cov_df_Fo %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Fo))))

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_Fo, aes(x = kmr, y = Covariate, fill = Cof_PModel_Fo), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_Fo, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Fo), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(axis.title.y=element_blank())
Tile_plot

Bar_plot <- ggplot(data = Cov_df_Fo_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
Bar_plot
TileBar_plot_Fo <- Tile_plot + Bar_plot + plot_layout(width = c(5,3.5))
ggsave("output/figures/TileBar_plot_Fo.png", TileBar_plot_Fo, width = 11, height = 8, dpi = 300, bg = "white")

## Infrastructure ----
Cov_df_In <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df2_In <- Cov_df_In %>% 
  uncount(weight = 3) %>% 
  mutate(x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Cov_df_In_sum <- Cov_df_In %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_In))))

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_In, aes(x = kmr, y = Covariate, fill = Cof_PModel_In), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_In, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_In), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(axis.title.y=element_blank())
Tile_plot

Bar_plot <- ggplot(data = Cov_df_In_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
Bar_plot
TileBar_plot_In <- Tile_plot + Bar_plot + plot_layout(width = c(5,3.5))
ggsave("output/figures/TileBar_plot_In.png", TileBar_plot_In, width = 11, height = 8, dpi = 300, bg = "white")


# Map for study area (KMR) ----

# Load spatial units for defining the study area 
SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")

# Load ABS urban areas shapefile
ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
  st_transform(st_crs(SUs_Ag$CC))

# Load Koala Modelling Regions shapefile, include a buffer zone
KMR_shp <- st_read("input/spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp") %>% 
  # calculate centroid coordinates for labeling
  mutate(x = st_coordinates(st_centroid(.))[,1], y = st_coordinates(st_centroid(.))[,2],
         # manipulate labels 
         KMR = c("(NC)", "(CC)", "(SC)", "(NT)", "(NS)", "(CST)", "(DRP)", "(FW)", "(R)"), KMRname2 = str_wrap(paste(KMRname, KMR, sep = " "), width = 6),
         # Modify coordinates for NS and NT for better visualisation
         y = if_else(KMR == "(NS)", y-100000, y), y = if_else(KMR == "(NT)", y+70000, y))

# Load state boundary shapefile
STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
  st_transform(st_crs(SUs_Ag$CC)) %>% st_crop(KMR_shp)

# Select urban areas for labelling
NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Lismore", "Port Macquarie",    #NC
                           "Sydney",     #CC
                           "Nowra - Bomaderry", "Bega",    #SC
                           "Armidale", "Tamworth",    #NT
                           "Narrabri", "Dubbo",    #NS
                           "Orange", "Wagga Wagga",    #CST
                           "Griffith", "Deniliquin",    #R
                           "Broken Hill", "Cobar")) %>%    #FW
  dplyr::select(UCL_NAME16, geometry) %>% 
  distinct(UCL_NAME16, .keep_all = TRUE) %>% 
  st_centroid(.) %>% 
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])

## Plot map----
KMR_map <- ggplot()+
  # State boundary grey background
  geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
  # Koala Modelling Regions
  geom_sf(data = KMR_shp, fill = "lightgoldenrod1", color = "grey10", lwd = 0.2, linetype= "dotdash" )+
  geom_shadowtext(data = KMR_shp, aes(x = x, y = y, label = KMRname2), 
                  fontface = "bold", nudge_y = -5, size = 5,
                  color = "black",     # text color
                  bg.color = "white", # shadow color
                  bg.r = 0.05)+
  # Urban areas points
  geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
  geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), nudge_y = -5, size = 4,
                  color = "grey50",     # text color
                  bg.color = "white", # shadow color
                  bg.r = 0.05)+          # shadow radius
  # North arrow and scale bar
  ggspatial::annotation_scale(location = "bl", pad_y = unit(1, "cm"), pad_x = unit(1, "cm"))+
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(1.5, "cm"), pad_x = unit(1, "cm"))+
  # Theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], 
           # datum = st_crs(KMR_shp),
           expand = FALSE)

# Export KMR map
ggsave("output/figures/KMR_map.png", KMR_map, width = 11, height = 11, dpi = 300)


# Deforestation Risk map ----

SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")

Pred_Ag <- qread("output/predictions/Pred_Ag.qs")
Pred_Fo <- qread("output/predictions/Pred_Fo.qs")
Pred_In <- qread("output/predictions/Pred_In.qs")

### Extract numbers for results ----
Pred_Ag_all_result <- Pred_Ag %>% mutate(Area = as.numeric(st_area(.)/1e4),
                                         RemWoodyHa = if_else( (Woody - Wdy_Clr) >0 , Woody - Wdy_Clr , 0) * 0.0625)
Pred_Fo %>% st_drop_geometry() %>% filter(is.na(PredAll))
sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.5,])

sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.25,])

sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.1,])

Pred_Fo_all_result <- do.call(rbind, Pred_Fo) %>% 
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else(NT-R>0, NT-R, 0) * 0.0625)
sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.5,])

sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.25,])

sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.1,])

Pred_In_all_result <- do.call(rbind, Pred_In) %>% 
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else(NT-R>0, NT-R, 0) * 0.0625)
sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.5,])

sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.25,])

sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.1,])

## Plot maps----
### load base layers and target data  layers
#### Data Layers
Pred_Ag <- qread("output/predictions/Pred_Ag.qs")
Pred_Fo <- qread("output/predictions/Pred_Fo.qs")
Pred_In <- qread("output/predictions/Pred_In.qs")

#### Base Layers
KMR_shp <- st_read("input/spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")

ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
  st_transform(st_crs(Pred_Ag))

STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
  st_transform(st_crs(Pred_Ag)) %>% 
  st_crop(KMR_shp)

NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Lismore", "Port Macquarie",
                         "Sydney",
                         "Nowra - Bomaderry", "Bega",
                         "Armidale", "Tamworth",
                         "Narrabri", "Dubbo",
                         "Orange", "Wagga Wagga",
                         "Griffith", "Deniliquin",
                         "Brewarrina (L)", 
                         "Broken Hill", "Ivanhoe (L)", "Cobar")) %>% 
  dplyr::select(UCL_NAME16, geometry) %>% 
  distinct(UCL_NAME16, .keep_all = TRUE) %>% 
  st_centroid(.) %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])

#### File names
FilenamePath_PNG_Ag <- "output/figures/Pred_Ag_map1.png"
FilenamePath_PNG_Fo <- "output/figures/Pred_Fo_map1.png"
FilenamePath_PNG_In <- "output/figures/Pred_In_map1.png"

PLOTMAP_risk_NSW(DATA = Pred_Ag, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 1, FilenamePath_PNG = "output/figures/Pred_Ag_map1.png")
PLOTMAP_risk_NSW(DATA = Pred_Fo, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 2, FilenamePath_PNG = "output/figures/Pred_Fo_map1.png")
PLOTMAP_risk_NSW(DATA = Pred_In, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 3, FilenamePath_PNG = "output/figures/Pred_In_map1.png")



# Pred_Ag_all <- Pred_Ag
# 
# # ggplot()+
# #   geom_point(data = NSW_urb2_pt, aes(x = x, y = y, color = "red3", size = 1))+
# #   geom_text(data = NSW_urb2_pt, aes(x = x, y = y, label = UCL_NAME16), size = 4)+
# #   geom_sf(data = KMR_shp, fill = NA, color = "black", lwd = 0.1)
# # 
# # st_bbox(KMR_shp)[c(1,3)]
# 
# Pred_Ag_map <- ggplot()+
#   
#   geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
#   
#   geom_sf(data = Pred_Ag_all, aes(fill = PredAll), color = NA)+
#   geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
#   scale_fill_gradientn(colours = hcl.colors(8, palette = "Blues 3" ,rev = TRUE), name = expression("Deforestation\nrisk"))+
#   
#   # start a new scale
#   new_scale_colour() +
#   
#   geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
#   geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), 
#                   fontface = "bold", nudge_y = -5, size = 3,
#                   color = "black",     # text color
#                   bg.color = "grey90", # shadow color
#                   bg.r = 0.05)+          # shadow radius
#   
#   ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
#   ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+
#   
#   theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#   theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(legend.position = c(0.9, 0.3))+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#   coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
# ggsave("output/figures/Pred_Ag_map1.png", Pred_Ag_map, width = 11, height = 11, dpi = 300, bg = "white")
# 
# 
# # Pred_Fo_all <- do.call(rbind, Pred_Fo)
# Pred_Fo_all <- Pred_Fo
# Pred_Fo_map <- ggplot()+
#   
#   geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
#   
#   geom_sf(data = Pred_Fo_all, aes(fill = PredAll), color = NA)+
#   geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
#   scale_fill_gradientn(colours = hcl.colors(8, palette = "Blues 3" ,rev = TRUE), name = expression("Deforestation\nrisk"))+
#   
#   # start a new scale
#   new_scale_colour() +
#   
#   geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
#   geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), 
#                   fontface = "bold", nudge_y = -5, size = 3,
#                   color = "black",     # text color
#                   bg.color = "grey90", # shadow color
#                   bg.r = 0.05)+          # shadow radius
#   
#   ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
#   ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+
#   
#   theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#   theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(legend.position = c(0.9, 0.3))+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#   coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
# ggsave("output/figures/Pred_Fo_map2.png", Pred_Fo_map, width = 11, height = 11, dpi = 300, bg = "white")
# 
# 
# 
# # Pred_In_all <- do.call(rbind, Pred_In)
# Pred_In_all <- Pred_In
# Pred_In_map <- ggplot() +
#   
#   geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
#   
#   geom_sf(data = Pred_In_all, aes(fill = PredAll), color = NA)+
#   geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
#   scale_fill_gradientn(colours = hcl.colors(8, palette = "Blues 3" ,rev = TRUE), name = expression("Deforestation\nrisk"))+
#   
#   # start a new scale
#   new_scale_colour() +
#   
#   geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
#   geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), 
#                   fontface = "bold", nudge_y = -5, size = 3,
#                   color = "black",     # text color
#                   bg.color = "grey90", # shadow color
#                   bg.r = 0.05)+          # shadow radius
#   
#   ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
#   ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+
#   
#   theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#   theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(legend.position = c(0.9, 0.3))+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#   coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
# ggsave("output/figures/Pred_In_map2.png", Pred_In_map, width = 11, height = 11, dpi = 300, bg = "white")
# 

# Koala Habitat loss risk----

# SUs_Ag <- qread("output/spatial_units/SUs_Ag.qs")
# 
# Khab_risk_Ag <- qread("output/predictions/Khab_risk_Ag.qs")
# Khab_risk_Fo <- qread("output/predictions/Khab_risk_Fo.qs")
# Khab_risk_In <- qread("output/predictions/Khab_risk_In.qs")
# 
# KMR_shp <- st_read("input/spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")
# 
# ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
#   st_transform(st_crs(SUs_Ag$CC))
# 
# STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
#   st_transform(st_crs(SUs_Ag$CC)) %>% 
#   st_crop(KMR_shp)
# 
# NSW_urb <- st_intersection(ABS_urb, KMR_shp) 
# NSW_urb2 <- NSW_urb %>%
#   filter(!str_detect(UCL_NAME16, "emain")) %>% 
#   group_by(KMRname) %>%
#   top_n(10, AREASQKM16)
# NSW_urb2_pt <- st_centroid(NSW_urb2)
# NSW_urb2_pt$x <- st_coordinates(NSW_urb2_pt)[,1]
# NSW_urb2_pt$y <- st_coordinates(NSW_urb2_pt)[,2]
# 
# NSW_urb_sel <- NSW_urb %>%
#   filter(UCL_NAME16 %in% c("Lismore", "Port Macquarie",
#                            "Sydney",
#                            "Nowra - Bomaderry", "Bega",
#                            "Armidale", "Tamworth",
#                            "Narrabri", "Dubbo",
#                            "Orange", "Wagga Wagga",
#                            "Griffith", "Deniliquin",
#                            "Brewarrina (L)", 
#                            "Broken Hill", "Ivanhoe (L)", "Cobar")) %>% 
#   dplyr::select(UCL_NAME16, geometry) %>% 
#   distinct(UCL_NAME16, .keep_all = TRUE)
# 
# NSW_urb_sel_pt <- st_centroid(NSW_urb_sel)
# NSW_urb_sel_pt$x <- st_coordinates(NSW_urb_sel_pt)[,1]
# NSW_urb_sel_pt$y <- st_coordinates(NSW_urb_sel_pt)[,2]
# 
# NSW_urb_pt <- st_centroid(NSW_urb) %>% 
#   mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% 
#   filter(!str_detect(UCL_NAME16, "emain"))
# 

## Plot maps----
### load base layers and target data  layers
#### Data Layers
Khab_risk_Ag <- qread("output/predictions/Khab_risk_Ag.qs")
Khab_risk_Fo <- qread("output/predictions/Khab_risk_Fo.qs")
Khab_risk_In <- qread("output/predictions/Khab_risk_In.qs")

#### Base Layers
KMR_shp <- st_read("input/spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")

ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
  st_transform(st_crs(Khab_risk_Ag))

STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
  st_transform(st_crs(Khab_risk_Ag)) %>% 
  st_crop(KMR_shp)

NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Lismore", "Port Macquarie",
                           "Sydney",
                           "Nowra - Bomaderry", "Bega",
                           "Armidale", "Tamworth",
                           "Narrabri", "Dubbo",
                           "Orange", "Wagga Wagga",
                           "Griffith", "Deniliquin",
                           "Brewarrina (L)", 
                           "Broken Hill", "Ivanhoe (L)", "Cobar")) %>% 
  dplyr::select(UCL_NAME16, geometry) %>% 
  distinct(UCL_NAME16, .keep_all = TRUE) %>% 
  st_centroid(.) %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])

NSW_urb_pt <- ABS_urb %>%
  filter(!str_detect(UCL_NAME16, "emain")) %>% 
  st_centroid() %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% 
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Mungindi (NSW Part) (L)" ~ "Mungindi (NSW Part)",
                                UCL_NAME16 == "Lismore" ~ "Lismore", 
                                UCL_NAME16 == "Bonalbo (L)" ~ "Bonalbo",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                UCL_NAME16 == "Jilliby (L)" ~ "Jilliby",
                                .default = UCL_NAME16))

NSW_urb_pt %>% filter(UCL_NAME16 %in% c("Oberon"))

Inset_BL_Ag <- data.frame(x = c(9376800, 9435500, 9464100), y = c(4726000, 4826500, 4911747))
Inset_BL_In <- data.frame(x = c(9599470, 9660498, 9769700), y = c(4352856, 4477888, 4593377))
Inset_BL_Fo <- data.frame(x = c(9779110, 9495375, 9327460), y = c(4903177, 4415825, 4196286))

#### File names
FilenamePath_PNG_Ag <- "output/figures/Khab_risk_Ag_map1.png"
FilenamePath_PNG_Fo <- "output/figures/Khab_risk_Fo_map1.png"
FilenamePath_PNG_In <- "output/figures/Khab_risk_In_map1.png"


Ag_risk_with_Insets<- PLOTMAP_risk_with_Insets(DATA = Khab_risk_Ag, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 1, 
                         Inset_BL = Inset_BL_Ag, Inset_dim = 100000, URB_PT_Main =NULL, URB_PT_SUB1 = "Coonamble", URB_PT_SUB2 = c("Collarenebri (L)", "Wee Waa"),  URB_PT_SUB3 = "Mungindi (NSW Part)", 
                         FilenamePath_PNG = FilenamePath_PNG_Ag, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

In_risk_with_Insets<- PLOTMAP_risk_with_Insets(DATA = Khab_risk_In, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 2,
                         Inset_BL = Inset_BL_In, Inset_dim = 100000, URB_PT_Main =NULL, 
                         URB_PT_SUB1 = c("Blue Mountains", "Sydney", "Galston", "The Oaks"), URB_PT_SUB2 = c("Singleton", "Newcastle", "Jilliby"),  URB_PT_SUB3 =  c("Port Macquarie", "Taree", "Forster - Tuncurry"),
                         FilenamePath_PNG = FilenamePath_PNG_In, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Fo_risk_with_Insets <- PLOTMAP_risk_with_Insets(DATA = Khab_risk_Fo, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 3,
                         Inset_BL = Inset_BL_Fo, Inset_dim = 100000, URB_PT_Main =NULL, URB_PT_SUB1 = c("Lismore", "Bonalbo"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
                         FilenamePath_PNG = FilenamePath_PNG_Fo, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

### Extract numbers for results ----
Khab_risk_Ag <- qread("output/predictions/Khab_risk_Ag.qs")
Khab_risk_Fo <- qread("output/predictions/Khab_risk_Fo.qs")
Khab_risk_In <- qread("output/predictions/Khab_risk_In.qs")

Khab_risk_Ag_all_result <- do.call(rbind, Khab_risk_Ag) %>% 
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.5,])

sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.25,])

sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.1,])

Khab_risk_Fo_all_result <- do.call(rbind, Khab_risk_Fo) %>% 
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.5,])

sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.25,])

sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.1,])

Khab_risk_In_all_result <- do.call(rbind, Khab_risk_In) %>% 
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.5,])

sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.25,])

sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.1,])