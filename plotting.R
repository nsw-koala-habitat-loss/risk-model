## Plotting model selection results

#Load library & functions ----
library(tidyverse)
library(qs)
library(qs2)
library(sf)
library(spdep)
library(INLA)
library(ggpubr)
library(colorspace)
library(ggspatial)
library(ggnewscale)
library(grid)
library(ggrepel)
library(cowplot)
library(patchwork)
library(shadowtext)
options(max.print = 50)
theme_set(theme_pubr())
source("functions.R")

# Prepare data ----
## Set directory
## Note: change the path to the input and output directory
INPUT_DIR <- file.path("R:/HABCLEAR22-Q5221/risk-model/input")
OUTPUT_DIR <- file.path("R:/HABCLEAR22-Q5221/risk-model/output")

INPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/input")
OUTPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/output")

# Extract numbers for results  -----


# Cov_Fullname_LU <- c(
#   AgProf = "Agricultural profit",
#   Area = "Parcel size",
#   DistCity = "Distance to urban center",
#   DistRoad = "Distance to road",
#   Drought = "Drought",
#   Drought1 = "Drought",
#   EcolCond = "Ecological condition",
#   Elev = "Elevation",
#   Fire = "Fire history",
#   Fire1 = "Fire history",
#   ForTen = "Forest tenure",
#   ForType = "Forest type",
#   LandTen = "Land tenure",
#   LandTen2 = "Land tenure (Leasehold)",
#   LandTen3 = "Land tenure (Crown purposes)",
#   LandTen4 = "Land tenure (Other crown land)",
#   LandUse = "Land use",
#   LandUse1 = "Land use (Production-Natural Env)",
#   LandUse2 = "Land use (Production-Dryland)",
#   LandUse3 = "Land use (Production-Irrigated)",
#   LandUse4 = "Land use  (Intensive Uses)",
#   NatVegReg = "Native vegetation regulation (NVR)",
#   NatVegReg1 = "NVR regulated land",
#   NatVegReg2 = "NVR regulated land",
#   PlanZone = "Planning zone",
#   PlanZone1 = "Planning zone (Environment)",
#   PlanZone2 = "Planning zone (Others)",
#   PlanZone3 = "Planning zone (Residential)",
#   PlanZone4 = "Planning zone (Rural)",
#   PolPref = "Political preference",
#   PopDen = "Population density",
#   PopGro = "Population Growth",
#   Precip = "Rainfall",
#   PropVal = "Land value",
#   Remoteness = "Remoteness",
#   ScEc_PC1 = "Socio-Economic PC1 (Lower income)",
#   ScEc_PC2 = "Socio-Economic PC2 (% Australia born, Eng. speaking)",
#   ScEc_PC3 = "Socio-Economic PC3 (% 1-parent fam. with/without children u15)",
#   ScEc_PC4 = "Socio-Economic PC4 (% coup. fam. no children u15)",
#   ScEc_PC5 = "Socio-Economic PC5 (% coup. fam. with children under 15)",
#   slope = "Slope",
#   Soil_PC1 = "Soil PC1 (High bulk density, sand content)",
#   Soil_PC2 = "Soil PC2 (High organic carbon, silt content)",
#   Soil_PC3 = "Soil PC3 (High total nitrogen, avail. water capacity)",
#   Temp = "Temperature",
#   TenType = "Forest tenure type"
# )

# # Checking for patterns for KMR along the coast
# # Cov_df <- Cov_df %>%
# #   filter(kmr %in% c("CC", "NC", "SC"))
# Cov_df2 <- Cov_df %>%
#   uncount(weight = 3) %>%
#   mutate(x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c)

# # Create long format for plotting
# Cov_df_long <- rbind(Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
#                      Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
#                      Cov_df %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In))

# Cov_df2_long <- rbind(Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x1, y1) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
#                       Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x1, y1) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
#                       Cov_df2 %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In, x1, y1) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In))

# qsave(Cov_df_long, file.path(OUTPUT_DIR, "data/Cov_df_long_ForPlotting.qs"))
# qsave(Cov_df2_long, file.path(OUTPUT_DIR, "data/Cov_df2_long_ForPlotting.qs"))

# Cov_df_long <- qread(file.path(OUTPUT_DIR, "data/Cov_df_long_ForPlotting.qs"))
# Cov_df2_long <- qread(file.path(OUTPUT_DIR, "data/Cov_df2_long_ForPlotting.qs"))


# Covariate selection barplot ----


####################################################################################
# PLot bar graph showing all covariates (standardized names) and all clearing types ----
Cov_df_all <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs")) %>%
  mutate(CovTypes = case_when(Covariate %in% c("DistRoad", "DistCity", "Area") ~ "Accessibility",
                           Covariate %in% c("Drought1", "Elevation", "Fire1",
                                            "ForType", "Precip", "slope",
                                            "Soil_PC1", "Soil_PC2", "Soil_PC3", "Temp") ~ "Biophysical",
                           Covariate %in% c("PolPref", "PopDen", "PopGro", "ScEc_PC1",
                                            "ScEc_PC2", "ScEc_PC3", "ScEc_PC4", "ScEc_PC5") ~ "Demography",
                           Covariate %in% c("PropVal", "AgProf") ~ "Economic",
                           Covariate %in% c("EcolCond") ~ "Environment",
                           Covariate %in% c("LandUse1", "LandUse2", "LandUse3", "LandUse4",
                                            "NatVegReg1", "NatVegReg2", "PlanZone1", "PlanZone2", "PlanZone3", "PlanZone4") ~ "Policy",
                           .default = NA_character_),
    Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Parcel size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire history",
                               Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
                               Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
                               Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
                               Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
                               Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
                               Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
                               Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
                               Covariate  == "NatVegReg1" ~ "NVR regulated land",
                               Covariate  == "NatVegReg2" ~ "NVR regulated land",
                               Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
                               Covariate  == "PlanZone2" ~ "Planning zone (Others)",
                               Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
                               Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
                               Covariate  == "PopDen" ~ "Population density",
                               Covariate  == "PopGro" ~ "Population Growth",
                               Covariate  == "Precip" ~ "Rainfall",
                               Covariate  == "PropVal" ~ "Land value",
                               Covariate  == "ScEc_PC1" ~ "SE-PC1: Lower income",
                               Covariate  == "ScEc_PC2" ~ "SE-PC2: Australia born, speak English",
                               Covariate  == "ScEc_PC3" ~ "SE-PC3: 1-parent fam., children <15",
                               Covariate  == "ScEc_PC4" ~ "SE-PC4: Coup. fam. no children <15",
                               Covariate  == "ScEc_PC5" ~ "SE-PC5: Coup. fam. children <15",
                               Covariate  == "slope" ~ "Slope",
                               Covariate  == "Soil_PC1" ~ "Soil PC1 Bulk density & sand",
                               Covariate  == "Soil_PC2" ~ "Soil PC2 Organic carbon & silt",
                               Covariate  == "Soil_PC3" ~ "Soil PC3 Nitrogen & water capacity",
                               Covariate  == "Temp" ~ "Temperature",
                               .default = Covariate))%>%
  filter(!Covariate == "(Intercept)")  

# Summarise data for plotting
Cov_df_sum_all <- Cov_df_all %>%
  dplyr::select(Covariate, kmr, CovTypes, Cov_Ag = Cof_PModel_Ag, Cov_Fo = Cof_PModel_Fo, Cov_In = Cof_PModel_In) %>%
  mutate(Cov_Ag = abs(as.numeric(Cov_Ag)),
         Cov_Fo = abs(as.numeric(Cov_Fo)),
         Cov_In = abs(as.numeric(Cov_In))) %>%
  pivot_longer(cols = c("Cov_Ag", "Cov_Fo", "Cov_In"), names_to = "ClearType", values_to = "Count") %>%
  mutate(ClearType = case_when(ClearType == "Cov_Ag" ~ "Agriculture",
                               ClearType == "Cov_Fo" ~ "Forestry",
                               ClearType == "Cov_In" ~ "Infrastructure"),
         Count = as.integer(Count), Count = if_else(!is.na(Count), Count, 0)) %>%
  group_by(Covariate, ClearType, CovTypes) %>%
  summarise(Count = sum(Count), .groups = "drop_last") %>%  ungroup()

# Define factor levels by total count across all clearing types and alphabetical order
ClearType_lvls <- c("Infrastructure", "Forestry", "Agriculture")
CovType_lvls <- c("Accessibility", "Biophysical", "Demography", "Economic", "Environment", "Policy")
# define factor levels for ClearType
Cov_lvls <- Cov_df_sum_all %>%
  arrange(desc(CovTypes), desc(Covariate)) %>%
  pull(Covariate) %>% unique()

Cov_df_sum_all$Covariate <- factor(Cov_df_sum_all$Covariate, levels = Cov_lvls)
Cov_df_sum_all$ClearType <- factor(Cov_df_sum_all$ClearType, levels = ClearType_lvls)
Cov_df_sum_all$CovTypes <- factor(Cov_df_sum_all$CovTypes, levels = CovType_lvls)
head(Cov_df_sum_all)

## Plot bar graph ----
Cov_df_sum_all_plot <- ggplot(Cov_df_sum_all, aes(x = Count, y = Covariate, fill = ClearType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Covariates", fill = "Deforestation\ndriver", 
       x = "Selected\ntimes")+
  scale_fill_manual(values = c("#009E73", "#E69F00", "grey30"),
                    breaks = c("Agriculture", "Forestry", "Infrastructure"),
                    guide = guide_legend(nrow = 2, byrow = TRUE))+
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)),
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme_pubr() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_text(size = 13), axis.text.x = element_text(size = 13, colour = "black"),
        legend.text = element_text(size = 13), legend.title = element_text(size = 13, face = "bold"),
        axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 30.5), expand = FALSE)
# # ggsave(filename = file.path(OUTPUT_DIR, "figures/Cov_df_sum_all_plot.png"), Cov_df_sum_all_plot, width = 10, height = 8, dpi = 300, bg = "white")

# # Manipulate data for plotting
# Cov_df_all2 <- Cov_df_all %>%
#   uncount(weight = 3) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls),
#          x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c)

# # Create long format for plotting
# Cov_df_long <- rbind(Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
#                      Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
#                      Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df2_long <- rbind(Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x1, y1) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
#                       Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x1, y1) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
#                       Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In, x1, y1) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))



# # Create tile plot
# Tile_plot <- ggplot() +
#   geom_tile(data = Cov_df_long, aes(x = kmr, y = Covariate, fill = Cof_PModel), color = "grey80")+

#   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
#   geom_polygon(data = Cov_df2_long, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel), color = "grey80")+
#   facet_wrap(vars(Model))+
#   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40", "grey80") ,
#                     breaks = c(-1, 1, 0, 99), guide = guide_legend(nrow = 2, byrow = TRUE),
#                     labels = c("Negative\ncoefficient", "Positive\ncoefficient", "Not\nSelected", "Excluded"), na.value = "grey80",
#                     name = "Covariate\nselection and\ncoefficient\ndirection")+
#   # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
#   labs(x = "Koala Modelling Regions") +
#   theme(
#     legend.position = "top",
#     legend.justification = "right",
#     legend.box.just = "right"
#   )+
#   theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 11, face = "bold"),
#         axis.text.y = element_text(size = 11, colour = "black"),
#         axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
#         legend.text = element_text(size = 11), legend.title = element_text(size = 11, face = "bold"))
#   #   inherit.aes = FALSE
#   # )
#   # coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)
# Tile_plot
# ggsave(filename = file.path(OUTPUT_DIR, "figures/Tile_plot_AllCovariates2.png"), plot = Tile_plot, width = 16, height = 17.5, dpi = 300)

## Tile plot ----
Cov_df_all_xy <- Cov_df_all%>% 
  mutate(Cof_PModel_Ag = ifelse(is.na(Cof_PModel_Ag), 99, Cof_PModel_Ag),
         Cof_PModel_Fo = ifelse(is.na(Cof_PModel_Fo), 99, Cof_PModel_Fo),
         Cof_PModel_In = ifelse(is.na(Cof_PModel_In), 99, Cof_PModel_In),
         Cof_NModel_Ag = ifelse(is.na(Cof_NModel_Ag), 99, Cof_NModel_Ag),
         Cof_NModel_Fo = ifelse(is.na(Cof_NModel_Fo), 99, Cof_NModel_Fo),
         Cof_NModel_In = ifelse(is.na(Cof_NModel_In), 99, Cof_NModel_In),
         CovTypes = factor(CovTypes, levels = CovType_lvls),
         Covariate = factor(Covariate, levels = Cov_lvls), 
         kmr = factor(kmr),
         x = as.numeric(kmr),
         y = as.numeric(Covariate))

Cov_df_all_xy2 <- Cov_df_all_xy %>% uncount(weight = 3) %>%
  mutate(x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

# Create long format for plotting
Cov_df_long_xy <- rbind(Cov_df_all_xy %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x, y) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                        Cov_df_all_xy %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x, y) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                        Cov_df_all_xy %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_In, Cof_NModel_In, x, y) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df2_long_xy <- rbind(Cov_df_all_xy2 %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x1, y1) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                         Cov_df_all_xy2 %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x1, y1) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                         Cov_df_all_xy2 %>% dplyr::select(Covariate, CovTypes, kmr, Cof_PModel_In, Cof_NModel_In, x1, y1) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

CovType_label_df <- Cov_df_long_xy %>%
  group_by(CovTypes) %>%
  summarise(x = -22, y = max(y), .groups = "drop_last") %>% ungroup() %>% 
  mutate(Model = "Agriculture")

CovType_label_seg_df <- Cov_df_long_xy %>%
  group_by(CovTypes) %>%
  summarise(ymax = max(y) , ymin = min(y), .groups = "drop_last") %>% ungroup() %>%
  uncount(weight = 3) %>% 
  mutate(Model = "Agriculture", REP = rep(1:3, nrow(.)/3),
         y = if_else(REP == 1, ymin - 0.45, if_else(REP == 2, ymax + 0.45, ymin - 0.45)),
         yend = if_else(REP == 1, ymin - 0.45, ymax + 0.45),
         x = -22, xend = if_else(REP != 3, -7.5, x))

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_long_xy, aes(x = x, y = y, fill = Cof_PModel), color = "grey80")+

  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_long_xy, aes(x = x1, y = y1, group = interaction(Covariate, kmr, CovTypes), fill = Cof_NModel), color = "grey80")+
  geom_text(data = CovType_label_df, aes(x = x, y = y, label = CovTypes), fontface = "bold", hjust = 0, size = 4.75,  inherit.aes = FALSE)+
  geom_segment(data = CovType_label_seg_df, aes(x = x, xend = xend, y = y, yend = yend), color = "grey10", inherit.aes = FALSE)+
  facet_wrap(vars(Model))+
  # facet_grid(rows = vars(CovTypes), cols = vars(Model), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = 1:length(levels(Cov_df_all_xy$kmr)), labels = levels(Cov_df_all_xy$kmr)) +
  scale_y_continuous(breaks = 1:length(levels(Cov_df_all_xy$Covariate)), labels = levels(Cov_df_all_xy$Covariate)) +
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40", "grey80") ,
                    breaks = c(-1, 1, 0, 99), guide = guide_legend(nrow = 2, byrow = TRUE),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected", "Excluded"), na.value = "grey80",
                    name = "Covariate selection and\ncoefficient direction")+
  # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  labs(x = "Koala Modelling Regions") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )+
  coord_cartesian(xlim = c(0.5, 9), ylim = c(1.7, 29.2), clip = "off")+
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # strip.placement = "outside", strip.switch.pad.grid = unit(-5, "lines"),
        # strip.background.y = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        # strip.text.y.left = element_text(size = 11, face = "bold", angle = 0, vjust = 0.95, hjust = 0),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(0, 0, 0, 20))
  
LegendKeyTxt_DF <- data.frame(x =     c(-0.075,                       0.075), 
                              y =     c(-0.051,                        0.051), 
                           label = c("Zero-inflation\ncomponent", "Amount\ncomponent"))
LegendKeyPol_DF <- data.frame(x = c(-0.05, -0.05, 0.045, -0.045, 0.05, 0.05), 
                             y = c(-0.045, 0.05, 0.05, -0.05, -0.05, 0.045),
                           group =c(1,1,1,             2,2,2))
LegendKeyArrow <- data.frame(x = c(-0.15, 0.15), 
                             y = c(0.04, -0.04), 
                             xend = c(-0.025, 0.025), 
                             yend = c(0.025, -0.025),
                             curve = c(0.5, -0.5), 
                             group = c(1,2))
LegendKey <- ggplot()+
  geom_polygon(data = LegendKeyPol_DF, aes(x = x, y = y, group = group), fill = "white", color = "black", linewidth = 0.5)+
  geom_text(data = LegendKeyTxt_DF, aes(x = x, y = y, label = label), hjust = c(1, 0), size = 4.1)+
  geom_curve(data = LegendKeyArrow, aes(x = x, y = y, xend = xend, yend = yend), curvature = -0.3, 
                arrow = arrow(length = unit(0.05, "npc"), type = "closed"), linewidth = 0.3)+
  coord_fixed(xlim = c(-0.85, 0.85), ylim = c(-0.18, 0.18), expand = FALSE, clip = "off")+
  # theme_bw()
  theme_void() +
  theme(plot.background  = element_rect(fill = NA, colour = NA), 
        panel.background = element_rect(fill = NA, colour = NA))
LegendKey

leg_left <- cowplot::get_legend(
  Tile_plot +
    theme(legend.title.position = "top", legend.direction = "horizontal",
          legend.position = "top", legend.justification = "left", legend.box.just = "left",
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"),)
)
leg_right <- cowplot::get_legend(
  Cov_df_sum_all_plot +
    theme(legend.title.position = "top", legend.direction = "horizontal",
          legend.position = "top", legend.justification = "left", legend.box.just = "left",
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
)

legend_row <- ggdraw() +
  draw_plot(LegendKey, x = -0.11, y = 0, width = 0.5, height = 1, hjust = 0, vjust = 0) +
  draw_grob(leg_left,  x = 0.27, y = 0, width = 1, height = 1, hjust = 0, vjust = 0) +
  draw_grob(leg_right, x = 0.71, y = 0, width = 1, height = 1, hjust = 0, vjust = 0)

tile_no_leg <- Tile_plot + theme(legend.position = "none")
bar_no_leg  <- Cov_df_sum_all_plot +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y  = element_blank())

main_row <- plot_grid(tile_no_leg, bar_no_leg, nrow = 1, rel_widths = c(6, 1), align = "h", axis = "tb")

final_plot <- ggdraw()+
  draw_plot(legend_row, x = 0, y = 0.85, width = 1, height = 0.15, hjust = 0, vjust = 0) +
  draw_plot(main_row, x = -0.02, y = 0, width = 1.02, height = 0.85, hjust = 0, vjust = 0)
ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_A42.png"), plot = final_plot, width = 8.27, height = 10, units = "in", dpi = 300, bg = "white")

#########################################################################################
# Export data for reporting
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")) %>%
  st_drop_geometry() %>% select(-GRIDCODE) %>% 
  mutate(KMR = c("NC", "CC", "SC", "NT", "NS", "CST", "DRP", "FW", "R"),
         KMRname = case_when(KMRname == "Central and Southern Tablelands" ~ "Central & Southern Tablelands", .default = KMRname))
KMR_LU <- setNames(KMR_shp$KMRname, KMR_shp$KMR)

Cov_CI <- qread(file.path(OUTPUT_DIR, "data/Cov_CI.qs")) %>% 
    mutate(CovTypes = case_when(Covariate %in% c("DistRoad", "DistCity", "Area") ~ "Accessibility",
                              Covariate %in% c("Drought1", "Elevation", "Fire1",
                                                "ForType", "Precip", "slope",
                                                "Soil_PC1", "Soil_PC2", "Soil_PC3", "Temp") ~ "Biophysical",
                              Covariate %in% c("PolPref", "PopDen", "PopGro", "ScEc_PC1",
                                                "ScEc_PC2", "ScEc_PC3", "ScEc_PC4", "ScEc_PC5") ~ "Demography",
                              Covariate %in% c("PropVal", "AgProf") ~ "Economic",
                              Covariate %in% c("EcolCond") ~ "Environment",
                              Covariate %in% c("LandUse1", "LandUse2", "LandUse3", "LandUse4",
                                                "NatVegReg1", "NatVegReg2", "PlanZone1", "PlanZone2", "PlanZone3", "PlanZone4") ~ "Policy",
                              .default = NA_character_),
    Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Parcel size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire history",
                               Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
                               Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
                               Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
                               Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
                               Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
                               Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
                               Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
                               Covariate  == "NatVegReg1" ~ "NVR regulated land",
                               Covariate  == "NatVegReg2" ~ "NVR regulated land",
                               Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
                               Covariate  == "PlanZone2" ~ "Planning zone (Others)",
                               Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
                               Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
                               Covariate  == "PopDen" ~ "Population density",
                               Covariate  == "PopGro" ~ "Population Growth",
                               Covariate  == "Precip" ~ "Rainfall",
                               Covariate  == "PropVal" ~ "Land value",
                               Covariate  == "ScEc_PC1" ~ "SE-PC1: Lower income",
                               Covariate  == "ScEc_PC2" ~ "SE-PC2: Australia born, speak English",
                               Covariate  == "ScEc_PC3" ~ "SE-PC3: 1-parent fam., children <15",
                               Covariate  == "ScEc_PC4" ~ "SE-PC4: Coup. fam. no children <15",
                               Covariate  == "ScEc_PC5" ~ "SE-PC5: Coup. fam. children <15",
                               Covariate  == "slope" ~ "Slope",
                               Covariate  == "Soil_PC1" ~ "Soil PC1 Bulk density & sand",
                               Covariate  == "Soil_PC2" ~ "Soil PC2 Organic carbon & silt",
                               Covariate  == "Soil_PC3" ~ "Soil PC3 Nitrogen & water capacity",
                               Covariate  == "Temp" ~ "Temperature",
                               .default = Covariate)) %>% 
    rename( "2.5%" = LowerCI, "97.5%" = UpperCI, region = KMR) %>% 
    mutate(region = KMR_LU[region], 
           Model = case_when(Model == "PModel" ~ "Zero-inflation component", Model == "NModel" ~ "Amount component", .default = Model)) %>% 
    select("Covariate types" = CovTypes, "Covariate" = Covariate, "Model component" = Model, "Modeling Region" = region, 
           "Coefficient estimate" = mean, "Standard deviation" = sd, "2.5%" = `2.5%`, "97.5%" = `97.5%`) %>%
    arrange(`Covariate types`, Covariate, `Model component`, `Modeling Region`) %>%
    drop_na()

write_csv(Cov_CI, file.path("Cov_CI_All.csv"))
##############################################################################################

# # model_key_plot <- ggplot() +
# #   theme_void() +
# #   theme(
# #     plot.background  = element_rect(fill = NA, colour = NA),
# #     panel.background = element_rect(fill = NA, colour = NA)
# #   ) +
# #   coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "on") +
# #   annotation_custom(model_key_grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
# legend_row <-
#   LegendKey +
#   wrap_elements(full = leg_left) +
#   wrap_elements(full = leg_right)+
#   plot_layout(nrow = 1, widths = c(2, 2, 2))
# legend_row
# # leg_both <- plot_grid(leg_left, leg_right, nrow = 1, align = "h")
# # leg_both_wrap <- patchwork::wrap_elements(full = leg_both)
# # leg_both_wrap
# ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_Legend.png"), plot = legend_row, width = 24, height = 4.5, units = "cm", dpi = 300)

# ## Combine tile plot and bar plot with inset legends
# # remove internal legends
# tile_no_leg <- Tile_plot + theme(legend.position = "none")
# bar_no_leg  <- Cov_df_sum_all_plot +
#   theme(
#     legend.position = "none",
#     axis.title.y = element_blank(),
#     axis.text.y  = element_blank()
#   )

# main_row <- tile_no_leg | bar_no_leg

# TileBar_plot2 <- legend_row / main_row +
#   plot_layout(
#     widths  = c(5, 1),  # barplot = 1/6
#     heights = c(1, 7)   # legend row = 1/8
#   )

# ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_AllCovariates2_v2.png"), plot = TileBar_plot2, width = 16, height = 24.5, units = "cm", dpi = 600)


# TileBar_plot2 <- (CovType_plot|(Tile_plot+theme(legend.position = "none")
#   )|(Cov_df_sum_all_plot +
#      theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
#      theme(legend.position = "none")))+
#   plot_layout(widths = c(1,5,1))
# TileBar_plot2
# ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_AllCovariates2_v2.png"), plot = TileBar_plot2, width = 16, height = 24.5, units = "cm", dpi = 600)


# TileBar_plot <- (((Tile_plot+theme(legend.position = "none")) +
#     patchwork::inset_element(
#     model_key_grob,
#     left = 0.09, right = 0.14,
#     bottom = 0.92, top = 0.98,
#     # left = 0.14, right = 0.19,
#     # bottom = 0.605, top = 0.655,
#     align_to = "full") +
#     patchwork::inset_element(
#     leg_both_wrap,
#     left = 0, right = 0.2,
#     bottom = 0.6, top = 0.9,
#     align_to = "full"      # coordinates relative to the full patchwork
#   ))|
#   (Cov_df_sum_all_plot +
#     theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
#     theme(legend.position = "none"))) +
#   plot_layout(widths = c(5,1))

# TileBar_plot

# ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_AllCovariates2.png"), plot = TileBar_plot, width = 11, height = 14, dpi = 600)

# (Tile_plot|Cov_df_sum_all_plot)+plot_layout(widths = c(5,1))

# ###################################################
# ## Barplot to show positive and negatively associated covariates. (Not used)
# Cov_df_sum2 <- Cov_df %>%
#   mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
#                                Covariate  == "Area" ~ "Parcel size",
#                                Covariate  == "DistCity" ~ "Distance to urban center",
#                                Covariate  == "DistRoad" ~ "Distance to road",
#                                Covariate  == "Drought1" ~ "Drought",
#                                Covariate  == "EcolCond" ~ "Ecological condition",
#                                Covariate  == "Fire1" ~ "Fire history",
#                                Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
#                                Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
#                                Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
#                                Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
#                                Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
#                                Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
#                                Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
#                                Covariate  == "NatVegReg1" ~ "NVR regulated land",
#                                Covariate  == "NatVegReg2" ~ "NVR regulated land",
#                                Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
#                                Covariate  == "PlanZone2" ~ "Planning zone (Others)",
#                                Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
#                                Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
#                                Covariate  == "PopDen" ~ "Population density",
#                                Covariate  == "PopGro" ~ "Population Growth",
#                                Covariate  == "Precip" ~ "Rainfall",
#                                Covariate  == "PropVal" ~ "Land value",
#                                Covariate  == "ScEc_PC1" ~ "Socio-Economic PC1",
#                                Covariate  == "ScEc_PC2" ~ "Socio-Economic PC2",
#                                Covariate  == "ScEc_PC3" ~ "Socio-Economic PC3",
#                                Covariate  == "ScEc_PC4" ~ "Socio-Economic PC4",
#                                Covariate  == "ScEc_PC5" ~ "Socio-Economic PC5",
#                                Covariate  == "slope" ~ "Slope",
#                                Covariate  == "Soil_PC1" ~ "Soil PC1",
#                                Covariate  == "Soil_PC2" ~ "Soil PC2",
#                                Covariate  == "Soil_PC3" ~ "Soil PC3",
#                                Covariate  == "Temp" ~ "Temperature",
#                                .default = Covariate)) %>%
#   distinct() %>%
#   group_by(Covariate) %>%
#   summarise(Cof_PModel_Ag_pos = sum(as.numeric(Cof_PModel_Ag == 1)),
#             Cof_PModel_Ag_neg = sum(as.numeric(Cof_PModel_Ag == -1))*-1,
#             Cof_NModel_Ag_pos = sum(as.numeric(Cof_NModel_Ag == 1)),
#             Cof_NModel_Ag_neg = sum(as.numeric(Cof_NModel_Ag == -1))*-1,
#             Cof_PModel_Ag_sum = sum(as.numeric(Cof_PModel_Ag == 1)) + sum(as.numeric(Cof_PModel_Ag == -1)),
#             Cof_PModel_Fo_pos = sum(as.numeric(Cof_PModel_Fo == 1)),
#             Cof_PModel_Fo_neg = sum(as.numeric(Cof_PModel_Fo == -1))*-1,
#             Cof_NModel_Fo_pos = sum(as.numeric(Cof_NModel_Fo == 1)),
#             Cof_NModel_Fo_neg = sum(as.numeric(Cof_NModel_Fo == -1))*-1,
#             Cof_PModel_Fo_sum = sum(as.numeric(Cof_PModel_Fo == 1)) + sum(as.numeric(Cof_PModel_Fo == -1)),
#             Cof_PModel_In_pos = sum(as.numeric(Cof_PModel_In == 1)),
#             Cof_PModel_In_neg = sum(as.numeric(Cof_PModel_In == -1))*-1,
#             Cof_NModel_In_pos = sum(as.numeric(Cof_NModel_In == 1)),
#             Cof_NModel_In_neg = sum(as.numeric(Cof_NModel_In == -1))*-1,
#             Cof_PModel_In_sum = sum(as.numeric(Cof_PModel_In == 1)) + sum(as.numeric(Cof_PModel_In == -1)))

# Cov_df_sum2 %>% arrange(desc(Cof_PModel_Fo_sum)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
# Cov_df_sum2 %>% arrange(Cof_PModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
# Cov_df_sum2 %>% arrange(desc(Cof_PModel_Fo_pos)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
# Cov_df_sum2 %>% arrange(Cof_PModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
# Cov_df_sum2 %>% arrange(desc(Cof_NModel_Fo_pos)) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)
# Cov_df_sum2 %>% arrange(Cof_NModel_Fo_neg) %>% select(Covariate, Cof_PModel_Fo_pos, Cof_PModel_Fo_neg, Cof_NModel_Fo_pos, Cof_NModel_Fo_neg)

# Cov_df_sum2 %>% arrange(desc(Cof_PModel_In_sum)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
# Cov_df_sum2 %>% arrange(Cof_PModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
# Cov_df_sum2 %>% arrange(desc(Cof_PModel_In_pos)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
# Cov_df_sum2 %>% arrange(Cof_PModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
# Cov_df_sum2 %>% arrange(desc(Cof_NModel_In_pos)) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)
# Cov_df_sum2 %>% arrange(Cof_NModel_In_neg) %>% select(Covariate, Cof_PModel_In_pos, Cof_PModel_In_neg, Cof_NModel_In_pos, Cof_NModel_In_neg)

# # Cov_df_sum_Ag2 <- Cov_df_sum2 %>%
# #   mutate(Cof = abs(Cof_PModel_Ag_pos) + abs(Cof_PModel_Ag_neg)) %>%
# #   arrange(desc(Cof)) %>%
# #   dplyr::select(Covariate,  Cof_PModel_Ag_pos,  Cof_PModel_Ag_neg,  Cof_NModel_Ag_pos,  Cof_NModel_Ag_neg) %>%
# #   drop_na()

# Cov_df_sum_Ag2 <- rbind(Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_pos) %>% mutate(Model = "PModel"),
#                         Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_neg) %>% mutate(Model = "PModel"),
#                         Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_pos) %>% mutate(Model = "NModel"),
#                         Cov_df_sum_Ag2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_neg) %>% mutate(Model = "NModel"))

# Cov_df_sum_Ag2 <- rbind(Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_pos, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "PModel"),
#                         Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_PModel_Ag_neg, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "PModel"),
#                         Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_pos, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "NModel"),
#                         Cov_df_sum2 %>% dplyr::select(Covariate, Cof_Ag = Cof_NModel_Ag_neg, Cof_Ag_sum = Cof_PModel_Ag_sum) %>% mutate(Model = "NModel")) %>%
#   drop_na()

# Cov_sum_Ag2_plot <- ggplot(data = Cov_df_sum_Ag2, aes(x = reorder(Covariate, -Cof_Ag_sum), y = Cof_Ag, fill = Model)) +
#   geom_bar(stat = "identity", position = "dodge")+
#   labs(title = "Agriculture",
#        y = "Number of times Covariate selected") +
#   scale_y_continuous(breaks = seq(-8, 8, by = 4),labels = abs(seq(-8, 8, by = 4)), guide = guide_axis(minor.ticks = TRUE))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.x=element_blank(),
#         plot.title = element_text(hjust=0.9 , vjust = -10))
# ggsave(filename = file.path(OUTPUT_DIR , "figures/Cov_sum_Ag2_plot.png"), Cov_sum_Ag2_plot, width = 11, height = 11, dpi = 300, bg = "white")

# Cov_sum_Ag_plot <- ggplot(data = Cov_df_sum_Ag %>%
#                             dplyr::select(Covariate, Cof_PModel_Ag_pos, Cof_PModel_Ag_neg, Cof_NModel_Ag_pos,Cof_PModel_Ag_neg) %>%
#                             drop_na(),
#                           aes(x = Covariate), y = Cov_Ag) +
#   geom_bar(stat = "identity") +
#   labs(title = "Agriculture",
#        y = "Number of times Covariate selected") +
#   scale_y_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), guide = guide_axis(minor.ticks = TRUE))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.x=element_blank(),
#         plot.title = element_text(hjust=0.9 , vjust = -10))
# Cov_sum_Ag_plot
###################################################

# # Combined tile and bar plot for each clear type----
# ## Prepare data ----
# Cov_df <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs")) %>%
#   mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
#                                Covariate  == "Area" ~ "Parcel size",
#                                Covariate  == "DistCity" ~ "Distance to urban center",
#                                Covariate  == "DistRoad" ~ "Distance to road",
#                                Covariate  == "Drought1" ~ "Drought",
#                                Covariate  == "EcolCond" ~ "Ecological condition",
#                                Covariate  == "Fire1" ~ "Fire history",
#                                Covariate  == "LandTen2" ~ "Land tenure (Leasehold)",
#                                Covariate  == "LandTen3" ~ "Land tenure (Crown purposes)",
#                                Covariate  == "LandTen4" ~ "Land tenure (Other crown land)",
#                                Covariate  == "LandUse1" ~ "Land use (Production-Natural Env)",
#                                Covariate  == "LandUse2" ~ "Land use (Production-Dryland)",
#                                Covariate  == "LandUse3" ~ "Land use (Production-Irrigated)",
#                                Covariate  == "LandUse4" ~ "Land use  (Intensive Uses)",
#                                Covariate  == "NatVegReg1" ~ "NVR regulated land",
#                                Covariate  == "NatVegReg2" ~ "NVR regulated land",
#                                Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
#                                Covariate  == "PlanZone2" ~ "Planning zone (Others)",
#                                Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
#                                Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
#                                Covariate  == "PopDen" ~ "Population density",
#                                Covariate  == "PopGro" ~ "Population Growth",
#                                Covariate  == "Precip" ~ "Rainfall",
#                                Covariate  == "PropVal" ~ "Land value",
#                                Covariate  == "ScEc_PC1" ~ "Socio-Economic PC1 (Lower income)",
#                                Covariate  == "ScEc_PC2" ~ "Socio-Economic PC2 (% Australia born, Eng. speaking)",
#                                Covariate  == "ScEc_PC3" ~ "Socio-Economic PC3 (% 1-parent fam. with/without children u15)",
#                                Covariate  == "ScEc_PC4" ~ "Socio-Economic PC4 (% coup. fam. no children u15)",
#                                Covariate  == "ScEc_PC5" ~ "Socio-Economic PC5 (% coup. fam. with children under 15)",
#                                Covariate  == "slope" ~ "Slope",
#                                Covariate  == "Soil_PC1" ~ "Soil PC1 (High bulk density, sand content)",
#                                Covariate  == "Soil_PC2" ~ "Soil PC2 (High organic carbon, silt content)",
#                                Covariate  == "Soil_PC3" ~ "Soil PC3 (High total nitrogen, avail. water capacity)",
#                                Covariate  == "Temp" ~ "Temperature",
#                                .default = Covariate))

# triangle_plot <- ggdraw() +
#   draw_grob(model_key_grob) +
#   theme(plot.margin = margin(0, 0, 0, 0))
# theme_set()

# ## Agriculture ----
# Cov_df_Ag <- Cov_df %>%
#   dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>%
#   drop_na() %>%
#   filter(Covariate != "(Intercept)")

# Cov_df_Ag_sum <- Cov_df_Ag %>%
#   group_by(Covariate) %>%
#   summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Ag))), .groups = "drop_last")

# Cov_lvls <- Cov_df_Ag_sum %>%
#   arrange(Cof_total, desc(Covariate)) %>%
#   pull(Covariate)

# Cov_df_Ag_sum <- Cov_df_Ag_sum %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df_Ag <- Cov_df_Ag %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df2_Ag <- Cov_df_Ag %>%
#   uncount(weight = 3) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls),
#          x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c)

# Tile_plot <- ggplot() +
#   geom_tile(data = Cov_df_Ag, aes(x = kmr, y = Covariate, fill = Cof_PModel_Ag), color = "grey80")+

#   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
#   geom_polygon(data = Cov_df2_Ag, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Ag), color = "grey80")+
#   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
#                     breaks = c(-1, 1, 0),
#                     labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
#                     name = "Covariate selection and coefficient direction")+
#   labs(x = "Koala Modelling Regions") +
#   theme(
#     legend.position = "top",
#     legend.justification = "right",
#     legend.box.just = "right"
#   )+
#   theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
#   coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

# Bar_plot <- ggplot(data = Cov_df_Ag_sum, aes(x = Cof_total, y = Covariate)) +
#   geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
#   labs(x = "Number of times covariate selected") +
#   scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)),
#                      minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
#   theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
#         axis.text  = element_text(size = 11.5, colour = "black")) +
#   coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

# TileBar_plot_Ag <- Tile_plot+ inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
#  Bar_plot + plot_layout(width = c(5,3.5))

# ggsave(file.path(OUTPUT_DIR, "figures/TileBar_plot_Ag.png"), TileBar_plot_Ag, width = 11, height = 8, dpi = 300, bg = "white")

# ## Forestry ----
# Cov_df_Fo <- Cov_df %>%
#   dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>%
#   drop_na() %>%
#   filter(Covariate != "(Intercept)")

# Cov_df_Fo_sum <- Cov_df_Fo %>%
#   group_by(Covariate) %>%
#   summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Fo))))

# Cov_lvls <- Cov_df_Fo_sum %>%
#   arrange(Cof_total, desc(Covariate)) %>%
#   pull(Covariate)

# Cov_df_Fo_sum <- Cov_df_Fo_sum %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df_Fo <- Cov_df_Fo %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df2_Fo <- Cov_df_Fo %>%
#   uncount(weight = 3) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls),
#          x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c)

# Tile_plot <- ggplot() +
#   geom_tile(data = Cov_df_Fo, aes(x = kmr, y = Covariate, fill = Cof_PModel_Fo), color = "grey80")+

#   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
#   geom_polygon(data = Cov_df2_Fo, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Fo), color = "grey80")+
#   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
#                     breaks = c(-1, 1, 0),
#                     labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
#                     name = "Covariate selection and coefficient direction")+
#   labs(x = "Koala Modelling Regions") +
#   theme(
#     legend.position = "top",
#     legend.justification = "right",
#     legend.box.just = "right"
#   )+
#   theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
#   coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

# Bar_plot <- ggplot(data = Cov_df_Fo_sum, aes(x = Cof_total, y = Covariate)) +
#   geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
#   labs(x = "Number of times covariate selected") +
#   scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)),
#                      minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
#   theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
#         axis.text  = element_text(size = 11.5, colour = "black")) +
#   coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 25.5), expand = FALSE)
# TileBar_plot_Fo <- Tile_plot + inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
#   Bar_plot + plot_layout(width = c(5,3.5))
# ggsave(filename = file.path(OUTPUT_DIR , "figures/TileBar_plot_Fo.png"), TileBar_plot_Fo, width = 11, height = 8, dpi = 300, bg = "white")

# ## Infrastructure ----
# Cov_df_In <- Cov_df %>%
#   dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>%
#   drop_na() %>%
#   filter(Covariate != "(Intercept)")

# Cov_df_In_sum <- Cov_df_In %>%
#   group_by(Covariate) %>%
#   summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_In))))

# Cov_lvls <- Cov_df_In_sum %>%
#   arrange(Cof_total, desc(Covariate)) %>%
#   pull(Covariate)

# Cov_df_In_sum <- Cov_df_In_sum %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Cov_df_In <- Cov_df_In %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls))


# Cov_df2_In <- Cov_df_In %>%
#   uncount(weight = 3) %>%
#   mutate(Covariate = factor(Covariate, levels = Cov_lvls),
#          x = as.numeric(as.factor(kmr)),
#          x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          y = as.numeric(as.factor(Covariate)),
#          y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
#          x1 = x + x_c, y1 = y + y_c)

# Tile_plot <- ggplot() +
#   geom_tile(data = Cov_df_In, aes(x = kmr, y = Covariate, fill = Cof_PModel_In), color = "grey80")+

#   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
#   geom_polygon(data = Cov_df2_In, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_In), color = "grey80")+
#   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
#                     breaks = c(-1, 1, 0),
#                     labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
#                     name = "Covariate selection and coefficient direction")+
#   labs(x = "Koala Modelling Regions") +
#   theme(
#     legend.position = "top",
#     legend.justification = "right",
#     legend.box.just = "right"
#   )+
#   theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
#   coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 24.5), expand = FALSE)

# Bar_plot <- ggplot(data = Cov_df_In_sum, aes(x = Cof_total, y = Covariate)) +
#   geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
#   labs(x = "Number of times covariate selected") +
#   scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)),
#                      minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
#   theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         panel.background = element_blank(), plot.background = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
#         axis.text  = element_text(size = 11.5, colour = "black")) +
#   coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 24.5), expand = FALSE)
# TileBar_plot_In <- Tile_plot + inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
#   Bar_plot + plot_layout(width = c(5,3.5))
# ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_In.png"), TileBar_plot_In, width = 11, height = 8, dpi = 300, bg = "white")



###########################################################################################################################################
# Plot map----
## Map for study area (KMR) ----
# Load SUs prediction shapefile
gpkg_path <- file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg")
SUs_Pred_SF <- st_read(dsn = gpkg_path, layer = "SUs_Predictions") %>%
  drop_na(Pred_All)
head(st_drop_geometry(SUs_Pred_SF))
# Load ABS urban areas shapefile
ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>%
  st_transform(st_crs("EPSG:4283"))

# Get number of SUs in each KMR
nSUs_KMR <- st_drop_geometry(SUs_Pred_SF)  %>%
  group_by(KMR) %>%
  summarise(nSUs = n())

# Load Koala Modelling Regions shapefile, include a buffer zone
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp")) %>%
  st_transform(st_crs("EPSG:4283")) %>% st_make_valid() %>%
  # calculate centroid coordinates for labeling
  mutate(x = st_coordinates(st_centroid(.))[,1], y = st_coordinates(st_centroid(.))[,2],
         # manipulate labels
         KMR = c("(NC)", "(CC)", "(SC)", "(NT)", "(NS)", "(CST)", "(DRP)", "(FW)", "(R)"),
         KMRname = case_when(KMRname == "Central and Southern Tablelands" ~ "Central & Southern Tablelands", .default = KMRname),
         KMRname2 = str_wrap(paste(KMRname, KMR, sep = " "), width = 10),
         # Modify coordinates for NS and NT for better visualisation
         y = if_else(KMR == "(NS)", y - .75, y),
         y = if_else(KMR == "(NT)", y + 0.5, y),
         y = if_else(KMR == "(CST)", y+ 0.61, y),
         x = if_else(KMR == "(CST)", x+ 0.3, x),
         x = if_else(KMR == "(NC)", x + .6, x),
         x = if_else(KMR == "(CC)", x + .5, x)) %>%
  # Add in number of SUs in each KMR
  left_join(., nSUs_KMR, by = c("KMRname" = "KMR")) %>%
  mutate(KMRname3 = paste0(KMRname2, "\n(", nSUs, ")") )

KMR_shp_dsvl <- st_union(KMR_shp)

# Load state boundary shapefile
STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>%
  st_transform(st_crs("EPSG:4283"))
STE1 <- STE %>%
  filter(!STATE_NAME %in% c("New South Wales", "Other Territories")) %>%
  mutate(x = c(145, 145, 135, NA,NA,NA, 147.2), y =  c(-36.5, -28.5, -44,NA,NA,NA, -36.5))  %>%
  mutate(STATE_NAME = case_when(STATE_NAME == "Australian Capital Territory" ~ "A.C.T.", .default = STATE_NAME))

URB_KMR <- st_intersection(st_centroid(ABS_urb), KMR_shp) %>%
  st_drop_geometry() %>%
  group_by(KMRname) %>%
  slice_max(n=2, order_by = AREASQKM16)

# Select urban areas for labelling
NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Casino", "Port Macquarie", "Sydney", "Bega", "Narrabri", "Dubbo", "Wagga Wagga", "Broken Hill", "Cobar")) %>% 
  mutate(UCL_NAME16 = case_when(
    UCL_NAME16 == "Wagga Wagga" ~ "Wagga\nWagga",
    UCL_NAME16 == "Broken Hill" ~ "Broken\nHill",
    .default = UCL_NAME16 )) %>%
  dplyr::select(UCL_NAME16, geometry) %>%
  distinct(UCL_NAME16, .keep_all = TRUE) %>%
  st_centroid(.) %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2],
         x2 = x, y2 = y,
         x2 = if_else(UCL_NAME16 %in% c("Port Macquarie"), x2 + 1.1, x2),
         x2 = if_else(UCL_NAME16 %in% c("Sydney"), x2 + 0.9, x2),
         x2 = if_else(UCL_NAME16 %in% c("Bega"), x2 + 0.6, x2),
         x2 = if_else(UCL_NAME16 %in% c("Broken\nHill", "Dubbo", "Wagga\nWagga"), x2 + 0.3, x2),
         x2 = if_else(UCL_NAME16 %in% c("Narrabri", "Casino"), x2 - 0.1, x2),
         y2 = if_else(UCL_NAME16 %in% c("Broken\nHill", "Wagga\nWagga"), y2 - 0.4, y2),
         y2 = if_else(UCL_NAME16 %in% c("Narrabri", "Dubbo", "Cobar", "Casino"), y2 - 0.2, y2))

# Create KMR map
KMR_map <- ggplot()+

  # Koala Modelling Regions
  geom_sf(data = KMR_shp, fill = "lightgoldenrod1", color = "lightgoldenrod3", lwd = 0.6, linetype= "solid", alpha = 0.9)+

  # State boundary grey background
  geom_sf(data = STE1, fill = "grey80", color = "black", lwd = 0.5)+
  geom_text(data = STE1, aes(x = x, y = y, label = STATE_NAME), colour = "grey40", size = 7)+
  annotate("segment", x = 147.7, xend = 149, y = -36.5, yend = -35.7, linewidth = 1, colour = "grey40", alpha = 0.9, lineend = "round")+
  geom_sf(data = STE, fill = "transparent", color = "black", lwd = 0.2)+

  # Urban areas points
  geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1.5)+
  # Labels for KMR and urban areas
  geom_shadowtext(data = KMR_shp, aes(x = x, y = y, label = KMRname2), fontface = "bold", nudge_y = 0, size = 6, color = "black", bg.color = "white", bg.r = 0.2)+
  geom_shadowtext(data = NSW_urb_sel_pt, aes(x = x2, y = y2 , label = UCL_NAME16), size = 5.5, color = "grey10", bg.color = "white",  bg.r = 0.05)+ 
  
  # North arrow and scale bar
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", height = unit(3, "cm"), width = unit(3, "cm"), pad_y = unit(1, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering(fill = c("black", "white"), text_size = 16, line_width = 2, line_col = "black", text_col = "black"))+
  ggspatial::annotation_scale(location = "bl", width_hint = 0.4, bar_height = unit(0.5, "cm"), line_width = 2, pad_x = unit(0.5, "cm"), text_cex = 1.2)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = "#C7E6F5"))+
  coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)]+c(-.2, 1.5), ylim = st_bbox(KMR_shp)[c(2,4)]+c(-.1, .1), expand = FALSE)
ggsave(file.path(OUTPUT_DIR, "figures/KMR_map_base.png"), KMR_map, width = 11, height = 11, dpi = 300)

STE_SA_plot <- ggplot()+
  geom_sf(data = KMR_shp_dsvl, fill = "grey20", color = "black", lwd = 0.2)+
  geom_sf(data = STE1, fill = "grey60", color = "grey10", lwd = 0.2)+
  coord_sf(xlim = c(111, 156), ylim = c(-44.5, -9.5), expand = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme_void()+theme(legend.position = "none")+
  theme(plot.background  = element_rect(fill = "#91daff", color = "black", linewidth = 0.7))

KMR_map_comb <- KMR_map +
  inset_element(STE_SA_plot, left = 0.72, bottom = 0.01, right = 0.99, top = 0.35)

# Export KMR map
ggsave(file.path(OUTPUT_DIR, "figures/KMR_map.png"), KMR_map_comb, width = 11, height = 11, dpi = 600, bg = "transparent")


## Deforestation risk maps----
### load base layers and target data  layers
#### Data Layers
SUs_Pred_SF <- st_read(file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg"), layer = "SUs_Predictions") %>%
  drop_na(Pred_All)
Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))
Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))

#### Base Layers
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp"))

ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>%
  st_transform(st_crs(Pred_Ag))

STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>%
  st_transform(st_crs(Pred_Ag))

STE_GDA20 <- st_read("R:/KOALACC22-Q5611/Data/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp") %>%
  filter(STE_CODE21 %in% c(1:8))
NSW_GDA2020 <- STE_GDA20 %>% filter(STE_NAME21 == "New South Wales")

NSW <- STE %>% filter(STATE_NAME == "New South Wales")
NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Casino", "Port Macquarie",
                           "Narrabri", "Dubbo",
                           "Newcastle", "Sydney",
                           "Oberon", "Wagga Wagga",
                           "Griffith", "Deniliquin",
                           "Brewarrina (L)", "Robinvale",
                           "Broken Hill", "Ivanhoe (L)", "Cobar")) %>%
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Brewarrina (L)" ~ "Brewarrina",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                .default = UCL_NAME16),
         UCL_NAME16 = str_wrap(UCL_NAME16, width = 8)) %>%
  dplyr::select(UCL_NAME16, geometry) %>%
  distinct(UCL_NAME16, .keep_all = TRUE) %>%
  st_centroid(.) %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])

# PLOTMAP_risk_with_3Insets
NSW_urb_pt <- ABS_urb %>%
  filter(!str_detect(UCL_NAME16, "emain")) %>%
  st_centroid() %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>%
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Mungindi (NSW Part) (L)" ~ "Mungindi (NSW Part)",
                                UCL_NAME16 == "Lismore" ~ "Lismore",
                                UCL_NAME16 == "Bonalbo (L)" ~ "Bonalbo",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                UCL_NAME16 == "Jilliby (L)" ~ "Jilliby",
                                UCL_NAME16 == "Collarenebri (L)" ~ "Collarenebri",
                                UCL_NAME16 == "Port Macquarie" ~ "Port\nMacquarie",
                                UCL_NAME16 == "Forster - Tuncurry" ~ "Forster-\nTuncurry",
                                UCL_NAME16 == "Blue Mountains" ~ "Blue\nMountains",
                                .default = UCL_NAME16))
############

# PLotting parametres
STE_SA_plot <- ggplot()+
  geom_sf(data = STE_GDA20, fill = "grey70", color = "grey10", lwd = 0.5)+
  geom_sf(data = NSW_GDA2020, fill = "grey20", color = "black", lwd = 0.5)+
  coord_sf(xlim = c(111, 156), ylim = c(-44.5, -9.5), expand = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme_void()+theme(legend.position = "none")+
  theme(plot.background  = element_rect(fill = "#91daff", color = "black", linewidth = 0.7))

# XLIM1 = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000)
# XLIM2 = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 150000)
# TOTAL_X = (XLIM2[2]-XLIM2[1]) + (XLIM1[2]-XLIM1[1])
# XLIM1_RATIO = (XLIM1[2]-XLIM1[1])/TOTAL_X
# XLIM2_RATIO = (XLIM2[2]-XLIM2[1])/TOTAL_X
# DATA = SUs_Pred_SF %>% filter(KMR == "Far West")
# FILL = Pred_All
# LEGEND_Title = "Deforestation\nrisk"
# COLOUR = hcl.colors(8, palette = "Reds 3" ,rev = TRUE)
# LEG_POS = c(0.895, 0.32)
# XLIM = XLIM2
# YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000)
# FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Test1.png")
# WH_Ratio = (XLIM[2]-XLIM[1])/(YLIM[2]-YLIM[1])
# PNG_width = XLIM1_RATIO * 16
# PNG_height = PNG_width / WH_Ratio
# PNG_dpi = 300

# SUs_Pred_range <- st_drop_geometry(SUs_Pred_SF) %>%
#   summarise(
#     min_Pred_Ag = min(Pred_Ag, na.rm = TRUE), max_Pred_Ag = max(Pred_Ag, na.rm = TRUE),
#     min_Pred_Fo = min(Pred_Fo, na.rm = TRUE), max_Pred_Fo = max(Pred_Fo, na.rm = TRUE),
#     min_Pred_In = min(Pred_In, na.rm = TRUE), max_Pred_In = max(Pred_In, na.rm = TRUE),
#     min_Pred_All = min(Pred_All, na.rm = TRUE), max_Pred_All = max(Pred_All, na.rm = TRUE),
#     min_KhabRisk_Ag = min(KhabRisk_Ag, na.rm = TRUE), max_KhabRisk_Ag = max(KhabRisk_Ag, na.rm = TRUE),
#     min_KhabRisk_Fo = min(KhabRisk_Fo, na.rm = TRUE), max_KhabRisk_Fo = max(KhabRisk_Fo, na.rm = TRUE),
#     min_KhabRisk_In = min(KhabRisk_In, na.rm = TRUE), max_KhabRisk_In = max(KhabRisk_In, na.rm = TRUE),
#     min_KhabRisk_All = min(KhabRisk_All, na.rm = TRUE), max_KhabRisk_All = max(KhabRisk_All, na.rm = TRUE),
#     QTL_25_Pred_Ag = quantile(Pred_Ag, 0.25, na.rm = TRUE), QTL_75_Pred_Ag = quantile(Pred_Ag, 0.75, na.rm = TRUE),
#     QTL_25_Pred_Fo = quantile(Pred_Fo, 0.25, na.rm = TRUE), QTL_75_Pred_Fo = quantile(Pred_Fo, 0.75, na.rm = TRUE),
#     QTL_25_Pred_In = quantile(Pred_In, 0.25, na.rm = TRUE), QTL_75_Pred_In = quantile(Pred_In, 0.75, na.rm = TRUE),
#     QTL_25_Pred_All = quantile(Pred_All, 0.25, na.rm = TRUE), QTL_75_Pred_All = quantile(Pred_All, 0.75, na.rm = TRUE),
#     QTL_25_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.25, na.rm = TRUE), QTL_75_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.75, na.rm = TRUE),
#     QTL_25_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.25, na.rm = TRUE), QTL_75_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.75, na.rm = TRUE),
#     QTL_25_KhabRisk_In = quantile(KhabRisk_In, 0.25, na.rm = TRUE), QTL_75_KhabRisk_In = quantile(KhabRisk_In, 0.75, na.rm = TRUE),
#     QTL_25_KhabRisk_All = quantile(KhabRisk_All, 0.25, na.rm = TRUE), QTL_75_KhabRisk_All = quantile(KhabRisk_All, 0.75, na.rm = TRUE),
#     IQR_Pred_Ag = QTL_75_Pred_Ag - QTL_25_Pred_Ag, IQR_Pred_Fo = QTL_75_Pred_Fo - QTL_25_Pred_Fo, IQR_Pred_In = QTL_75_Pred_In - QTL_25_Pred_In, IQR_Pred_All = QTL_75_Pred_All - QTL_25_Pred_All,
#     IQR_KhabRisk_Ag = QTL_75_KhabRisk_Ag - QTL_25_KhabRisk_Ag, IQR_KhabRisk_Fo = QTL_75_KhabRisk_Fo - QTL_25_KhabRisk_Fo, IQR_KhabRisk_In = QTL_75_KhabRisk_In - QTL_25_KhabRisk_In, IQR_KhabRisk_All = QTL_75_KhabRisk_All - QTL_25_KhabRisk_All,
#     Lower_Lim_Pred_Ag = QTL_25_Pred_Ag - 1.5 * IQR_Pred_Ag, Upper_Lim_Pred_Ag = QTL_75_Pred_Ag + 1.5 * IQR_Pred_Ag,
#     Lower_Lim_Pred_Fo = QTL_25_Pred_Fo - 1.5 * IQR_Pred_Fo, Upper_Lim_Pred_Fo = QTL_75_Pred_Fo + 1.5 * IQR_Pred_Fo,
#     Lower_Lim_Pred_In = QTL_25_Pred_In - 1.5 * IQR_Pred_In, Upper_Lim_Pred_In = QTL_75_Pred_In + 1.5 * IQR_Pred_In,
#     Lower_Lim_Pred_All = QTL_25_Pred_All - 1.5 * IQR_Pred_All, Upper_Lim_Pred_All = QTL_75_Pred_All + 1.5 * IQR_Pred_All,
#     Lower_Lim_KhabRisk_Ag = QTL_25_KhabRisk_Ag - 1.5 * IQR_KhabRisk_Ag, Upper_Lim_KhabRisk_Ag = QTL_75_KhabRisk_Ag + 1.5 * IQR_KhabRisk_Ag,
#     Lower_Lim_KhabRisk_Fo = QTL_25_KhabRisk_Fo - 1.5 * IQR_KhabRisk_Fo, Upper_Lim_KhabRisk_Fo = QTL_75_KhabRisk_Fo + 1.5 * IQR_KhabRisk_Fo,
#     Lower_Lim_KhabRisk_In = QTL_25_KhabRisk_In - 1.5 * IQR_KhabRisk_In, Upper_Lim_KhabRisk_In = QTL_75_KhabRisk_In + 1.5 * IQR_KhabRisk_In,
#     Lower_Lim_KhabRisk_All = QTL_25_KhabRisk_All - 1.5 * IQR_KhabRisk_All, Upper_Lim_KhabRisk_All = QTL_75_KhabRisk_All + 1.5 * IQR_KhabRisk_All,
#     QTL_90_Pred_Ag = quantile(Pred_Ag, 0.90, na.rm = TRUE), QTL_95_Pred_Ag = quantile(Pred_Ag, 0.95, na.rm = TRUE), QTL_99_Pred_Ag = quantile(Pred_Ag, 0.99, na.rm = TRUE), QTL_995_Pred_Ag = quantile(Pred_Ag, 0.995, na.rm = TRUE),
#     QTL_90_Pred_Fo = quantile(Pred_Fo, 0.90, na.rm = TRUE), QTL_95_Pred_Fo = quantile(Pred_Fo, 0.95, na.rm = TRUE), QTL_99_Pred_Fo = quantile(Pred_Fo, 0.99, na.rm = TRUE), QTL_995_Pred_Fo = quantile(Pred_Fo, 0.995, na.rm = TRUE),
#     QTL_90_Pred_In = quantile(Pred_In, 0.90, na.rm = TRUE), QTL_95_Pred_In = quantile(Pred_In, 0.95, na.rm = TRUE), QTL_99_Pred_In = quantile(Pred_In, 0.99, na.rm = TRUE), QTL_995_Pred_In = quantile(Pred_In, 0.995, na.rm = TRUE),
#     QTL_90_Pred_All = quantile(Pred_All, 0.90, na.rm = TRUE), QTL_95_Pred_All = quantile(Pred_All, 0.95, na.rm = TRUE), QTL_99_Pred_All = quantile(Pred_All, 0.99, na.rm = TRUE), QTL_995_Pred_All = quantile(Pred_All, 0.995, na.rm = TRUE),
#     QTL_90_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.90, na.rm = TRUE), QTL_95_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.95, na.rm = TRUE), QTL_99_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.99, na.rm = TRUE), QTL_995_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.995, na.rm = TRUE),
#     QTL_90_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.90, na.rm = TRUE), QTL_95_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.95, na.rm = TRUE), QTL_99_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.99, na.rm = TRUE), QTL_995_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.995, na.rm = TRUE),
#     QTL_90_KhabRisk_In = quantile(KhabRisk_In, 0.90, na.rm = TRUE), QTL_95_KhabRisk_In = quantile(KhabRisk_In, 0.95, na.rm = TRUE), QTL_99_KhabRisk_In = quantile(KhabRisk_In, 0.99, na.rm = TRUE), QTL_995_KhabRisk_In = quantile(KhabRisk_In, 0.995, na.rm = TRUE),
#     QTL_90_KhabRisk_All = quantile(KhabRisk_All, 0.90, na.rm = TRUE), QTL_95_KhabRisk_All = quantile(KhabRisk_All, 0.95, na.rm = TRUE), QTL_99_KhabRisk_All = quantile(KhabRisk_All, 0.99, na.rm = TRUE), QTL_995_KhabRisk_All = quantile(KhabRisk_All, 0.995, na.rm = TRUE),
#     QTL_98_Pred_Ag = quantile(Pred_Ag, 0.98, na.rm = TRUE), QTL_98_Pred_Fo = quantile(Pred_Fo, 0.98, na.rm = TRUE), QTL_98_Pred_In = quantile(Pred_In, 0.98, na.rm = TRUE), QTL_98_Pred_All = quantile(Pred_All, 0.98, na.rm = TRUE),
#     QTL_98_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.98, na.rm = TRUE), QTL_98_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.98, na.rm = TRUE), QTL_98_KhabRisk_In = quantile(KhabRisk_In, 0.98, na.rm = TRUE), QTL_98_KhabRisk_All = quantile(KhabRisk_All, 0.98, na.rm = TRUE),
#     QTL_99_Pred_Ag = quantile(Pred_Ag, 0.99, na.rm = TRUE), QTL_99_Pred_Fo = quantile(Pred_Fo, 0.99, na.rm = TRUE), QTL_99_Pred_In = quantile(Pred_In, 0.99, na.rm = TRUE), QTL_99_Pred_All = quantile(Pred_All, 0.99, na.rm = TRUE),
#     QTL_99_KhabRisk_Ag = quantile(KhabRisk_Ag, 0.99, na.rm = TRUE), QTL_99_KhabRisk_Fo = quantile(KhabRisk_Fo, 0.99, na.rm = TRUE), QTL_99_KhabRisk_In = quantile(KhabRisk_In, 0.99, na.rm = TRUE), QTL_99_KhabRisk_All = quantile(KhabRisk_All, 0.99, na.rm = TRUE)
#   ) %>%
#   select(starts_with("min_"), starts_with("max_"),starts_with("Upper_Lim_"), starts_with("Lower_Lim_"), starts_with("QTL_90_"), starts_with("QTL_95_"), starts_with("QTL_98_"), starts_with("QTL_99_"), starts_with("QTL_995_"))
# SUs_Pred_range

# Lim_Pred_MIN_IQR <- max(min(SUs_Pred_range$Lower_Lim_Pred_Ag, SUs_Pred_range$Lower_Lim_Pred_Fo, SUs_Pred_range$Lower_Lim_Pred_In, SUs_Pred_range$Lower_Lim_Pred_All), 0)
# Lim_Pred_MAX_IQR <- min(max(SUs_Pred_range$Upper_Lim_Pred_Ag, SUs_Pred_range$Upper_Lim_Pred_Fo, SUs_Pred_range$Upper_Lim_Pred_In, SUs_Pred_range$Upper_Lim_Pred_All), 1)
# Lim_KhabRisk_MIN_IQR <- max(min(SUs_Pred_range$Lower_Lim_KhabRisk_Ag, SUs_Pred_range$Lower_Lim_KhabRisk_Fo, SUs_Pred_range$Lower_Lim_KhabRisk_In, SUs_Pred_range$Lower_Lim_KhabRisk_All), 0)
# Lim_KhabRisk_MAX_IQR <- min(max(SUs_Pred_range$Upper_Lim_KhabRisk_Ag, SUs_Pred_range$Upper_Lim_KhabRisk_Fo, SUs_Pred_range$Upper_Lim_KhabRisk_In, SUs_Pred_range$Upper_Lim_KhabRisk_All), 1)

# Lim_Pred_MAX_90Q <- min(max(SUs_Pred_range$QTL_90_Pred_Ag, SUs_Pred_range$QTL_90_Pred_Fo, SUs_Pred_range$QTL_90_Pred_In, SUs_Pred_range$QTL_90_Pred_All), 1)
# Lim_KhabRisk_MAX_90Q <- min(max(SUs_Pred_range$QTL_90_KhabRisk_Ag, SUs_Pred_range$QTL_90_KhabRisk_Fo, SUs_Pred_range$QTL_90_KhabRisk_In, SUs_Pred_range$QTL_90_KhabRisk_All), 1)

# Lim_Pred_MAX_95Q <- min(max(SUs_Pred_range$QTL_95_Pred_Ag, SUs_Pred_range$QTL_95_Pred_Fo, SUs_Pred_range$QTL_95_Pred_In, SUs_Pred_range$QTL_95_Pred_All), 1)
# Lim_KhabRisk_MAX_95Q <- min(max(SUs_Pred_range$QTL_95_KhabRisk_Ag, SUs_Pred_range$QTL_95_KhabRisk_Fo, SUs_Pred_range$QTL_95_KhabRisk_In, SUs_Pred_range$QTL_95_KhabRisk_All), 1)

# Lim_Pred_MAX_98Q <- min(max(SUs_Pred_range$QTL_98_Pred_Ag, SUs_Pred_range$QTL_98_Pred_Fo, SUs_Pred_range$QTL_98_Pred_In, SUs_Pred_range$QTL_98_Pred_All), 1)
# Lim_KhabRisk_MAX_98Q <- min(max(SUs_Pred_range$QTL_98_KhabRisk_Ag, SUs_Pred_range$QTL_98_KhabRisk_Fo, SUs_Pred_range$QTL_98_KhabRisk_In, SUs_Pred_range$QTL_98_KhabRisk_All), 1)

# Lim_Pred_MAX_99Q <- min(max(SUs_Pred_range$QTL_99_Pred_Ag, SUs_Pred_range$QTL_99_Pred_Fo, SUs_Pred_range$QTL_99_Pred_In, SUs_Pred_range$QTL_99_Pred_All), 1)
# Lim_KhabRisk_MAX_99Q <- min(max(SUs_Pred_range$QTL_99_KhabRisk_Ag, SUs_Pred_range$QTL_99_KhabRisk_Fo, SUs_Pred_range$QTL_99_KhabRisk_In, SUs_Pred_range$QTL_99_KhabRisk_All), 1)

# Lim_Pred_MAX_995Q <- min(max(SUs_Pred_range$QTL_995_Pred_Ag, SUs_Pred_range$QTL_995_Pred_Fo, SUs_Pred_range$QTL_995_Pred_In, SUs_Pred_range$QTL_995_Pred_All), 1)
# Lim_KhabRisk_MAX_995Q <- min(max(SUs_Pred_range$QTL_995_KhabRisk_Ag, SUs_Pred_range$QTL_995_KhabRisk_Fo, SUs_Pred_range$QTL_995_KhabRisk_In, SUs_Pred_range$QTL_995_KhabRisk_All), 1)

# Lim_Pred_MAX_100 <- min(max(SUs_Pred_range$max_Pred_Ag, SUs_Pred_range$max_Pred_Fo, SUs_Pred_range$max_Pred_In, SUs_Pred_range$max_Pred_All), 1)
# Lim_KhabRisk_MAX_100 <- min(max(SUs_Pred_range$max_KhabRisk_Ag, SUs_Pred_range$max_KhabRisk_Fo, SUs_Pred_range$max_KhabRisk_In, SUs_Pred_range$max_KhabRisk_All), 1)

# Lim_Pred_MIN <- 0
# Lim_KhabRisk_MIN <- 0



#########################
# Plot 4 combined maps----


# ## 99 percentile limit plots
# Ag_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# Fo_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# In_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In, LEGEND_Title = "Deforestation\nrisk",
#                          COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = c(0,0.1, 0.2, 0.3, 0.44), SCALE_LABELS = c("0", "0.1", "0.2", "0.3", "≥0.44"),
#                          North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# All_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF, FILL = Pred_All, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

# Combine_plot <- ((Ag_risk_plot | In_risk_plot) / (Fo_risk_plot | All_risk_plot) +
#   plot_annotation(tag_levels = 'A') &
#   theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
#         plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
#   inset_element(
#     STE_SA_plot + theme(plot.tag = element_blank()),
#     left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)
# ggsave(file.path(OUTPUT_DIR, "figures/Pred_map4_Lim_99Q.png"), Combine_plot, width = 14.5, height = 11.2, dpi = 400)

# Ag_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# Fo_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# In_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In, LEGEND_Title = "Koala habitat\nloss risk",
#                          COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.2439), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.24"),
#                          North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# All_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF, FILL = KhabRisk_All, LEGEND_Title = NULL,
#                          COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
#                          SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
#                          North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
#                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
# Combine_KoalaHab_plot <- ((Ag_KoalaHab_plot | In_KoalaHab_plot) / (Fo_KoalaHab_plot | All_KoalaHab_plot) +
#   plot_annotation(tag_levels = 'A') &
#   theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
#         plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
#   inset_element(
#     STE_SA_plot + theme(plot.tag = element_blank()),
#     left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)
# ggsave(file.path(OUTPUT_DIR, "figures/KoalaHab_map4_Lim_99Q.png"), Combine_KoalaHab_plot, width = 14.5, height = 11.2, dpi = 400)

## Prepare city labels for the map
## This is needed for PLOTMAP_risk2() function
NSW_urb_sel_pt2 <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Casino", "Port Macquarie", "Narrabri", "Dubbo", "Newcastle", "Sydney", 
                           "Oberon", "Wagga Wagga", "Deniliquin","Brewarrina (L)",
                           "Robinvale", "Broken Hill", "Ivanhoe (L)", "Cobar")) %>%
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Brewarrina (L)" ~ "Brewarrina",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                UCL_NAME16 == "Broken Hill" ~ "Broken\nHill",
                                UCL_NAME16 == "Port Macquarie" ~ "Port\nMacquarie",
                                .default = UCL_NAME16)) %>%
  dplyr::select(UCL_NAME16, geometry) %>% distinct(UCL_NAME16, .keep_all = TRUE) %>%
  st_centroid(.) %>% mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% 
  mutate(x = case_when(UCL_NAME16 == "Casino" ~ x + 135000,
                       UCL_NAME16 == "Port\nMacquarie" ~ x +80000,
                       UCL_NAME16 == "Narrabri" ~ x - 160000,
                       UCL_NAME16 == "Dubbo" ~ x +55000,
                       UCL_NAME16 == "Newcastle" ~ x +115000,
                       UCL_NAME16 == "Sydney" ~ x +110000,
                       UCL_NAME16 == "Oberon" ~ x -35000,
                       UCL_NAME16 == "Wagga Wagga" ~ x -85000,
                       UCL_NAME16 == "Deniliquin" ~ x - 170000,
                       UCL_NAME16 == "Brewarrina" ~ x -155000,
                       UCL_NAME16 == "Robinvale" ~ x - 65000,
                       UCL_NAME16 == "Cobar" ~ x -15000,
                       UCL_NAME16 == "Broken\nHill" ~ x + 35000, 
                       .default = x),
         y = case_when(UCL_NAME16 == "Casino" ~ y - 30000,
                       UCL_NAME16 == "Narrabri" ~ y + 180000,
                       UCL_NAME16 == "Port\nMacquarie" ~ y -10000,
                       UCL_NAME16 == "Newcastle" ~ y - 15000,
                       UCL_NAME16 == "Dubbo" ~ y - 30000,
                       UCL_NAME16 == "Oberon" ~ y + 30000,
                       UCL_NAME16 == "Wagga Wagga" ~ y - 130000,
                       UCL_NAME16 == "Deniliquin" ~ y - 30000,
                       UCL_NAME16 == "Brewarrina" ~ y + 130000,
                       UCL_NAME16 == "Robinvale" ~ y - 35000,
                       UCL_NAME16 == "Ivanhoe" ~ y -30000,
                       UCL_NAME16 == "Cobar" ~ y + 35000,
                       UCL_NAME16 == "Broken\nHill" ~ y + 65000, 
                       .default = y))
                      
# Segments to connect city labels to their points
NSW_urb_sel_pt2_seg_df <- NSW_urb_sel_pt2 %>% 
  mutate(xend = st_coordinates(.)[,1], yend = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% as.data.frame() %>% 
  filter(UCL_NAME16 %in% c("Brewarrina", "Narrabri", "Wagga Wagga", "Deniliquin", "Casino")) %>% 
  mutate(yend = case_when(UCL_NAME16 == "Brewarrina" ~ yend + 10000,
                          UCL_NAME16 == "Narrabri" ~ yend + 10000,
                          UCL_NAME16 == "Deniliquin" ~ yend - 7500,
                          UCL_NAME16 == "Wagga Wagga" ~ yend - 10000,
                          UCL_NAME16 == "Casino" ~ yend - 7500, 
                          .default = yend),
         xend = if_else(UCL_NAME16 == "Casino", xend + 10000, 
                        if_else(UCL_NAME16 == "Deniliquin", xend - 10000, xend + 2500)),
         y = if_else(UCL_NAME16 == "Wagga Wagga", y + 15000, y),
         x = case_when(UCL_NAME16 == "Brewarrina" ~ x + 115000,
                       UCL_NAME16 == "Narrabri" ~ x + 90000,
                       UCL_NAME16 == "Wagga Wagga" ~ x + 120000, 
                       UCL_NAME16 == "Deniliquin" ~ x + 100000,
                       UCL_NAME16 == "Casino" ~ x - 75000,
                       .default = x))

# Identify the 0.25 limit for Pred_All
ecdf(SUs_Pred_SF$Pred_All) %>% plot()
F_hat <- ecdf(SUs_Pred_SF$Pred_All)
F_hat(0.25)

Ag_risk_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Fo_risk_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
In_risk_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In, LEGEND_Title = "Deforestation\nrisk",
                         COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = c(0,0.05, 0.1, 0.15, 0.2, 0.25), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.25"),
                         North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
All_risk_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF, FILL = Pred_All, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Combine_plot <- ((Ag_risk_plot | In_risk_plot) / (Fo_risk_plot | All_risk_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
        plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
  inset_element(
    STE_SA_plot + theme(plot.tag = element_blank()),
    left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)
ggsave(file.path(OUTPUT_DIR, "figures/Pred_map4_Lim_0p25.png"), Combine_plot, width = 14.5, height = 11.2, dpi = 400)

# Identify the 0.25 limit for KhabRisk_All
# ecdf(SUs_Pred_SF$KhabRisk_All) %>% plot()
F_hat_khab <- ecdf(SUs_Pred_SF$KhabRisk_All)
F_hat_khab(0.25)

Ag_KoalaHab_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Fo_KoalaHab_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                          XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
In_KoalaHab_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In, LEGEND_Title = "Koala habitat\nloss risk",
                         COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), 
                         SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.25"),
                         North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
All_KoalaHab_plot <- PLOTMAP_risk2(DATA = SUs_Pred_SF, FILL = KhabRisk_All, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
                         SCALE_LIMIT = c(0, 0.25), SCALE_BREAKS = NULL, SCALE_LABELS = NULL,
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Combine_KoalaHab_plot <- ((Ag_KoalaHab_plot | In_KoalaHab_plot) / (Fo_KoalaHab_plot | All_KoalaHab_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
        plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
  inset_element(
    STE_SA_plot + theme(plot.tag = element_blank()),
    left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)

ggsave(file.path(OUTPUT_DIR, "figures/KoalaHab_map4_Lim_0p25.png"), Combine_KoalaHab_plot, width = 14.5, height = 11.2, dpi = 400)



# # Plot maps for each max limit types-----
# Lim_Pred_MAX <- c(Lim_Pred_MAX_IQR, Lim_Pred_MAX_95Q, Lim_Pred_MAX_99Q, Lim_Pred_MAX_100, 0.25, 0.5)
# Lim_KhabRisk_MAX <- c(Lim_KhabRisk_MAX_IQR, Lim_KhabRisk_MAX_95Q, Lim_KhabRisk_MAX_99Q, Lim_KhabRisk_MAX_100, 0.25, 0.5)

# for(i in 4:6){
#   Ag_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_Pred_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   Fo_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_Pred_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   In_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In, LEGEND_Title = "Deforestation\nrisk",
#                           COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_Pred_MAX[i]),
#                           North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   All_risk_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF, FILL = Pred_All, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_Pred_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

#   Combine_plot <- ((Ag_risk_plot | In_risk_plot) / (Fo_risk_plot | All_risk_plot) +
#     plot_annotation(tag_levels = 'A') &
#     theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
#           plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
#     inset_element(
#       STE_SA_plot + theme(plot.tag = element_blank()),
#       left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)
#   Output_FPATH <- file.path(OUTPUT_DIR, "figures", paste0("Pred_map4_Lim_", c("IQR", "95Q", "99Q", "Max", "25", "50")[i], ".png"))
#   ggsave(Output_FPATH, Combine_plot, width = 14.5, height = 11.2, dpi = 300)


#   Ag_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_KhabRisk_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   Fo_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_KhabRisk_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   In_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In, LEGEND_Title = "Koala habitat\nloss risk",
#                           COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_KhabRisk_MAX[i]),
#                           North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   All_KoalaHab_plot <- PLOTMAP_risk(DATA = SUs_Pred_SF, FILL = KhabRisk_All, LEGEND_Title = NULL,
#                           COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE), SCALE_LIMIT = c(0, Lim_KhabRisk_MAX[i]),
#                           North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
#                           XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
#                           YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
#                           FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
#   Combine_KoalaHab_plot <- ((Ag_KoalaHab_plot | In_KoalaHab_plot) / (Fo_KoalaHab_plot | All_KoalaHab_plot) +
#     plot_annotation(tag_levels = 'A') &
#     theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
#           plot.tag.location = "panel", plot.tag.position = c(0.03,0.96))) +
#     inset_element(
#       STE_SA_plot + theme(plot.tag = element_blank()),
#       left = 0.85, bottom = 0.01, right = 0.99, top = 0.34)
#   Output_FPATH_Koala <- file.path(OUTPUT_DIR, "figures", paste0("KoalaHab_map4_Lim_", c("IQR", "95Q", "99Q", "Max", "25", "50")[i], ".png"))
#   ggsave(Output_FPATH_Koala, Combine_KoalaHab_plot, width = 14.5, height = 11.2, dpi = 300)
# }




############
# Plot one map with 3 insets----
## These maps were used in supplementary materials, but not in the main text.
#### File names
All_risk_with_Insets_FPath  <- file.path(OUTPUT_DIR, "figures/Pred_All_map1.png")
All_risk_with_Insets_FPath6 <- file.path(OUTPUT_DIR, "figures/Pred_All_map6.png")
Ag_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Ag_map1.png")
Fo_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Fo_map1.png")
In_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_In_map1.png")

# Inset bounding box locations
Inset_BL_Ag <- data.frame(x = c(9395500, 9376800, 9381722), y = c(4836000, 4726000, 4560187))
Inset_BL_In <- data.frame(x = c(9650498, 9769700, 9599470), y = c(4507888, 4593377, 4352856))
Inset_BL_Fo <- data.frame(x = c(9810110, 9503821, 9327460), y = c(4913177, 4379838, 4196286))
Inset_BL_All <- bind_rows(Inset_BL_Ag[1,], Inset_BL_Fo[1,], Inset_BL_In[1,])

# Generate maps with insets
All_risk_with_Insets <- PLOTMAP_risk_with_6Insets(
  DATA_ALL = SUs_Pred_SF, FILL_ALL = Pred_All ,
  DATA_Agri = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL_Agri = Pred_Ag ,
  DATA_Fort = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL_Fort = Pred_Fo ,
  DATA_Infr = SUs_Pred_SF %>% drop_na(Pred_In), FILL_Infr = Pred_In ,
  COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_risk_with_Insets_FPath6, PNG_width = 11, PNG_height = 11, PNG_dpi = 400)

All_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF, FILL = Pred_All , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = c(0,0.1, 0.2, 0.3, 0.44), SCALE_LABELS = c("0", "0.1", "0.2", "0.3", "≥0.44"),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 400)

Ag_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = c(0,0.1, 0.2, 0.3, 0.44), SCALE_LABELS = c("0", "0.1", "0.2", "0.3", "≥0.44"),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000,
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo",
  FilenamePath_PNG = Ag_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 400)

Fo_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = c(0,0.1, 0.2, 0.3, 0.44), SCALE_LABELS = c("0", "0.1", "0.2", "0.3", "≥0.44"),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000,
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 400)

In_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_Pred_MAX_99Q), SCALE_BREAKS = c(0,0.1, 0.2, 0.3, 0.44), SCALE_LABELS = c("0", "0.1", "0.2", "0.3", "≥0.44"),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_In, Inset_dim = 100000,
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 400)

### Koala habitat loss risk maps----
# Define file paths
All_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_All_map1.png")
All_Koala_risk_Plot_FPath6 <- file.path(OUTPUT_DIR, "figures/Pred_Koala_All_map6.png")
Ag_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_Ag_map1.png")
Fo_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_Fo_map1.png")
In_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_In_map1.png")

# Generate maps with insets
All_Koala_risk_Plot <- PLOTMAP_risk_with_6Insets(
  DATA_ALL = SUs_Pred_SF, FILL_ALL = KhabRisk_All ,
  DATA_Agri = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL_Agri = KhabRisk_Ag ,
  DATA_Fort = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL_Fort = KhabRisk_Fo ,
  DATA_Infr = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL_Infr = KhabRisk_In ,
  COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_Koala_risk_Plot_FPath6, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

All_Koala_risk_Plot <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_All), FILL = KhabRisk_All , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.2439), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.24"),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"),  URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Ag_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.2439), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.24"),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000,
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo",
  FilenamePath_PNG = Ag_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Fo_Koala_risk_Plot <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.2439), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.24"),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000,
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

In_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  SCALE_LIMIT = c(0, Lim_KhabRisk_MAX_99Q), SCALE_BREAKS = c(0, 0.05, 0.1, 0.15, 0.2, 0.2439), SCALE_LABELS = c("0", "0.05", "0.1", "0.15", "0.2", "≥0.24"),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_In, Inset_dim = 100000,
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

################################################################################################################################
# Plot random effect
SA1 <- qs_read(file.path(OUTPUT_DIR, "spatial_units", "sa1s.qs2")) %>% 
  bind_rows()

Rand_CI_Ag <- qs_read(file.path(OUTPUT_DIR, "data/Rand_CI_Ag.qs2"))
Rand_CI_Fo <- qs_read(file.path(OUTPUT_DIR, "data/Rand_CI_Fo.qs2"))
Rand_CI_In <- qs_read(file.path(OUTPUT_DIR, "data/Rand_CI_In.qs2"))

Rand_CI_Ag_SF <- SA1 %>% left_join(Rand_CI_Ag, by = c("SA1_CODE21" = "SA1_CODE21")) %>% drop_na() %>% dplyr::select("BsgEFF_P", "IIDEFF_P", "BsgEFF_N", "IIDEFF_N")
Rand_CI_Fo_SF <- SA1 %>% left_join(Rand_CI_Fo, by = c("SA1_CODE21" = "SA1_CODE21")) %>% drop_na() %>% dplyr::select("BsgEFF_P", "IIDEFF_P", "BsgEFF_N", "IIDEFF_N")
Rand_CI_In_SF <- SA1 %>% left_join(Rand_CI_In, by = c("SA1_CODE21" = "SA1_CODE21")) %>% drop_na() %>% dplyr::select("BsgEFF_P", "IIDEFF_P", "BsgEFF_N", "IIDEFF_N")

NSW_urb_sel_pt2 <- NSW_urb_sel_pt2 %>% 
  mutate(x = case_when(UCL_NAME16 == "Newcastle" ~ x-100000, .default = x),
         y = case_when(UCL_NAME16 == "Newcastle" ~ y-10000, .default = y))

Ag_BsgEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_Ag_SF, FILL = BsgEFF_P, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Broc")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Broc")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Ag_IIDEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_Ag_SF, FILL = IIDEFF_P, LEGEND_Title = "Spatial unstructured\neffect (IIDs)",
                         COL_Low = hcl.colors(15, palette = "Broc")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Broc")[2],
                         North_arrow = TRUE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Ag_BsgEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_Ag_SF, FILL = BsgEFF_N, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Broc")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Broc")[2],
                         North_arrow = FALSE, Scale_bar = TRUE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Ag_IIDEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_Ag_SF, FILL = IIDEFF_N, LEGEND_Title = "Spatial\nunstructured\neffect\n(IIDs)",
                         COL_Low = hcl.colors(15, palette = "Broc")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Broc")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = c(0.92, 0.28),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Ag_RAND_plot <- ((Ag_BsgEFF_P_plot | Ag_IIDEFF_P_plot) / (Ag_BsgEFF_N_plot | Ag_IIDEFF_N_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
        plot.tag.location = "panel", plot.tag.position = c(0.03, 0.96)))
ggsave(file.path(OUTPUT_DIR, "figures/Ag_RAND_plot.png"), Ag_RAND_plot, width = 15.5, height = 11.2, dpi = 400)


In_BsgEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_In_SF, FILL = BsgEFF_P, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Blue-Red 3")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Blue-Red 3")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

In_IIDEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_In_SF, FILL = IIDEFF_P, LEGEND_Title = "Spatial unstructured\neffect (IIDs)",
                         COL_Low = hcl.colors(15, palette = "Blue-Red 3")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Blue-Red 3")[2],
                         North_arrow = TRUE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

In_BsgEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_In_SF, FILL = BsgEFF_N, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Blue-Red 3")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Blue-Red 3")[2],
                         North_arrow = FALSE, Scale_bar = TRUE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

In_IIDEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_In_SF, FILL = IIDEFF_N, LEGEND_Title = "Spatial\nunstructured\neffect\n(IIDs)",
                         COL_Low = hcl.colors(15, palette = "Blue-Red 3")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Blue-Red 3")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = c(0.92, 0.28),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

In_RAND_plot <- ((In_BsgEFF_P_plot | In_IIDEFF_P_plot) / (In_BsgEFF_N_plot | In_IIDEFF_N_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
        plot.tag.location = "panel", plot.tag.position = c(0.03, 0.96)))
ggsave(file.path(OUTPUT_DIR, "figures/In_RAND_plot.png"), In_RAND_plot, width = 15.5, height = 11.2, dpi = 400)


Fo_BsgEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_Fo_SF, FILL = BsgEFF_P, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Cork")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Cork")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Fo_IIDEFF_P_plot <- PLOTMAP_Rand(DATA = Rand_CI_Fo_SF, FILL = IIDEFF_P, LEGEND_Title = "Spatial unstructured\neffect (IIDs)",
                         COL_Low = hcl.colors(15, palette = "Cork")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Cork")[2],
                         North_arrow = TRUE, Scale_bar = FALSE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Fo_BsgEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_Fo_SF, FILL = BsgEFF_N, LEGEND_Title = "Spatial structured\neffect (Besag)",
                         COL_Low = hcl.colors(15, palette = "Cork")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Cork")[2],
                         North_arrow = FALSE, Scale_bar = TRUE, LABEL_CITYS = FALSE, LEG_POS = c(0.85, 0.25),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Fo_IIDEFF_N_plot <- PLOTMAP_Rand(DATA = Rand_CI_Fo_SF, FILL = IIDEFF_N, LEGEND_Title = "Spatial\nunstructured\neffect\n(IIDs)",
                         COL_Low = hcl.colors(15, palette = "Cork")[14], COL_MID = "white", COL_High = hcl.colors(15, palette = "Cork")[2],
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = c(0.92, 0.28),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 5000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Fo_RAND_plot <- ((Fo_BsgEFF_P_plot | Fo_IIDEFF_P_plot) / (Fo_BsgEFF_N_plot | Fo_IIDEFF_N_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25, colour = "black"),
        plot.tag.location = "panel", plot.tag.position = c(0.03, 0.96)))
ggsave(file.path(OUTPUT_DIR, "figures/Fo_RAND_plot.png"), Fo_RAND_plot, width = 15.5, height = 11.2, dpi = 400)



######################################################################################################
# PLOTMAP_risk_NSW(DATA = Pred_Ag, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 1, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_Ag_map1.png"))
# PLOTMAP_risk_NSW(DATA = Pred_Fo, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 2, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_Fo_map1.png"))
# PLOTMAP_risk_NSW(DATA = Pred_In, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 3, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_In_map1.png"))


## Koala habitat deforestation risk maps----
### load base layers and target data layers
#### Data Layers
# Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
# Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))
# Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))

# #### Base Layers
# KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp"))

# ## Data Source: https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.004July%202016?OpenDocument
# ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>%
#   st_transform(st_crs(Khab_risk_Ag))

# ## Data source https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.001July%202016?OpenDocument
# STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>%
#   st_transform(st_crs(Khab_risk_Ag))

# NSW_urb_sel_pt <- ABS_urb %>%
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
#   distinct(UCL_NAME16, .keep_all = TRUE) %>%
#   st_centroid(.) %>%
#   mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])

# NSW_urb_pt <- ABS_urb %>%
#   filter(!str_detect(UCL_NAME16, "emain")) %>%
#   st_centroid() %>%
#   mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>%
#   mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Mungindi (NSW Part) (L)" ~ "Mungindi (NSW Part)",
#                                 UCL_NAME16 == "Lismore" ~ "Lismore",
#                                 UCL_NAME16 == "Bonalbo (L)" ~ "Bonalbo",
#                                 UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
#                                 UCL_NAME16 == "Jilliby (L)" ~ "Jilliby",
#                                 .default = UCL_NAME16))

# NSW_urb_pt %>% filter(UCL_NAME16 %in% c("Oberon"))

# Inset_BL_Ag <- data.frame(x = c(9376800, 9435500, 9464100), y = c(4726000, 4826500, 4911747))
# Inset_BL_In <- data.frame(x = c(9599470, 9660498, 9769700), y = c(4352856, 4477888, 4593377))
# Inset_BL_Fo <- data.frame(x = c(9779110, 9503821, 9327460), y = c(4903177, 4379838, 4196286))

# #### File names
# FilenamePath_PNG_Ag <- file.path(OUTPUT_DIR, "figures/Khab_risk_Ag_map1.png")
# FilenamePath_PNG_Fo <- file.path(OUTPUT_DIR, "figures/Khab_risk_Fo_map1.png")
# FilenamePath_PNG_In <- file.path(OUTPUT_DIR, "figures/Khab_risk_In_map1.png")


# Ag_risk_with_Insets<- PLOTMAP_risk_with_Insets(DATA = Khab_risk_Ag, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 1,
#                                                Inset_BL = Inset_BL_Ag, Inset_dim = 100000, URB_PT_Main =NULL, URB_PT_SUB1 = "Coonamble", URB_PT_SUB2 = c("Collarenebri (L)", "Wee Waa"),  URB_PT_SUB3 = "Mungindi (NSW Part)",
#                                                FilenamePath_PNG = FilenamePath_PNG_Ag, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

# In_risk_with_Insets<- PLOTMAP_risk_with_Insets(DATA = Khab_risk_In, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 2,
#                                                Inset_BL = Inset_BL_In, Inset_dim = 100000, URB_PT_Main =NULL,
#                                                URB_PT_SUB1 = c("Blue Mountains", "Sydney", "Galston", "The Oaks"), URB_PT_SUB2 = c("Singleton", "Newcastle", "Jilliby"),  URB_PT_SUB3 =  c("Port Macquarie", "Taree", "Forster - Tuncurry"),
#                                                FilenamePath_PNG = FilenamePath_PNG_In, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

# Fo_risk_with_Insets <- PLOTMAP_risk_with_Insets(DATA = Khab_risk_Fo, FILL = KhabRisk , LEGEND_Title = "Koala habitat\nloss risk", ClearType = 3,
#                                                 Inset_BL = Inset_BL_Fo, Inset_dim = 100000, URB_PT_Main =NULL, URB_PT_SUB1 = c("Lismore", "Bonalbo"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
#                                                 FilenamePath_PNG = FilenamePath_PNG_Fo, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

# ## plot combined map ----

# gpkg_path <- file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg")
# SUs_Pred_SF <- st_read(dsn = gpkg_path, layer = "SUs_Predictions")

# summary(SUs_Pred_SF$Pred_All)

# Plot <- ggplot()+

#     geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+

#     geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
#     geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
#     scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+

#     # start a new scale
#     new_scale_colour() +

#     geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
#     geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16),
#                     fontface = "bold", nudge_y = -5, size = 3,
#                     color = "black",     # text color
#                     bg.color = "grey90", # shadow color
#                     bg.r = 0.05)+          # shadow radius

#     ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
#     ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+

#     theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#     theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(legend.position = c(0.9, 0.3))+
#     theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#     coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
# ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_combined.png"), Plot, width = 11, height = 11, dpi = 300)



# Inset_BL_In <- data.frame(x = c(9599470, 9660498, 9769700), y = c(4352856, 4477888, 4593377))

# Plot_All <- ggplot()+
#     geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
#     geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
#     scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+

#     theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#     theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(legend.position = "right")+
#     theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#     coord_sf(xlim = c(Inset_BL_In[2,1], Inset_BL_In[2,1]+100000), ylim = c(Inset_BL_In[2,2], Inset_BL_In[2,2]+100000), expand = FALSE)
# ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_ALL_b1.png"), Plot_All, width = 11, height = 11, dpi = 300)

# SUs_Pred_SF[is.na(SUs_Pred_SF$Pred_All), ]

# Plot_All <- ggplot()+
#     geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
#     geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
#     scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+

#     theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
#     theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(legend.position = "right")+
#     theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
#     coord_sf(xlim = c(Inset_BL_In[2,1], Inset_BL_In[2,1]+100000), ylim = c(Inset_BL_In[2,2], Inset_BL_In[2,2]+100000), expand = FALSE)
# ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_ALL_b.png"), Plot_All, width = 11, height = 11, dpi = 300)



# Extract numbers for reporting results ----
# Load covariate selection results
Cov_df <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs"))

## Get the covariate that is selected most times across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(!is.na(.), ., "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

### Occurance Positively associated across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(is.na(.), "0", if_else(.=="1", "1", "0")))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

### Occurance Negatively associated across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(is.na(.), "0", if_else(.=="-1", "1", "0")))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))


## Amount positive association across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(is.na(.), "0", if_else(.=="1", "1", "0")))),
         Cof_total = abs(Cof_NModel_Ag) + abs(Cof_NModel_Fo) + abs(Cof_NModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Amount negative association across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(is.na(.), "0", if_else(.=="-1", "1", "0")))),
         Cof_total = abs(Cof_NModel_Ag) + abs(Cof_NModel_Fo) + abs(Cof_NModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Occurance Positively associated across all KMR and specific Cleartype
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(.=="1", "1", "0"))),
         Cof_total = abs(Cof_NModel_Fo)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Occurance negatively associated across all KMR and specific Cleartype
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(.=="-1", "1", "0"))),
         Cof_total = abs(Cof_NModel_Fo)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Deforestation risk ----
SUs_Pred_SF <- st_read(file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg"), layer = "SUs_Predictions") %>%
  drop_na(Pred_All)
Pred_Ag <- qs_read(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs2"))
Pred_Fo <- qs_read(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs2"))
Pred_In <- qs_read(file.path(OUTPUT_DIR, "predictions/Pred_In.qs2"))

Pred_sum_df <- SUs_Pred_SF %>%
  mutate(Area = as.numeric(st_area(.)/1e4)) %>%
  st_drop_geometry() %>%
  mutate(RemWoodyHa = if_else((Woody - WoodClr)>0, Woody - WoodClr , 0) * 0.0625,
         WoodyHa = Woody * 0.0625)

Pred_sum_df2 <- data.frame(
  ClearType = c("All", "Agriculture", "Forestry", "Infrastructure"),
  Risk_0.75_WoodyHa = c(
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_All > 0.75], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Ag > 0.75], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Fo > 0.75], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_In > 0.75], na.rm = TRUE)
  ),
  Risk_0.5_WoodyHa = c(
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_All > 0.5], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Ag > 0.5], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Fo > 0.5], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_In > 0.5], na.rm = TRUE)
  ),
  Risk_0.25_WoodyHa = c(
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_All > 0.25], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Ag > 0.25], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Fo > 0.25], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_In > 0.25], na.rm = TRUE)
  ),
  Risk_0.1_WoodyHa = c(
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_All > 0.1], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Ag > 0.1], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_Fo > 0.1], na.rm = TRUE),
    sum(Pred_sum_df$Woody[Pred_sum_df$Pred_In > 0.1], na.rm = TRUE)
    ),
  Risk_0.75_Num_SUs = c(
    Pred_sum_df %>% drop_na(Pred_All) %>% filter(Pred_All > 0.75) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Ag) %>% filter(Pred_Ag > 0.75) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Fo) %>% filter(Pred_Fo > 0.75) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_In) %>% filter(Pred_In > 0.75) %>% nrow()
  ),
  Risk_0.5_Num_SUs = c(
    Pred_sum_df %>% drop_na(Pred_All) %>% filter(Pred_All > 0.5) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Ag) %>% filter(Pred_Ag > 0.5) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Fo) %>% filter(Pred_Fo > 0.5) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_In) %>% filter(Pred_In > 0.5) %>% nrow()
  ),
  Risk_0.25_Num_SUs = c(
    Pred_sum_df %>% drop_na(Pred_All) %>% filter(Pred_All > 0.25) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Ag) %>% filter(Pred_Ag > 0.25) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Fo) %>% filter(Pred_Fo > 0.25) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_In) %>% filter(Pred_In > 0.25) %>% nrow()
  ),
  Risk_0.1_Num_SUs = c(
    Pred_sum_df %>% drop_na(Pred_All) %>% filter(Pred_All > 0.1) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Ag) %>% filter(Pred_Ag > 0.1) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_Fo) %>% filter(Pred_Fo > 0.1) %>% nrow(),
    Pred_sum_df %>% drop_na(Pred_In) %>% filter(Pred_In > 0.1) %>% nrow()
  )
)
write.csv(Pred_sum_df2, file.path(OUTPUT_DIR, "figures/Deforestation_risk_summary_by_ClearType.csv"), row.names = FALSE)

Pred_Ag_all_result <- Pred_Ag %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else((Woody - Woody_Clrtype)>0, Woody - Woody_Clrtype , 0) * 0.0625)
sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.5,])

sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.25,])

sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.1,])

Ag_sum_df <- data.frame(Risk = c(0.75, 0.5, 0.25, 0.1),
                         RemWoodyHa = c(
                           sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.75], na.rm = TRUE),
                           sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.5], na.rm = TRUE),
                           sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.25], na.rm = TRUE),
                           sum(Pred_Ag_all_result$RemWoodyHa[Pred_Ag_all_result$PredAll > 0.1], na.rm = TRUE)
                         ),
                         Num_SUs = c(
                           nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.75,]),
                           nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.5,]),
                           nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.25,]),
                           nrow(Pred_Ag_all_result[Pred_Ag_all_result$PredAll > 0.1,])
                         )) %>%
  mutate(ClearType = "Agriculture")


Pred_Fo_all_result <- Pred_Fo %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else((Woody - Woody_Clrtype)>0, Woody - Woody_Clrtype , 0) * 0.0625)
sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.5,])

sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.25,])

sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.1,])

Fo_sum_df <- data.frame(Risk = c(0.75, 0.5, 0.25, 0.1),
                         RemWoodyHa = c(
                           sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.75], na.rm = TRUE),
                           sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.5], na.rm = TRUE),
                           sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.25], na.rm = TRUE),
                           sum(Pred_Fo_all_result$RemWoodyHa[Pred_Fo_all_result$PredAll > 0.1], na.rm = TRUE)
                         ),
                          Num_SUs = c(
                            nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.75,]),
                            nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.5,]),
                            nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.25,]),
                            nrow(Pred_Fo_all_result[Pred_Fo_all_result$PredAll > 0.1,])
                          )) %>%
  mutate(ClearType = "Forestry")
Fo_sum_df

Pred_In_all_result <- Pred_In %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else((Woody - Woody_Clrtype)>0, Woody - Woody_Clrtype , 0) * 0.0625)
sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.5], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.5,])

sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.25], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.25,])

sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.1], na.rm = TRUE)
nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.1,])

In_sum_df <- data.frame(Risk = c(0.75, 0.5, 0.25, 0.1),
                         RemWoodyHa = c(
                           sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.75], na.rm = TRUE),
                           sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.5], na.rm = TRUE),
                           sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.25], na.rm = TRUE),
                           sum(Pred_In_all_result$RemWoodyHa[Pred_In_all_result$PredAll > 0.1], na.rm = TRUE)
                         ),
                         Num_SUs = c(
                           nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.75,]),
                           nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.5,]),
                           nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.25,]),
                           nrow(Pred_In_all_result[Pred_In_all_result$PredAll > 0.1,])
                         )) %>%
  mutate(ClearType = "Infrastructure")
In_sum_df

## High-quality koala habitat clearing risk ----
Khab_risk_Ag <- qs_read(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs2"))
Khab_risk_Fo <- qs_read(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs2"))
Khab_risk_In <- qs_read(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs2"))

Khab_risk_Ag_all_result <- Khab_risk_Ag %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.5,])

sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.25,])

sum(Khab_risk_Ag_all_result$KhabHa[Khab_risk_Ag_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_Ag_all_result[Khab_risk_Ag_all_result$KhabRisk > 0.1,])

Khab_risk_Fo_all_result <- Khab_risk_Fo %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.5,])

sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.25,])

sum(Khab_risk_Fo_all_result$KhabHa[Khab_risk_Fo_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_Fo_all_result[Khab_risk_Fo_all_result$KhabRisk > 0.1,])

Khab_risk_In_all_result <- Khab_risk_In %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         KhabHa = Khab_P * Area)
sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.5], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.5,])

sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.25], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.25,])

sum(Khab_risk_In_all_result$KhabHa[Khab_risk_In_all_result$KhabRisk > 0.1], na.rm = TRUE)
nrow(Khab_risk_In_all_result[Khab_risk_In_all_result$KhabRisk > 0.1,])

Khab_risk_sum_df <- SUs_Pred_SF %>%
  mutate(Area = as.numeric(st_area(.)/1e4),
         RemWoodyHa = if_else((Woody - WoodClr)>0, Woody - WoodClr , 0) * 0.0625,
         WoodyHa = Woody * 0.0625) %>%
  st_drop_geometry() %>%
  mutate(KhabHa_All = KhabRisk_All * WoodyHa,
         KhabHa_Ag = KhabRisk_Ag * WoodyHa,
         KhabHa_Fo = KhabRisk_Fo * WoodyHa,
         KhabHa_In = KhabRisk_In * WoodyHa)

Khab_risk_sum_df2 <- data.frame(
  ClearType = c("All", "Agriculture", "Forestry", "Infrastructure"),
  Risk_75_KhabHa = c(
    sum(Khab_risk_sum_df$KhabHa_All[Khab_risk_sum_df$KhabRisk_All > 0.75], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Ag[Khab_risk_sum_df$KhabRisk_Ag > 0.75], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Fo[Khab_risk_sum_df$KhabRisk_Fo > 0.75], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_In[Khab_risk_sum_df$KhabRisk_In > 0.75], na.rm = TRUE)
  ),
  Risk_50_KhabHa = c(
    sum(Khab_risk_sum_df$KhabHa_All[Khab_risk_sum_df$KhabRisk_All > 0.5], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Ag[Khab_risk_sum_df$KhabRisk_Ag > 0.5], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Fo[Khab_risk_sum_df$KhabRisk_Fo > 0.5], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_In[Khab_risk_sum_df$KhabRisk_In > 0.5], na.rm = TRUE)
  ),
  Risk_25_KhabHa = c(
    sum(Khab_risk_sum_df$KhabHa_All[Khab_risk_sum_df$KhabRisk_All > 0.25], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Ag[Khab_risk_sum_df$KhabRisk_Ag > 0.25], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Fo[Khab_risk_sum_df$KhabRisk_Fo > 0.25], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_In[Khab_risk_sum_df$KhabRisk_In > 0.25], na.rm = TRUE)
  ),
  Risk_10_KhabHa = c(
    sum(Khab_risk_sum_df$KhabHa_All[Khab_risk_sum_df$KhabRisk_All > 0.1], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Ag[Khab_risk_sum_df$KhabRisk_Ag > 0.1], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_Fo[Khab_risk_sum_df$KhabRisk_Fo > 0.1], na.rm = TRUE),
    sum(Khab_risk_sum_df$KhabHa_In[Khab_risk_sum_df$KhabRisk_In > 0.1], na.rm = TRUE)
  ),
  Risk_75_NumSUs = c(
    Khab_risk_sum_df %>% drop_na(KhabRisk_All) %>% filter(KhabRisk_All > 0.75) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Ag) %>% filter(KhabRisk_Ag > 0.75) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Fo) %>% filter(KhabRisk_Fo > 0.75) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_In) %>% filter(KhabRisk_In > 0.75) %>% nrow()
  ),
  Risk_50_NumSUs = c(
    Khab_risk_sum_df %>% drop_na(KhabRisk_All) %>% filter(KhabRisk_All > 0.5) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Ag) %>% filter(KhabRisk_Ag > 0.5) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Fo) %>% filter(KhabRisk_Fo > 0.5) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_In) %>% filter(KhabRisk_In > 0.5) %>% nrow()
  ),
  Risk_25_NumSUs = c(
    Khab_risk_sum_df %>% drop_na(KhabRisk_All) %>% filter(KhabRisk_All > 0.25) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Ag) %>% filter(KhabRisk_Ag > 0.25) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Fo) %>% filter(KhabRisk_Fo > 0.25) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_In) %>% filter(KhabRisk_In > 0.25) %>% nrow()
  ),
  Risk_10_NumSUs = c(
    Khab_risk_sum_df %>% drop_na(KhabRisk_All) %>% filter(KhabRisk_All > 0.1) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Ag) %>% filter(KhabRisk_Ag > 0.1) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_Fo) %>% filter(KhabRisk_Fo > 0.1) %>% nrow(),
    Khab_risk_sum_df %>% drop_na(KhabRisk_In) %>% filter(KhabRisk_In > 0.1) %>% nrow()
  )
)
Khab_risk_sum_df2

# export as table
write.csv(Khab_risk_sum_df2, file.path(OUTPUT_DIR, "figures", "Khab_risk_summary.csv"), row.names = FALSE)


#########################################################################
# Plot covariate correlation matrix----
# Load correlation matrix data
Corr_Cont <- qread(file.path(OUTPUT_DIR, "collinearity/Corr_Cont_matrix.qs"))
# colnames(Corr_Cont) <- Cov_Fullname_LU[colnames(Corr_Cont)]
# rownames(Corr_Cont) <- Cov_Fullname_LU[rownames(Corr_Cont)]
scales::show_col(hcl.colors(15, palette = "plasma")[c(7,10,13)])
cov_levels <- sort(Cov_Fullname_LU[colnames(Corr_Cont)])

Corr_Cont_df <- as.data.frame(as.table(Corr_Cont)) %>%
  rename(Cov1 = Var1, Cov2 = Var2, Corr = Freq) %>%
  mutate(Cov1 = as.character(Cov1), Cov2 = as.character(Cov2),
         Cov1_full = Cov_Fullname_LU[Cov1],
         Cov2_full = Cov_Fullname_LU[Cov2],
         Corr = round(Corr, 2),
         corr_cat = ifelse(abs(Corr) > 0.6, "strong", ifelse(abs(Corr) > 0.5, "moderate", "weak")),
         Cov1_f = factor(Cov1_full, levels = cov_levels),
         Cov2_f = factor(Cov2_full, levels = cov_levels),
         Cov1_int = as.integer(Cov1_f),
         Cov2_int = as.integer(Cov2_f)) %>%
  dplyr::filter(Cov1_int > Cov2_int) %>%
  mutate(Cov1_lab = str_wrap(Cov1_full, width = max(nchar(Cov1_full))/2),
         Cov2_lab = str_wrap(Cov2_full, width = max(nchar(Cov2_full))/2))

Y_lab <- Corr_Cont_df

Corr_plot <- ggplot(Corr_Cont_df, aes(x = Cov1_lab, y = Cov2_lab)) +
  geom_tile(aes(fill = corr_cat), color = "white") +
  geom_text(aes(label = Corr), color = "black", size = 4) +
  scale_fill_manual(
    breaks = c("strong", "moderate", "weak"),
    values = hcl.colors(15, palette = "plasma")[c(7,10,12)]) +
  # reverse x-axis items
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y.right = element_text(hjust = 0 , size = 10), legend.position = "none",
        axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed()
Corr_plot
ggsave(filename = file.path(OUTPUT_DIR , "figures/Corr_Cont_plot.png"),
       Corr_plot, width = 10, height = 10, dpi = 300, bg = "white")

# Catogorical covariate correlation plot----

Corr_Categ <- qread(file.path(OUTPUT_DIR, "collinearity/CramersV_matrix.qs"))
Corr_Categ_levels <- sort(Cov_Fullname_LU[colnames(Corr_Categ)])
Corr_Categ_df <- as.data.frame(as.table(Corr_Categ)) %>%
  rename(Cov1 = Var1, Cov2 = Var2, Corr = Freq) %>%
  mutate(Cov1 = as.character(Cov1), Cov2 = as.character(Cov2),
         Cov1_full = Cov_Fullname_LU[Cov1],
         Cov2_full = Cov_Fullname_LU[Cov2],
         Corr = round(Corr, 2),
         corr_cat = ifelse(abs(Corr) > 0.6, "strong", ifelse(abs(Corr) > 0.5, "moderate", "weak")),
         Cov1_f = factor(Cov1_full, levels = Corr_Categ_levels),
         Cov2_f = factor(Cov2_full, levels = Corr_Categ_levels),
         Cov1_int = as.integer(Cov1_f),
         Cov2_int = as.integer(Cov2_f)) %>%
  dplyr::filter(Cov1_int > Cov2_int) %>%
  #remove forest tenure type because not applicable
  dplyr::filter(!(Cov1 %in% c("TenType", "ForTen") |
                    Cov2 %in% c("TenType", "ForTen"))) %>%
  mutate(Cov1_lab = str_wrap(Cov1_full, width = max(nchar(Cov1_full))/2),
         Cov2_lab = str_wrap(Cov2_full, width = max(nchar(Cov2_full))/2))
Corr_Categ_plot <- ggplot(Corr_Categ_df, aes(x = Cov1_lab, y = Cov2_lab)) +
  geom_tile(aes(fill = corr_cat), color = "white") +
  geom_text(aes(label = Corr), color = "black", size = 4) +
  scale_fill_manual(
    breaks = c("strong", "moderate", "weak"),
    values = hcl.colors(15, palette = "plasma")[c(7,10,12)]) +
  # reverse x-axis items
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1, size = 10),
        axis.text.y.right = element_text(hjust = 0 , size = 10), legend.position = "none",
        axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed()
Corr_Categ_plot
ggsave(filename = file.path(OUTPUT_DIR , "figures/Corr_Categ_plot.png"), Corr_Categ_plot, width = 4, height = 4, dpi = 300, bg = "white")

Corr_Cont_all_plot <- GGally::ggcorr(data = NULL, geom= "blank", cor_matrix = Corr_Cont, label = TRUE, hjust = 1, layout.exp = 2)+
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient)> 0.5))+  ## highlight variables with correlations  > 0.5 OR < -0.5
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = "none", alpha = "none")
Corr_Cont_df$Cov2_lab
ggplot(Corr_Cont_df, aes(x = Cov1_lab, y = Cov2_lab)) +
  geom_tile(aes(fill = corr_cat), color = "white") +
  geom_text(aes(label = Corr), color = "black", size = 3) +
  # geom_text(data = Y_lab, aes(x = Cov1, y = Cov2, label = label),
  #           inherit.aes = FALSE, angle = 0,
  #           hjust = 1, vjust = 0.5, size = 3) +
  scale_fill_manual(values = c("#D24E71", "#EB9619", "transparent"))
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 10),
        axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  coord_fixed()
