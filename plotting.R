## Plotting model selection results

#Load library & functions ----
library(tidyverse)
library(terra)
library(tidyterra)
# library(qs)
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


KMRshort_LU <- setNames(c("CC", "CST", "DRP", "FW", "NC", "NS", "NT", "R", "SC"),
  c("Central Coast", "Central & Southern Tablelands", "Darling Riverine Plains", "Far West", 
  "North Coast", "Northwest Slopes", "Northern Tablelands", "Riverina", "South Coast"))

LogLinEff_ALL <- qs_read(file.path(OUTPUT_DIR, "data/LogLinEff_ALL.qs2")) %>% 
  mutate(CovTypes = case_when(Covariate %in% c("DistRoad", "DistCity", "Area") ~ "Accessibility",
                              Covariate %in% c("Drought1", "Elevation", "EcolCond", "Fire1", "ForType", "Precip", "slope", "Soil_PC1", "Soil_PC2", "Soil_PC3", "Temp") ~ "Biophysical",
                              Covariate %in% c("PolPref", "PopDen", "PopGro", "ScEc_PC1", "ScEc_PC2", "ScEc_PC3", "ScEc_PC4", "ScEc_PC5") ~ "Demography",
                              Covariate %in% c("PropVal", "AgProf") ~ "Economic",
                              Covariate %in% c("LandUse1", "LandUse2", "LandUse3", "LandUse4", "NatVegReg1", "NatVegReg2", "PlanZone1", "PlanZone2", "PlanZone3", "PlanZone4") ~ "Land-use\nand\nplanning",
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
                          Covariate  == "LandUse1" ~ "Land use: Natural",
                          Covariate  == "LandUse2" ~ "Land use: Dryland",
                          Covariate  == "LandUse3" ~ "Land use: Irrigated",
                          Covariate  == "LandUse4" ~ "Land use: Intensive",
                          Covariate  == "NatVegReg1" ~ "Native vegetation regulated land",
                          Covariate  == "NatVegReg2" ~ "Native vegetation regulated land",
                          Covariate  == "PlanZone1" ~ "Planning zone (Environment)",
                          Covariate  == "PlanZone2" ~ "Planning zone (Others)",
                          Covariate  == "PlanZone3" ~ "Planning zone (Residential)",
                          Covariate  == "PlanZone4" ~ "Planning zone (Rural)",
                          Covariate  == "PopDen" ~ "Population density",
                          Covariate  == "PopGro" ~ "Population Growth",
                          Covariate  == "Precip" ~ "Rainfall",
                          Covariate  == "PropVal" ~ "Land value",
                          Covariate  == "ScEc_PC1" ~ "Socioeconomic PC1:Lower income",
                          Covariate  == "ScEc_PC2" ~ "Socioeconomic PC2:English-speaking",
                          Covariate  == "ScEc_PC3" ~ "Socioeconomic PC3:Single parent",
                          Covariate  == "ScEc_PC4" ~ "Socioeconomic PC4:Couple,no children",
                          Covariate  == "ScEc_PC5" ~ "Socioeconomic PC5:Couple,children",
                          Covariate  == "slope" ~ "Slope",
                          Covariate  == "Soil_PC1" ~ "Soil PC1: Bulk density & sand",
                          Covariate  == "Soil_PC2" ~ "Soil PC2: Organic carbon & silt",
                          Covariate  == "Soil_PC3" ~ "Soil PC3: Nitrogen & water capacity",
                          Covariate  == "Temp" ~ "Temperature",
                          .default = Covariate),
    KMRshort = recode(KMR, !!!KMRshort_LU)) 
# View(LogLinEff_ALL)
# LogLinEff_ALL %>% group_by(Covariate) %>% summarise(n_Positve = sum(mean_LogLinComb > 0, na.rm = TRUE), n_Negative = sum(mean_LogLinComb < 0, na.rm = TRUE), n_Total = sum(!is.na(mean_LogLinComb))) %>% 
#   arrange(desc(n_Negative))

# Summarise data for plotting
LogLinEff_ALL_sum <- expand_grid(Covariate = unique(LogLinEff_ALL$Covariate), ClearType = c("Ag", "Fo", "In"), KMRshort = unique(LogLinEff_ALL$KMRshort)) %>% 
    left_join(LogLinEff_ALL %>% select(Covariate, ClearType, KMRshort, mean_LogLinComb), by = c("Covariate", "ClearType", "KMRshort")) %>% 
    left_join(LogLinEff_ALL %>% select(Covariate, CovTypes) %>% distinct(), by = "Covariate")  %>% 
    mutate(Count = if_else(!is.na(mean_LogLinComb), 1L, 0L), 
           ClearType = case_when(ClearType == "Ag" ~ "Agriculture",
                                 ClearType == "Fo" ~ "Forestry",
                                 ClearType == "In" ~ "Infrastructure")) %>%
    group_by(Covariate, ClearType, CovTypes) %>%
    summarise(Count = sum(Count), .groups = "drop_last") %>% ungroup()

# Define factor levels by total count across all clearing types and alphabetical order
ClearType_lvls <- c("Infrastructure", "Forestry", "Agriculture")
CovType_lvls <- c("Accessibility", "Biophysical", "Demography", "Economic", "Land-use\nand\nplanning")
Cov_lvls <- LogLinEff_ALL_sum %>% arrange(desc(CovTypes), desc(Covariate)) %>% pull(Covariate) %>% unique()

LogLinEff_ALL_sum$Covariate <- factor(LogLinEff_ALL_sum$Covariate, levels = Cov_lvls)
LogLinEff_ALL_sum$ClearType <- factor(LogLinEff_ALL_sum$ClearType, levels = ClearType_lvls)
LogLinEff_ALL_sum$CovTypes <- factor(LogLinEff_ALL_sum$CovTypes, levels = CovType_lvls)
# head(LogLinEff_ALL_sum)

# plot bar graph ----
LogLinEff_ALL_sum_plot <- ggplot(LogLinEff_ALL_sum, aes(x = Count, y = Covariate, fill = ClearType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Covariates", fill = "Deforestation driver", 
       x = "Selected\ntimes")+
  scale_fill_manual(values = c("#009E73", "#E69F00", "#5F1415"),
                    breaks = c("Agriculture", "Forestry", "Infrastructure"),
                    guide = guide_legend(nrow = 2, byrow = TRUE, 
                                         keywidth = unit(0.6, "cm"), keyheight = unit(0.6, "cm"),
                                         title.position = "top", theme = theme(legend.box.just = "top")))+
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

# Tile plot ----
## Create a data frame with all combinations of Covariate and ClearType, and mark whether each combination is included or excluded
LogLinEff_ALL_IncExc <- expand_grid(Covariate = unique(LogLinEff_ALL$Covariate), ClearType = c("Ag", "Fo", "In")) %>% 
    left_join(LogLinEff_ALL %>% select(Covariate, ClearType) %>% distinct() %>% mutate(IncExc = "Inc"), by = c("Covariate", "ClearType")) %>% 
    mutate(IncExc = if_else(is.na(IncExc), "Exc", IncExc))

## Create a data frame for plotting the tile plot, including the mean log-linear effect and whether each combination is included or excluded
LogLinEff_ALL_xy <- expand_grid(Covariate = unique(LogLinEff_ALL$Covariate), ClearType = c("Ag", "Fo", "In"), KMRshort = unique(LogLinEff_ALL$KMRshort)) %>% 
    left_join(LogLinEff_ALL %>% select(Covariate, ClearType, KMRshort, mean_LogLinComb), by = c("Covariate", "ClearType", "KMRshort")) %>% 
    left_join(LogLinEff_ALL %>% select(Covariate, CovTypes) %>% distinct(), by = "Covariate")  %>% 
    left_join(LogLinEff_ALL_IncExc, by = c("Covariate", "ClearType")) %>% 
    mutate(TileFill = case_when(IncExc == "Exc" ~ "Excluded",
                            mean_LogLinComb > 0 ~ "Positive effect",
                            mean_LogLinComb < 0 ~ "Negative effect",
                            .default = "Not selected"),
           ClearType = case_when(ClearType == "Ag" ~ "Agriculture",
                                 ClearType == "Fo" ~ "Forestry",
                                 ClearType == "In" ~ "Infrastructure"),
            CovTypes = factor(CovTypes, levels = CovType_lvls),
            Covariate = factor(Covariate, levels = Cov_lvls),
            kmr = factor(KMRshort),
            x = as.numeric(kmr),
            y = as.numeric(Covariate))
## Calculate the lower and upper bounds for the color scale based on the interquartile range of the mean log-linear effect
LogLinEff_Q1 <- quantile(LogLinEff_ALL_xy$mean_LogLinComb, 0.25, na.rm = TRUE)
LogLinEff_Q3 <- quantile(LogLinEff_ALL_xy$mean_LogLinComb, 0.75, na.rm = TRUE)
LogLinEff_IQR <- IQR(LogLinEff_ALL_xy$mean_LogLinComb, na.rm = TRUE)
LowerBB <- LogLinEff_Q1 - 1.5 * LogLinEff_IQR
UpperBB <- LogLinEff_Q3 + 1.5 * LogLinEff_IQR

## Create a data frame for plotting the covariate type labels and segments
CovType_label_df2 <- LogLinEff_ALL_xy %>%
  group_by(CovTypes) %>%
  summarise(x = -22, y = max(y)+ 0.25, vjust = 1, .groups = "drop_last") %>% ungroup() %>% 
  mutate(ClearType = "Agriculture")

CovType_label_seg_df2 <- LogLinEff_ALL_xy %>%
  group_by(CovTypes) %>%
  summarise(ymax = max(y) , ymin = min(y), .groups = "drop_last") %>% ungroup() %>%
  uncount(weight = 3) %>% 
  mutate(ClearType = "Agriculture", REP = rep(1:3, nrow(.)/3),
         y = if_else(REP == 1, ymin - 0.45, if_else(REP == 2, ymax + 0.45, ymin - 0.45)),
         yend = if_else(REP == 1, ymin - 0.45, ymax + 0.45),
         x = -22.3, xend = if_else(REP != 3, -7.5, x))

## Create the tile plot with covariate type labels and segments
LogLinEff_Tile <- ggplot() +
  geom_text(data = CovType_label_df2, aes(x = x, y = y, vjust = 1, label = CovTypes), fontface = "bold", hjust = 0, size = 4.75,  inherit.aes = FALSE)+
  geom_segment(data = CovType_label_seg_df2, aes(x = x, xend = xend, y = y, yend = yend), color = "grey10", inherit.aes = FALSE)+
  facet_wrap(vars(ClearType))+
  scale_x_continuous(breaks = 1:length(levels(LogLinEff_ALL_xy$kmr)), labels = levels(LogLinEff_ALL_xy$kmr)) +
  scale_y_continuous(breaks = 1:length(levels(LogLinEff_ALL_xy$Covariate)), labels = levels(LogLinEff_ALL_xy$Covariate)) +
  geom_tile(data = LogLinEff_ALL_xy, aes(x = x, y = y, fill = TileFill), color = "grey80")+
  scale_fill_manual(values = c("grey30", "grey80"), breaks = c("Not selected", "Excluded"),
                    limits = c("Not selected", "Excluded"), name = "Covariate selection",
                    guide = guide_legend(order = 1, keywidth = unit(0.6, "cm"), keyheight = unit(0.6, "cm"),
                                         title.position = "top", theme = theme(legend.box.just = "top")))+

  ggnewscale::new_scale_fill()+

  geom_tile(data = LogLinEff_ALL_xy, aes(x = x, y = y, fill = mean_LogLinComb), color = "grey80")+
  
  scale_fill_gradient2(low = hcl.colors(palette = "Blue-Red", n = 11)[10], mid = "white", 
                       high = hcl.colors(palette = "Blue-Red", n = 10)[2], midpoint = 0, na.value = "transparent",
                       limits = c(LowerBB, UpperBB), breaks = c(LowerBB, 0, UpperBB), oob = scales::squish,
                       labels = c("≤-1.18", "0", "≥1.13"), name = "Mean log-linear effect",
                       guide = guide_colorbar(order = 2, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"), 
                                             title.position = "top", label.position = "bottom", 
                                             theme = theme(legend.box.just = "top")))+
  labs(x = "Koala Modelling Regions") +
  theme(legend.position = "top", legend.justification = "left", legend.spacing.x = unit(1.5, "cm"))+
  coord_cartesian(xlim = c(0.5, 9), ylim = c(1.7, 29.2), clip = "off")+
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(0, 0, 0, 20))
# LogLinEff_Tile
# ggsave(filename = file.path(OUTPUT_DIR, "figures/LogLinEff_Tile.png"), plot = LogLinEff_Tile, width = 7, height = 10, units = "in", dpi = 300)

# Create custome legend for the tile plot
LegendKeyTxt_DF2 <- data.frame(x =     0.25, 
                              y =     0.0, 
                           label = "Selected")
LegendKeyPol_DF2 <- data.frame(x = c(-0.05, -0.05, 0.05, -0.05, 0.05, 0.05), 
                             y = c(-0.05, 0.05, 0.05, -0.05, -0.05, 0.05),
                           group =c(1,1,1,             2,2,2)) %>% 
  mutate(group = factor(group, levels = c(1, 2)))

# Create the legend key for the tile plot
LegendKey <- ggplot()+
  geom_polygon(data = LegendKeyPol_DF2, aes(x = x, y = y, group = group, fill = group),color = NA)+
  scale_fill_manual(values = c(hcl.colors(palette = "Blue-Red", n = 11)[10], hcl.colors(palette = "Blue-Red", n = 10)[2]), guide = "none")+
  geom_text(data = LegendKeyTxt_DF2, aes(x = x, y = y, label = label),size = 4.5)+
  coord_fixed(xlim = c(-0.05, 0.85), ylim = c(-0.18, 0.18), expand = FALSE, clip = "off")+
  # theme_bw()
  theme_void() +
  theme(plot.background  = element_rect(fill = NA, colour = NA), 
        panel.background = element_rect(fill = NA, colour = NA))

# Extract legends for the tile plot and bar plot
LogLinEff_Tile_leg <- cowplot::get_legend(
  LogLinEff_Tile +
    theme(legend.title.position = "top", legend.direction = "horizontal",
          legend.position = "top", legend.justification = "left", legend.box.just = "top",
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
)
LogLinEff_ALL_sum_plot_leg  <- cowplot::get_legend(
  LogLinEff_ALL_sum_plot +
    theme(legend.title.position = "top", legend.direction = "horizontal",
          legend.position = "top", legend.justification = "left", legend.box.just = "left",
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
)
# Combine legends and plots into a final figure
Legend_row <- ggdraw() +
  draw_plot(LegendKey, x = -0.007, y = -0.28, width = 0.3, height = 0.9, hjust = 0, vjust = 0) +
  draw_grob(LogLinEff_Tile_leg,  x = 0.009, y = 0.07, width = 1, height = 1, hjust = 0, vjust = 0) +
  draw_grob(LogLinEff_ALL_sum_plot_leg, x = 0.7, y = 0, width = 1, height = 1, hjust = 0, vjust = 0)

LogLinEff_Tile_xleg <- LogLinEff_Tile + theme(legend.position = "none")
LogLinEff_ALL_sum_plot_xleg <- LogLinEff_ALL_sum_plot + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y  = element_blank())

main_row <- plot_grid(LogLinEff_Tile_xleg, LogLinEff_ALL_sum_plot_xleg, nrow = 1, rel_widths = c(6, 1), align = "h", axis = "tb")

final_plot <- ggdraw()+
  draw_plot(Legend_row, x = 0, y = 0.91, width = 1, height = 0.09, hjust = 0, vjust = 0) +
  draw_plot(main_row, x = -0.02, y = 0, width = 1.02, height = 0.91, hjust = 0, vjust = 0)
ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_LogLinEff_ALL.png"), plot = final_plot, width = 8.35, height = 10, units = "in", dpi = 300, bg = "white")

#########################################################################################
# Export data for reporting ----
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
# Predictive integral transform (PIT) histogram for model checking

PIT_ls <- qs_read(file.path(OUTPUT_DIR, "data/PIT.qs2"))

ClearType_LU <- c("Ag" = "Agriculture", "Fo" = "Forestry", "In" = "Infrastructure")

KMR_LU <- c("CC" = "Central Coast", 
            "CST" = "Central and Southern Tablelands",
            "DRP" = "Darling Riverine Plains",
            "FW" = "Far West",
            "NC" = "North Coast",
            "NS" = "Northwest Slopes",
            "NT" = "Northern Tablelands",
            "R" = "Riverina",
            "SC" = "South Coast")
KMR_LU_wrap <- str_wrap(KMR_LU, width = 19)
PIT_df <- enframe(PIT_ls, name = "Model", value = "PIT") %>%
  unnest(PIT) %>% 
  mutate(KMRs = str_extract(Model, "^[A-Z]+"),
         ClrTyps = str_extract(Model, "[A-Za-z]+$"),
         KMR_long = KMR_LU_wrap[KMRs], 
         KMR_long = factor(KMR_long, levels = KMR_LU_wrap),
         ClearType_long = ClearType_LU[ClrTyps]) %>% 
  select(KMR_long, ClearType_long, PIT)
str(PIT_df)
PIT_plot <- ggplot(PIT_df, aes(x = PIT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, boundary = 0, closed = "left" ) +
  geom_rug()+
  facet_grid(KMR_long ~ ClearType_long, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs( x = "Probability", y = "Density") + theme_bw()
ggsave(filename = file.path(OUTPUT_DIR, "figures/PIT_histogram3.png"), plot = PIT_plot, width = 8.27, height = 11.69, units = "in", dpi = 300)


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
  mutate(
    # x = st_coordinates(st_centroid(.))[,1], y = st_coordinates(st_centroid(.))[,2],
    KMR = c("(NC)", "(CC)", "(SC)", "(NT)", "(NS)", "(CST)", "(DRP)", "(FW)", "(R)"),
    KMRname = case_when(KMRname == "Central and Southern Tablelands" ~ "Central & Southern Tablelands", .default = KMRname),
    KMRname2 = str_wrap(paste(KMRname, KMR, sep = " "), width = 10),
    # Manual definition of coordinates for KMR labels to avoid overlap and improve visualization
    x = c(153.2559, 151.5240, 150.0433, 151.6160, 149.8047, 149.2512, 147.6922, 144.3541, 144.9514),
    y = c(-30.33956, -33.07577, -35.99513, -29.76898, -31.50941, -34.25783, -30.44102, -31.95294, -35.03332),
    x = c(153.8559,  152.4240, 150.8433,   151.6160 - 0.3,   149.8047 - 0.2,  149.2512 + 1.9,  147.6922 - 0.2,  144.3541,  144.9514),
    y = c(-30.33956, -33.17577, -37.15, -29.76898 + 1.3, -31.50941 + 2.1, -34.25783 - 1.5, -30.44102 + 0.7, -31.95294, -35.03332+0.3),
  ) %>%
  # Add in number of SUs in each KMR
  left_join(., nSUs_KMR, by = c("KMRname" = "KMR")) %>%
  mutate(KMRname3 = paste0(KMRname2, "\n(", nSUs, ")") )

KMR_shp_dsvl <- st_union(KMR_shp)

# Load state boundary shapefile
STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>%
  st_transform(st_crs("EPSG:4283"))
# Create a data frame for state boundary labeling, excluding New South Wales and Other Territories, and adding coordinates for label placement
STE1 <- STE %>%
  filter(!STATE_NAME %in% c("New South Wales", "Other Territories")) %>%
  mutate(x = c(145, 145, 135, NA,NA,NA, 147.2), y =  c(-36.5, -28.5, -44,NA,NA,NA, -36.5))  %>%
  mutate(STATE_NAME = case_when(STATE_NAME == "Australian Capital Territory" ~ "A.C.T.", .default = STATE_NAME))

# Extract urban areas that intersect with KMRs and select the two largest urban areas in each KMR based on area
URB_KMR <- st_intersection(st_centroid(ABS_urb), KMR_shp) %>%
  st_drop_geometry() %>%
  group_by(KMRname) %>%
  slice_max(n=2, order_by = AREASQKM16)

# Load Koala Habitat Suitability Map and aggregate to reduce resolution for plotting
khab <- rast(file.path(INPUT_DIR, "covariates/KoalaHabitatSuitability.tif"))
khab_only <- aggregate(khab, fact = 20, fun = "mean", cores= 16, na.rm = TRUE) %>% round() %>% as.factor()


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
         x2 = if_else(UCL_NAME16 %in% c("Casino"), x2 + 0.7, x2),
         x2 = if_else(UCL_NAME16 %in% c("Narrabri"), x2 - 0.3, x2),
         y2 = if_else(UCL_NAME16 %in% c("Broken\nHill"), y2 - 0.4, y2),
         y2 = if_else(UCL_NAME16 %in% c("Wagga\nWagga"), y2 + 0.4, y2),
         y2 = if_else(UCL_NAME16 %in% c("Narrabri", "Cobar"), y2 + 0.2, y2),
         y2 = if_else(UCL_NAME16 %in% c("Dubbo"), y2 - 0.2, y2))

# Create KMR map
KMR_map <- ggplot()+

  # Koala Modelling Regions
  geom_sf(data = KMR_shp, fill = "lightgoldenrod1", color = NA, alpha = 0.9)+
  geom_spatraster(data = khab_only, aes(fill = KHSMclass), maxcell = ncell(khab_only))+
  scale_fill_manual(values = c("1" = scales::alpha("#aa983f", 0.9)),
                    breaks = c("1"),
                    labels = c("High-quality\nkoala habitat "),
                    na.value = "transparent", name = NULL)+
  geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.3, linetype= "solid", alpha = 0.85)+

  # State boundary grey background
  geom_sf(data = STE1, fill = "grey80", color = "black", lwd = 0.5)+
  geom_text(data = STE1, aes(x = x, y = y, label = STATE_NAME), colour = "grey40", size = 7)+
  annotate("segment", x = 147.7, xend = 149, y = -36.5, yend = -35.7, linewidth = 0.6, colour = "grey40", alpha = 0.9, lineend = "round")+
  annotate("segment", x = c(151.4, 150.45), xend = c(151.8, 149.5), y = c(-28.9, -35.5), yend = c(-29.3, -34.9), linewidth = 0.6, colour = "grey40", alpha = 0.9, lineend = "round")+
  geom_sf(data = STE, fill = "transparent", color = "black", lwd = 0.2)+

  # Urban areas points
  geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1.5)+
  # Labels for KMR and urban areas
  geom_shadowtext(data = KMR_shp, aes(x = x, y = y, label = KMRname2),
                 fontface = "bold", nudge_y = 0, size = 5.5, lineheight = 0.8, 
                 color = "black", bg.color = "white", bg.r = 0.08)+
  geom_shadowtext(data = NSW_urb_sel_pt, aes(x = x2, y = y2 , label = UCL_NAME16), 
                 size = 5.2, lineheight = 0.8, color = "grey10", bg.color = "white",  bg.r = 0.05)+ 
  
  # North arrow and scale bar
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", height = unit(3, "cm"), width = unit(3, "cm"), pad_y = unit(1, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering(fill = c("black", "white"), text_size = 16, line_width = 2, line_col = "black", text_col = "black"))+
  ggspatial::annotation_scale(location = "bl", width_hint = 0.4, line_width = 2, pad_x = unit(0.5, "cm"), text_cex = 1.2)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = c(0.9, 0.32), legend.text = element_text(size = 16))+
  theme(panel.background = element_rect(fill = "#C7E6F5"))+
  theme(legend.key = element_rect(color = "black"), legend.background = element_rect(fill = "transparent", color = "transparent"))+
  coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)]+c(-.2, 1.5), ylim = st_bbox(KMR_shp)[c(2,4)]+c(-.1, .1), expand = FALSE)

# Create inset map for state of New South Wales
STE_SA_plot <- ggplot()+
  geom_sf(data = KMR_shp_dsvl, fill = "grey20", color = "black", lwd = 0.2)+
  geom_sf(data = STE1, fill = "grey60", color = "grey10", lwd = 0.2)+
  coord_sf(xlim = c(111, 156), ylim = c(-44.5, -9.5), expand = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme_void()+theme(legend.position = "none")+
  theme(plot.background  = element_rect(fill = "#91daff", color = "black", linewidth = 0.7))

KMR_map_comb <- KMR_map +
  inset_element(STE_SA_plot, left = 0.78, bottom = 0.01, right = 0.99, top = 0.32)

# Export KMR map
ggsave(file.path(OUTPUT_DIR, "figures/KMR_map2.png"), KMR_map_comb, width = 11, height = 11, dpi = 600, bg = "transparent")


## Deforestation risk maps----
### load base layers and target data  layers
#### Data Layers
SUs_Pred_SF <- st_read(file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg"), layer = "SUs_Predictions") %>%
  drop_na(Pred_All)

#### Base Layers
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp"))

ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>%
  st_transform(st_crs(SUs_Pred_SF))

STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>%
  st_transform(st_crs(SUs_Pred_SF))

STE_GDA20 <- st_read("R:/KOALACC22-Q5611/Data/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp") %>%
  filter(STE_CODE21 %in% c(1:8))
NSW_GDA2020 <- STE_GDA20 %>% filter(STE_NAME21 == "New South Wales")

NSW <- STE %>% filter(STATE_NAME == "New South Wales") %>% 
  rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE) %>% st_make_valid()
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
# Create inset map for state of New South Wales
STE_SA_plot <- ggplot()+
  geom_sf(data = STE_GDA20, fill = "grey70", color = "grey10", lwd = 0.5)+
  geom_sf(data = NSW_GDA2020, fill = "grey20", color = "black", lwd = 0.5)+
  coord_sf(xlim = c(111, 156), ylim = c(-44.5, -9.5), expand = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme_void()+theme(legend.position = "none")+
  theme(plot.background  = element_rect(fill = "#91daff", color = "black", linewidth = 0.7))

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

# Create tags for each sub plots
TAG_DF <- data.frame(TAG = paste0(LETTERS[1:4], "."), 
                     ClearTypes = c("Agriculture", "Infrastructure", "Forestry", "All combined"),
                     x = rep(st_bbox(KMR_shp)[1], 4), 
                     y = rep(st_bbox(KMR_shp)[4], 4) + 20000,
                     Pred_Type = c("Pred_Ag", "Pred_In", "Pred_Fo", "Pred_All"),
                     KoalaHab_Type = c("KhabRisk_Ag", "KhabRisk_In", "KhabRisk_Fo", "KhabRisk_All")) %>% 
          pivot_longer(cols = c(Pred_Type, KoalaHab_Type), names_to = "Type", values_to = "Var")

# Identify the 0.25 limit for Pred_All
ecdf(SUs_Pred_SF$Pred_All) %>% plot()
F_hat <- ecdf(SUs_Pred_SF$Pred_All)
F_hat(0.25)
SUs_Pred_DF <- st_drop_geometry(SUs_Pred_SF)

GLOBAL_MAX_PRED <- max(c(SUs_Pred_SF$Pred_Ag, SUs_Pred_SF$Pred_Fo, SUs_Pred_SF$Pred_In, SUs_Pred_SF$Pred_All), na.rm = TRUE)
Upper99 <- map(c("Pred_Ag", "Pred_Fo", "Pred_In", "Pred_All"), ~ quantile(SUs_Pred_SF[[.x]], 0.99, na.rm = TRUE)) %>% unlist() %>% mean()
SUs_Pred_SF2 <- SUs_Pred_SF %>% 
  sample_n(3000)

BREAKS_Pred <- c(seq(0, 0.25, length.out = 6), Inf)
LABELS_Pred <- c(paste(sprintf("%.2f", BREAKS_Pred[1:5]), sprintf("%.2f", BREAKS_Pred[2:6]), sep = "-"), paste0("≥", sprintf("%.2f", BREAKS_Pred[6])))

# Plot map for deforestation risk for each land clearing type and combined (main figure)
Ag_risk_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Fo_risk_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
In_risk_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In, LEGEND_Title = "Deforestation\nrisk",
                         COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred,
                         TAG_CT = TRUE,
                         North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
All_risk_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(Pred_All), FILL = Pred_All, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)

Combine_plot2 <- (Ag_risk_plot | In_risk_plot) / (Fo_risk_plot | All_risk_plot) +
  inset_element(
    STE_SA_plot + theme(plot.tag = element_blank()),
    left = 0.845, bottom = 0.001, right = 0.995, top = 0.34)
ggsave(file.path(OUTPUT_DIR, "figures/Pred_map4_class2.png"), Combine_plot2, width = 14.5, height = 11.85, dpi = 400)

# Plot map for koala habitat loss risk for each land clearing type and combined (main figure)
Ag_KoalaHab_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Fo_KoalaHab_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 10000),
                          YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                          FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
In_KoalaHab_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In, LEGEND_Title = "Koala habitat\nloss risk",
                         COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk,
                         TAG_CT = TRUE,
                         North_arrow = TRUE, Scale_bar = TRUE, LEG_POS = c(0.895, 0.32),
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
All_KoalaHab_plot <- PLOTMAP_risk3(DATA = SUs_Pred_SF %>% drop_na(KhabRisk_All), FILL = KhabRisk_All, LEGEND_Title = NULL,
                         COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
                         SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk,
                         TAG_CT = TRUE,
                         North_arrow = FALSE, Scale_bar = FALSE, LABEL_CITYS = TRUE, LEG_POS = "none",
                         XLIM = st_bbox(KMR_shp)[c(1,3)] + c(-10000, 175000),
                         YLIM = st_bbox(KMR_shp)[c(2,4)] + c(-5000, 70000),
                         FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300)
Combine_KoalaHab_plot2 <- (Ag_KoalaHab_plot | In_KoalaHab_plot) / (Fo_KoalaHab_plot | All_KoalaHab_plot) +
  inset_element(
    STE_SA_plot + theme(plot.tag = element_blank()),
    left = 0.845, bottom = 0.001, right = 0.995, top = 0.34)
ggsave(file.path(OUTPUT_DIR, "figures/KoalaHab_map4_class.png"), Combine_KoalaHab_plot2, width = 14.5, height = 11.85, dpi = 400)


############
# Plot one map with 3 insets (supplementary materials)----
## These maps were used in supplementary materials, but not in the main text.
#### File names
All_risk_with_Insets_FPath  <- file.path(OUTPUT_DIR, "figures/Pred_All_map3.png")
All_risk_with_Insets_FPath6 <- file.path(OUTPUT_DIR, "figures/Pred_All_map6.png")
Ag_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Ag_map3.png")
Fo_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Fo_map3.png")
In_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/Pred_In_map3.png")

# Inset bounding box locations
Inset_BL_Ag <- data.frame(x = c(9395500, 9376800, 9381722), y = c(4836000, 4726000, 4560187))
Inset_BL_In <- data.frame(x = c(9650498, 9769700, 9599470), y = c(4507888, 4593377, 4352856))
Inset_BL_Fo <- data.frame(x = c(9810110, 9503821, 9327460), y = c(4913177, 4379838, 4196286))
Inset_BL_All <- bind_rows(Inset_BL_Ag[1,], Inset_BL_Fo[1,], Inset_BL_In[1,])

GLOBAL_MAX_PRED <- max(c(SUs_Pred_SF$Pred_Ag, SUs_Pred_SF$Pred_Fo, SUs_Pred_SF$Pred_In, SUs_Pred_SF$Pred_All), na.rm = TRUE)
Upper99 <- map(c("Pred_Ag", "Pred_Fo", "Pred_In", "Pred_All"), ~ quantile(SUs_Pred_SF[[.x]], 0.99, na.rm = TRUE)) %>% unlist() %>% mean()
Upper99
SUs_Pred_SF2 <- SUs_Pred_SF %>% 
  sample_n(30000)

BREAKS_Pred <- c(seq(0, 0.25, length.out = 6), Inf)
LABELS_Pred <- c(paste(sprintf("%.2f", BREAKS_Pred[1:5]), sprintf("%.2f", BREAKS_Pred[2:6]), sep = "-"), paste0("≥", sprintf("%.2f", BREAKS_Pred[6])))
LABELS_Pred %>% rev()

# Generate maps with insets
All_risk_with_Insets <- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(Pred_All), FILL = Pred_All , COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_risk_with_Insets_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600)

Ag_risk_with_Insets <- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag , COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000,
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo",
  FilenamePath_PNG = Ag_risk_with_Insets_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)

Fo_risk_with_Insets <- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo , COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000,
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_risk_with_Insets_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)

In_risk_with_Insets<- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In , COLOUR = hcl.colors(6, palette = "Purples 3" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_Pred, SCALE_LABELS = LABELS_Pred, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_In, Inset_dim = 100000,
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_risk_with_Insets_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)


##############################################################################################################

### Koala habitat loss risk maps----
# Define file paths
All_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_All_map3.png")
All_Koala_risk_Plot_FPath6 <- file.path(OUTPUT_DIR, "figures/Pred_Koala_All_map6.png")
Ag_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_Ag_map3.png")
Fo_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_Fo_map3.png")
In_Koala_risk_Plot_FPath <- file.path(OUTPUT_DIR, "figures/Pred_Koala_In_map3.png")


Upper99_KhabRisk <- map(c("KhabRisk_Ag", "KhabRisk_Fo", "KhabRisk_In", "KhabRisk_All"), ~ quantile(SUs_Pred_SF[[.x]], 0.995, na.rm = TRUE)) %>% unlist() %>% mean()
Upper99_KhabRisk
BREAKS_KhabRisk <- c(seq(0, 0.25, length.out = 6), Inf)
LABELS_KhabRisk <- c(paste(sprintf("%.2f", BREAKS_KhabRisk[1:5]), sprintf("%.2f", BREAKS_KhabRisk[2:6]), sep = "-"), paste0("≥", sprintf("%.2f", BREAKS_KhabRisk[6])))

All_Koala_risk_Plot <- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_All), FILL = KhabRisk_All , COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_All, Inset_dim = 100000,
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"),  URB_PT_SUB3 = c("Singleton", "Newcastle"),
  FilenamePath_PNG = All_Koala_risk_Plot_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600)

Ag_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag , COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000,
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo",
  FilenamePath_PNG = Ag_Koala_risk_Plot_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)

Fo_Koala_risk_Plot <- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo , COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000,
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_Koala_risk_Plot_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)

In_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets2(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In , COLOUR = hcl.colors(6, palette = "Reds 2" ,rev = FALSE),
  SCALE_BREAKS = BREAKS_KhabRisk, SCALE_LABELS = LABELS_KhabRisk, LEG_POS = c(0.875, 0.17),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_In, Inset_dim = 100000,
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_Koala_risk_Plot_FPath, PNG_width = 19.25, PNG_height = 22, PNG_dpi = 600
)



###########################################################################################################
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
# Extract numbers for reporting results ----
# Load covariate selection results
Cov_df <- qs_read(file.path(OUTPUT_DIR, "data/Cov_df.qs2"))

## Get the covariate that is selected most times across all KMR and Clear Type
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(!is.na(.), ., "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>%
  group_by(Covariate) %>%
  summarise(Cof_total = sum(Cof_total)) %>%
  arrange(desc(Cof_total))

## Get the covariate that is selected most times across all KMR
Cov_df %>%
  mutate(across(3:8, ~as.integer(if_else(!is.na(.), ., "0"))),
         Cof_total = abs(Cof_PModel_Ag) + abs(Cof_PModel_Fo) + abs(Cof_PModel_In)) %>%
  group_by(Covariate, kmr) %>%
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
