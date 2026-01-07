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
library(grid)
library(ggrepel)
library(cowplot)
library(patchwork)
library(shadowtext)
options(max.print = 500)
theme_set(theme_pubr())
source("functions.R")

# Prepare data ----
## Set directory
## Note: change the path to the input and output directory
INPUT_DIR <- file.path("R:/HABCLEAR22-Q5221/risk-model/input")
OUTPUT_DIR <- file.path("R:/HABCLEAR22-Q5221/risk-model/output")


INPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/input")
OUTPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/output")


## Extract long format dataframe of covariate selection results
Cov_ls_Ag_long <- Get_cov_coeff_long(ClearType = 1, OUTPUT_DIR = OUTPUT_DIR)
Cov_ls_In_long <- Get_cov_coeff_long(ClearType = 2, OUTPUT_DIR = OUTPUT_DIR)
Cov_ls_Fo_long <- Get_cov_coeff_long(ClearType = 3, OUTPUT_DIR = OUTPUT_DIR)

# Combine results from all clearing types ----
Cov_df <- full_join(Cov_ls_Ag_long, Cov_ls_Fo_long, by = c("Covariate", "kmr")) %>% 
  full_join(Cov_ls_In_long, by = c("Covariate", "kmr")) %>% 
  filter(!Covariate == "LandUse5") %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, Cof_PModel_Fo, Cof_NModel_Fo, Cof_PModel_In, Cof_NModel_In) %>% 
  arrange(kmr, Covariate)

qsave(Cov_df, file.path(OUTPUT_DIR, "data/Cov_df.qs"))
Cov_df <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs"))

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

qsave(Cov_df_long, file.path(OUTPUT_DIR, "data/Cov_df_long_ForPlotting.qs"))
qsave(Cov_df2_long, file.path(OUTPUT_DIR, "data/Cov_df2_long_ForPlotting.qs"))

Cov_df_long <- qread(file.path(OUTPUT_DIR, "data/Cov_df_long_ForPlotting.qs"))
Cov_df2_long <- qread(file.path(OUTPUT_DIR, "data/Cov_df2_long_ForPlotting.qs"))


# Cov_plot <- ggplot() +
#   geom_tile(data = Cov_df_long, aes(x = kmr, y = Covariate, fill = Cof_PModel), color = "grey80")+
  
#   # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
#   geom_polygon(data = Cov_df2_long, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel), color = "grey80")+
#   facet_wrap(vars(Model))+
#   scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
#                     breaks = c(-1, 1, 0),
#                     labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
#                     name = "Covariate selection and coefficient direction")+
#   scale_y_discrete(labels = c("(Intercept)", 
#                               "Agricultural profit", 
#                               "Parcel size", 
#                               "Distance to urban center", 
#                               "Distance to road", 
#                               "Drought", 
#                               "Ecological condition", 
#                               "Fire history", 
#                               expression("Land tenure\n(Leasehold)"), 
#                               expression("Land tenure\n(Crown purposes)"), 
#                               "Land tenure (Other crown land)", 
#                               "Land use (Production-Natural Env)", 
#                               "Land use (Production-Dryland)", 
#                               "Land use (Production-Irrigated)", 
#                               "Land use (Intensive Uses)", 
#                               "NVR regulated land", 
#                               "NVR regulated land", 
#                               "Planning zone (Environment)", 
#                               "Planning zone (Others)", 
#                               "Planning zone (Residential)", 
#                               "Planning zone (Rural)", 
#                               "Population density", 
#                               "Population Growth", 
#                               "Rainfall", 
#                               "Land value", 
#                               expression("Socio-Economic PC1\n(Lower income)"), 
#                               expression("Socio-Economic PC2\n(% Australia born, Eng. speaking)"), 
#                               expression("Socio-Economic PC3\n(% 1-parent fam. with/without children u15)"), 
#                               expression("Socio-Economic PC4\n(% coup. fam. no children u15, large household)"), 
#                               expression("Socio-Economic PC5\n(% coup. fam. with children under 15)"), 
#                               "Slope", 
#                               expression("Soil PC1\n(High bulk density, sand content)"), 
#                               expression("Soil PC2\n(High organic carbon, silt content)"), 
#                               expression("Soil PC3\n(High total nitrogen, avail. water capacity)"), 
#                               "Temperature"))+
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# Cov_plot

# ggsave(filename = file.path(OUTPUT_DIR, "figures/Cov_plot.png"), plot = Cov_plot, width = 16, height = 17.5, dpi = 300)


# Extract numbers for results  -----

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

# Covariate selection barplot ----

## Manipulate data for plotting
Cov_df_sum <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cov_Ag = Cof_PModel_Ag, Cov_Fo = Cof_PModel_Fo, Cov_In = Cof_PModel_In) %>%
  mutate(Cov_Ag = abs(as.numeric(Cov_Ag)),
         Cov_Fo = abs(as.numeric(Cov_Fo)),
         Cov_In = abs(as.numeric(Cov_In)),
         Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
                               Covariate  == "Area" ~ "Parcel size",
                               Covariate  == "DistCity" ~ "Distance to urban center",
                               Covariate  == "DistRoad" ~ "Distance to road",
                               Covariate  == "Drought1" ~ "Drought",
                               Covariate  == "EcolCond" ~ "Ecological condition",
                               Covariate  == "Fire1" ~ "Fire history",
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
                               Covariate  == "PropVal" ~ "Land value",
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
        # plot.title = element_text(hjust=0.9 , vjust = -10)
        )
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
ggsave(filename = file.path(OUTPUT_DIR, "figures/Cov_sum_plot.png"), Cov_sum_plot, width = 11, height = 6, dpi = 300, bg = "white")

####################################################################################
# PLot bar graph showing all covariates (standardized names) and all clearing types ----
Cov_df_all <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs")) %>% 
  mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
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
                               .default = Covariate))%>%
  filter(!Covariate == "(Intercept)")

# Summarise data for plotting
Cov_df_sum_all <- Cov_df_all %>% 
  dplyr::select(Covariate, kmr, Cov_Ag = Cof_PModel_Ag, Cov_Fo = Cof_PModel_Fo, Cov_In = Cof_PModel_In) %>%
  mutate(Cov_Ag = abs(as.numeric(Cov_Ag)),
         Cov_Fo = abs(as.numeric(Cov_Fo)),
         Cov_In = abs(as.numeric(Cov_In))) %>% 
  pivot_longer(cols = c("Cov_Ag", "Cov_Fo", "Cov_In"), names_to = "ClearType", values_to = "Count") %>% 
  mutate(ClearType = case_when(ClearType == "Cov_Ag" ~ "Agriculture",
                               ClearType == "Cov_Fo" ~ "Forestry",
                               ClearType == "Cov_In" ~ "Infrastructure"),
         Count = as.integer(Count), Count = if_else(!is.na(Count), Count, 0)) %>%
  group_by(Covariate, ClearType) %>%
  summarise(Count = sum(Count), .groups = "drop_last") %>%  ungroup()

# Define factor levels by total count across all clearing types and alphabetical order
Cov_lvls <- Cov_df_sum_all %>%
  group_by(Covariate) %>%
  summarise(Total = sum(Count)) %>%
  arrange(Total, desc(Covariate)) %>%
  pull(Covariate)

# define factor levels for ClearType
ClearType_lvls <- c("Infrastructure", "Forestry", "Agriculture")
Cov_df_sum_all$Covariate <- factor(Cov_df_sum_all$Covariate, levels = Cov_lvls)
Cov_df_sum_all$ClearType <- factor(Cov_df_sum_all$ClearType, levels = ClearType_lvls) 

# Plot bar graph
Cov_df_sum_all_plot <- ggplot(Cov_df_sum_all, aes(x = Count, y = Covariate, fill = ClearType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Covariates", x = "Number of times\ncovariate selected", fill = "Clearing Type")+
  scale_fill_manual(values = c("#009E73", "#E69F00", "grey30"),
                    breaks = c("Agriculture", "Forestry", "Infrastructure"))+
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme_pubr() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text  = element_text(size = 11.5, colour = "black")) +
  coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 30.5), expand = FALSE)
Cov_df_sum_all_plot

ggsave(filename = file.path(OUTPUT_DIR, "figures/Cov_df_sum_all_plot.png"), Cov_df_sum_all_plot, width = 10, height = 8, dpi = 300, bg = "white")

###################################################
# Create customised tile plot legend 
model_key_grob <- grobTree(
  # left triangle (e.g., Model A)
  polygonGrob(
    x = unit(c(0.2, 0.2, 0.8), "npc"),
    y = unit(c(0.2, 0.8, 0.8), "npc"),
    gp = gpar(fill = "white", col = "black", lwd = 2.5)
  ),
  # right triangle (e.g., Model B)
  polygonGrob(
    x = unit(c(0.2, 0.8, 0.8), "npc"),
    y = unit(c(0.2, 0.2, 0.8), "npc"),
    gp = gpar(fill = "white", col = "black", lwd = 2.5)
  ),
  textGrob("Zero-inflation\nComponent", x = unit(-1.7, "npc"),
    y = unit(0.5, "npc"),
    just = c("left", "centre"),
    gp = gpar(cex = 0.85)),
  textGrob("Amount\nComponent", x = unit(1, "npc"),
    y = unit(0.5, "npc"),
    just = c("left", "centre"),
    gp = gpar(cex = 0.85))
)
# Manipulate data for plotting
Cov_df_all2 <- Cov_df_all %>%
  uncount(weight = 3) %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls),
         x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

# Create long format for plotting
Cov_df_long <- rbind(Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                     Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                     Cov_df_all %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>% 
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df2_long <- rbind(Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag, x1, y1) %>% mutate(Model = "Agriculture") %>% rename(Cof_PModel = Cof_PModel_Ag, Cof_NModel = Cof_NModel_Ag),
                      Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo, x1, y1) %>% mutate(Model = "Forestry") %>% rename(Cof_PModel = Cof_PModel_Fo, Cof_NModel = Cof_NModel_Fo),
                      Cov_df_all2 %>% dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In, x1, y1) %>% mutate(Model = "Infrastructure") %>% rename(Cof_PModel = Cof_PModel_In, Cof_NModel = Cof_NModel_In)) %>% 
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

# Create tile plot
Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_long, aes(x = kmr, y = Covariate, fill = Cof_PModel), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_long, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel), color = "grey80")+
  facet_wrap(vars(Model), )+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and\ncoefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )+  
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        
        axis.text.y = element_text(size = 10.5, colour = "black"), 
        axis.text.x = element_text(size = 10.5, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        legend.text = element_text(size = 10.5))
  # coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)
Tile_plot
ggsave(filename = file.path(OUTPUT_DIR, "figures/Tile_plot_AllCovariates.png"), plot = Tile_plot, width = 16, height = 17.5, dpi = 300)

## Combining both plots
### Extract legends for customised arrangement
leg_left <- cowplot::get_legend(
  Tile_plot +
    theme(legend.direction = "vertical",
      legend.position = "top",
          legend.justification = "left",
          legend.box.just = "left")
)

leg_right <- cowplot::get_legend(
  Cov_df_sum_all_plot +
    theme(legend.direction = "vertical",
    legend.position = "top",
          legend.justification = "left",
          legend.box.just = "left")
)
leg_both <- cowplot::plot_grid(leg_left, leg_right, ncol = 1, align = "v")
leg_both_wrap <- patchwork::wrap_elements(full = leg_both)
ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_Legend.png"), plot = leg_both_wrap, width = 2.5, height = 2.5, dpi = 300)

## Combine tile plot and bar plot with inset legends
TileBar_plot <- (((Tile_plot+theme(legend.position = "none")) + 
    patchwork::inset_element(
    model_key_grob,
    left = 0.09, right = 0.14, 
    bottom = 0.92, top = 0.98,  
    # left = 0.14, right = 0.19, 
    # bottom = 0.605, top = 0.655, 
    align_to = "full") +
    patchwork::inset_element(
    leg_both_wrap,
    left = 0, right = 0.2,
    bottom = 0.6, top = 0.9,
    align_to = "full"      # coordinates relative to the full patchwork
  ))| 
  (Cov_df_sum_all_plot + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    theme(legend.position = "none"))) +
  plot_layout(widths = c(5,1))
TileBar_plot
ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_AllCovariates.png"), plot = TileBar_plot, width = 11, height = 8, dpi = 600)


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

# Combined tile and bar plot for each clear type----
## Prepare data ----
Cov_df <- qread(file.path(OUTPUT_DIR, "data/Cov_df.qs")) %>% 
  mutate(Covariate = case_when(Covariate  == "AgProf" ~ "Agricultural profit",
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

triangle_plot <- ggdraw() +
  draw_grob(model_key_grob) +
  theme(plot.margin = margin(0, 0, 0, 0))
theme_set()

## Agriculture ----
Cov_df_Ag <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Ag, Cof_NModel_Ag) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df_Ag_sum <- Cov_df_Ag %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Ag))), .groups = "drop_last") 

Cov_lvls <- Cov_df_Ag_sum %>%
  arrange(Cof_total, desc(Covariate)) %>%
  pull(Covariate)

Cov_df_Ag_sum <- Cov_df_Ag_sum %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df_Ag <- Cov_df_Ag %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df2_Ag <- Cov_df_Ag %>% 
  uncount(weight = 3) %>% 
  mutate(Covariate = factor(Covariate, levels = Cov_lvls),
         x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_Ag, aes(x = kmr, y = Covariate, fill = Cof_PModel_Ag), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_Ag, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Ag), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )+  
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
  coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

Bar_plot <- ggplot(data = Cov_df_Ag_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
        axis.text  = element_text(size = 11.5, colour = "black")) +
  coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

TileBar_plot_Ag <- Tile_plot+ inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
 Bar_plot + plot_layout(width = c(5,3.5))

ggsave(file.path(OUTPUT_DIR, "figures/TileBar_plot_Ag.png"), TileBar_plot_Ag, width = 11, height = 8, dpi = 300, bg = "white")

## Forestry ----
Cov_df_Fo <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_Fo, Cof_NModel_Fo) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df_Fo_sum <- Cov_df_Fo %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_Fo))))

Cov_lvls <- Cov_df_Fo_sum %>%
  arrange(Cof_total, desc(Covariate)) %>%
  pull(Covariate)

Cov_df_Fo_sum <- Cov_df_Fo_sum %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df_Fo <- Cov_df_Fo %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df2_Fo <- Cov_df_Fo %>% 
  uncount(weight = 3) %>% 
  mutate(Covariate = factor(Covariate, levels = Cov_lvls),
         x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_Fo, aes(x = kmr, y = Covariate, fill = Cof_PModel_Fo), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_Fo, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_Fo), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )+  
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
  coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 25.5), expand = FALSE)

Bar_plot <- ggplot(data = Cov_df_Fo_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
        axis.text  = element_text(size = 11.5, colour = "black")) +
  coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 25.5), expand = FALSE)
TileBar_plot_Fo <- Tile_plot + inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
  Bar_plot + plot_layout(width = c(5,3.5))
ggsave(filename = file.path(OUTPUT_DIR , "figures/TileBar_plot_Fo.png"), TileBar_plot_Fo, width = 11, height = 8, dpi = 300, bg = "white")

## Infrastructure ----
Cov_df_In <- Cov_df %>% 
  dplyr::select(Covariate, kmr, Cof_PModel_In, Cof_NModel_In) %>% 
  drop_na() %>% 
  filter(Covariate != "(Intercept)")

Cov_df_In_sum <- Cov_df_In %>% 
  group_by(Covariate) %>%
  summarise(Cof_total = sum(abs(as.numeric(Cof_PModel_In))))

Cov_lvls <- Cov_df_In_sum %>%
  arrange(Cof_total, desc(Covariate)) %>%
  pull(Covariate)

Cov_df_In_sum <- Cov_df_In_sum %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))

Cov_df_In <- Cov_df_In %>%
  mutate(Covariate = factor(Covariate, levels = Cov_lvls))


Cov_df2_In <- Cov_df_In %>% 
  uncount(weight = 3) %>% 
  mutate(Covariate = factor(Covariate, levels = Cov_lvls),
         x = as.numeric(as.factor(kmr)),
         x_c = rep(c(-0.5, 0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         y = as.numeric(as.factor(Covariate)),
         y_c = rep(c(-0.5, -0.5, 0.5), length(unique(Covariate))*length(unique(kmr))),
         x1 = x + x_c, y1 = y + y_c)

Tile_plot <- ggplot() +
  geom_tile(data = Cov_df_In, aes(x = kmr, y = Covariate, fill = Cof_PModel_In), color = "grey80")+
  
  # geom_polygon(aes(x=c(0.5,1.5,1.5), y=c(0.5,0.5,1.5)))
  geom_polygon(data = Cov_df2_In, aes(x = x1, y = y1, group = interaction(Covariate, kmr), fill = Cof_NModel_In), color = "grey80")+
  scale_fill_manual(values = c(diverge_hcl(5, palette = "Blue_Red3")[4:5], "grey40") ,
                    breaks = c(-1, 1, 0),
                    labels = c("Negative coefficient", "Positive coefficient", "Not Selected"), na.value = "grey80",
                    name = "Covariate selection and coefficient direction")+
  labs(x = "Koala Modelling Regions") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )+
  theme(axis.title.y=element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text  = element_text(size = 11.5, colour = "black"), legend.text = element_text(size = 11.5))+
  coord_cartesian(xlim = c(0.4, 9.5), ylim = c(0.4, 24.5), expand = FALSE)

Bar_plot <- ggplot(data = Cov_df_In_sum, aes(x = Cof_total, y = Covariate)) +
  geom_bar(stat = "identity", fill = diverge_hcl(7, palette = "Blue_Red3")[6]) +
  labs(x = "Number of times covariate selected") +
  scale_x_continuous(breaks = seq(0, 9 , by = 3),labels = abs(seq(0, 9 , by = 3)), 
                     minor_breaks = 0:9, guide = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", colour = "grey80"), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey80"),
        axis.text  = element_text(size = 11.5, colour = "black")) +
  coord_cartesian(xlim = c(0, 9.5), ylim = c(0.4, 24.5), expand = FALSE)
TileBar_plot_In <- Tile_plot + inset_element(triangle_plot, 1.36,1.027,1.52,1.077)+
  Bar_plot + plot_layout(width = c(5,3.5))
ggsave(filename = file.path(OUTPUT_DIR, "figures/TileBar_plot_In.png"), TileBar_plot_In, width = 11, height = 8, dpi = 300, bg = "white")

###########################################################################################################################################
# Plot map----
## Map for study area (KMR) ----
# Load SUs prediction shapefile
gpkg_path <- file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg")
SUs_Pred_SF <- st_read(dsn = gpkg_path, layer = "SUs_Predictions") %>% 
  drop_na(Pred_All)

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
         y = if_else(KMR == "(NS)", y-.5, y), 
         y = if_else(KMR == "(NT)", y+0.5, y),
         y = if_else(KMR == "(CST)", y+0.6, y),
         x = if_else(KMR == "(NC)", x + .5, x)) %>% 
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
  mutate(UCL_NAME16 = str_wrap(UCL_NAME16, width = 8)) %>%
  dplyr::select(UCL_NAME16, geometry) %>% 
  distinct(UCL_NAME16, .keep_all = TRUE) %>% 
  st_centroid(.) %>% 
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])



KMR_map <- ggplot()+
 
  # Koala Modelling Regions
  geom_sf(data = KMR_shp, fill = "lightgoldenrod1", color = "grey10", lwd = 0.2, linetype= "dotdash" )+
  
  
  # State boundary grey background
  geom_sf(data = STE1, fill = "grey80", color = "black", lwd = 0.5)+
  geom_text(data = STE1, aes(x = x, y = y, label = STATE_NAME), colour = "grey40", size = 7)+
  annotate("segment", x = 147.7, xend = 149, y = -36.5, yend = -35.7, linewidth = 1, 
           colour = "grey40", alpha = 0.9, lineend = "round")+

  geom_sf(data = STE, fill = "transparent", color = "black", lwd = 0.2)+
  
  # Urban areas points
  geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
  
  geom_shadowtext(data = KMR_shp, aes(x = x, y = y, label = KMRname2), 
                  fontface = "bold", nudge_y = 0, size = 5,
                  color = "black",     # text color
                  bg.color = "white", # shadow color
                  bg.r = 0.2)+

  geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), nudge_y = -0.15, size = 4,
                  color = "grey50",     # text color
                  bg.color = "white", # shadow color
                  bg.r = 0.15)+          # shadow radius
  # North arrow and scale bar
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    height = unit(3, "cm"),
                                    width = unit(3, "cm"),
                                    pad_y = unit(1, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering(
                                      fill = c("black", "white"),
                                      text_size = 16,
                                      line_width = 2,
                                      line_col = "black",
                                      text_col = "black"))+
  ggspatial::annotation_scale(location = "bl",
                              width_hint = 0.4,     # makes scale bar wider
                              bar_height = unit(0.5, "cm"),
                              line_width = 2,
                              pad_x = unit(0.5, "cm"),
                              text_cex = 1.2)+
  # Theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = "#C7E6F5"))+
  coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)]+c(-.2, 1.5), ylim = st_bbox(KMR_shp)[c(2,4)]+c(-.1, .1), 
           # datum = st_crs(KMR_shp),
           expand = FALSE)
ggsave(file.path(OUTPUT_DIR, "figures/KMR_map_base.png"), KMR_map, width = 11, height = 11, dpi = 300)

STE_SA_plot <- ggplot()+
  geom_sf(data = KMR_shp_dsvl, fill = "grey20", color = "black", lwd = 0.2)+
  geom_sf(data = STE1, fill = "grey60", color = "grey10", lwd = 0.2)+
  coord_sf(xlim = c(111, 156), ylim = c(-44.5, -9.5), 
           expand = FALSE)+
  # annotate("rect", xmin = st_bbox(KMR)[1]-.5, xmax = st_bbox(KMR)[3]+1,
  #               ymin = st_bbox(KMR)[2]-.5, ymax = st_bbox(KMR)[4]+.5,
  #               color = "red", linewidth = 1, fill = NA)+
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

NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Coffs Harbour",
                           "Nowra - Bomaderry", "Bega",
                           "Armidale", "Tamworth",
                           "Narrabri", "Dubbo",
                           "Orange", "Wagga Wagga",
                           "Griffith", "Deniliquin",
                           "Brewarrina (L)", 
                           "Broken Hill", "Ivanhoe (L)", "Cobar")) %>% 
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Brewarrina (L)" ~ "Brewarrina",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                .default = UCL_NAME16),
         UCL_NAME16 = str_wrap(UCL_NAME16, width = 8)) %>%
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
                                UCL_NAME16 == "Collarenebri (L)" ~ "Collarenebri",
                                UCL_NAME16 == "Port Macquarie" ~ "Port\nMacquarie",
                                UCL_NAME16 == "Forster - Tuncurry" ~ "Forster-\nTuncurry",
                                UCL_NAME16 == "Blue Mountains" ~ "Blue\nMountains",
                                .default = UCL_NAME16))

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
  FilenamePath_PNG = All_risk_with_Insets_FPath6, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

All_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF, FILL = Pred_All , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"), 
  FilenamePath_PNG = All_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Ag_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL = Pred_Ag , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000, 
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo", 
  FilenamePath_PNG = Ag_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Fo_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL = Pred_Fo , COLOUR = hcl.colors(8, palette = "Purples 2" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

In_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(Pred_In), FILL = Pred_In , COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_In, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

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
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_All, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"),  URB_PT_SUB3 = c("Singleton", "Newcastle"), 
  FilenamePath_PNG = All_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Ag_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Ag), FILL = KhabRisk_Ag , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Ag, Inset_dim = 100000, 
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo", 
  FilenamePath_PNG = Ag_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

Fo_Koala_risk_Plot <- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_Fo), FILL = KhabRisk_Fo , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = Fo_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

In_Koala_risk_Plot<- PLOTMAP_risk_with_3Insets(
  DATA = SUs_Pred_SF %>% drop_na(KhabRisk_In), FILL = KhabRisk_In , COLOUR = hcl.colors(8, palette = "Reds 2" ,rev = TRUE),
  LEGEND_Title = "Koala habitat\nloss risk", Inset_BL = Inset_BL_In, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Singleton", "Newcastle"), URB_PT_SUB2 = c("Port\nMacquarie", "Forster-\nTuncurry"),  URB_PT_SUB3 =  c("Blue\nMountains", "Sydney"),
  FilenamePath_PNG = In_Koala_risk_Plot_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)


PLOTMAP_risk_NSW(DATA = Pred_Ag, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 1, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_Ag_map1.png"))
PLOTMAP_risk_NSW(DATA = Pred_Fo, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 2, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_Fo_map1.png"))
PLOTMAP_risk_NSW(DATA = Pred_In, FILL = PredAll, LEGEND_Title = "Deforestation\nrisk", ClearType = 3, FilenamePath_PNG = file.path(OUTPUT_DIR, "figures/Pred_In_map1.png"))


## Koala habitat deforestation risk maps----
### load base layers and target data layers
#### Data Layers
Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))
Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))

#### Base Layers
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp"))

## Data Source: https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.004July%202016?OpenDocument
ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
  st_transform(st_crs(Khab_risk_Ag))

## Data source https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.001July%202016?OpenDocument
STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
  st_transform(st_crs(Khab_risk_Ag)) 

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
Inset_BL_Fo <- data.frame(x = c(9779110, 9503821, 9327460), y = c(4903177, 4379838, 4196286))

#### File names
FilenamePath_PNG_Ag <- file.path(OUTPUT_DIR, "figures/Khab_risk_Ag_map1.png")
FilenamePath_PNG_Fo <- file.path(OUTPUT_DIR, "figures/Khab_risk_Fo_map1.png")
FilenamePath_PNG_In <- file.path(OUTPUT_DIR, "figures/Khab_risk_In_map1.png")


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

## plot combined map ----

gpkg_path <- file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg")
SUs_Pred_SF <- st_read(dsn = gpkg_path, layer = "SUs_Predictions")

summary(SUs_Pred_SF$Pred_All)

Plot <- ggplot()+
    
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    
    geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
    geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+
    
    # start a new scale
    new_scale_colour() +
    
    geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
    geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), 
                    fontface = "bold", nudge_y = -5, size = 3,
                    color = "black",     # text color
                    bg.color = "grey90", # shadow color
                    bg.r = 0.05)+          # shadow radius
    
    ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
    ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = c(0.9, 0.3))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_combined.png"), Plot, width = 11, height = 11, dpi = 300)



Inset_BL_In <- data.frame(x = c(9599470, 9660498, 9769700), y = c(4352856, 4477888, 4593377))

Plot_All <- ggplot()+
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "right")+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = c(Inset_BL_In[2,1], Inset_BL_In[2,1]+100000), ylim = c(Inset_BL_In[2,2], Inset_BL_In[2,2]+100000), expand = FALSE)
ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_ALL_b1.png"), Plot_All, width = 11, height = 11, dpi = 300)

SUs_Pred_SF[is.na(SUs_Pred_SF$Pred_All), ]

Plot_All <- ggplot()+
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    geom_sf(data = SUs_Pred_SF, aes(fill = Pred_All), color = NA)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Purples 3" ,rev = TRUE), name = "Deforestation\nrisk")+
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "right")+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = c(Inset_BL_In[2,1], Inset_BL_In[2,1]+100000), ylim = c(Inset_BL_In[2,2], Inset_BL_In[2,2]+100000), expand = FALSE)
ggsave(file.path(OUTPUT_DIR, "figures/Deforestation_risk_map_ALL_b.png"), Plot_All, width = 11, height = 11, dpi = 300)



# Extract numbers for reporting results ----
## Deforestation risk ----
SUs_Pred_SF <- st_read(file.path(OUTPUT_DIR, "predictions/SUs_Predictions.gpkg"), layer = "SUs_Predictions") %>% 
  drop_na(Pred_All)
Pred_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Ag.qs"))
Pred_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Pred_Fo.qs"))
Pred_In <- qread(file.path(OUTPUT_DIR, "predictions/Pred_In.qs"))

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

Ag_sum_df2 <- st_drop_geometry(Pred_Ag_all_result) %>% 
  mutate(RiskBin = round(PredAll*1000, 0)/1000) %>%
  group_by(RiskBin) %>%
  summarise(RemWoodyHa = sum(RemWoodyHa, na.rm = TRUE),
            Num_SUs = n(), .groups = "drop_last") %>%
  ungroup() %>%
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

Fo_sum_df2 <- st_drop_geometry(Pred_Fo_all_result) %>% 
  mutate(RiskBin = round(PredAll*1000, 0)/1000) %>%
  group_by(RiskBin) %>%
  summarise(RemWoodyHa = sum(RemWoodyHa, na.rm = TRUE),
            Num_SUs = n(), .groups = "drop_last") %>% 
  ungroup() %>%
  mutate(ClearType = "Forestry")

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

In_sum_df2 <- st_drop_geometry(Pred_In_all_result) %>% 
  mutate(RiskBin = round(PredAll*1000, 0)/1000) %>%
  group_by(RiskBin) %>%
  summarise(RemWoodyHa = sum(RemWoodyHa, na.rm = TRUE),
            Num_SUs = n(), .groups = "drop_last") %>% 
  ungroup() %>% 
  mutate(ClearType = "Infrastructure")



Defor_sum_df <- bind_rows(Ag_sum_df, Fo_sum_df, In_sum_df)
Defor_sum_df

Defor_sum_df2 <- bind_rows(Ag_sum_df2, Fo_sum_df2, In_sum_df2) %>% 
  filter(RiskBin > 0)

Defor_sum_df3 <- bind_rows(st_drop_geometry(Pred_Ag_all_result) %>% mutate(ClearType = "Agriculture") %>% sample_frac(0.3),
                           st_drop_geometry(Pred_Fo_all_result) %>% mutate(ClearType = "Forestry") %>% sample_frac(0.3),
                           st_drop_geometry(Pred_In_all_result) %>% mutate(ClearType = "Infrastructure") %>% sample_frac(0.3)) 
nrow(Defor_sum_df3)                  

#Violin plot of RemWoodyHa by ClearType
ggplot(Defor_sum_df3, aes(x = ClearType, y = RemWoodyHa, fill = ClearType))+
  geom_violin(trim = FALSE, alpha = 0.5)+
  scale_y_continuous(labels = scales::comma)+
  labs(x = "Clearing type", y = "Removable woody biomass (ha)")+
  theme_pubr()

# Cumulative Proportion of parcels exceeding any risk threshold
Defor_sum_df4 <- Defor_sum_df2 %>%
  arrange(ClearType, RiskBin) %>%
  group_by(ClearType) %>%
  mutate(CumRemWoodyHa = cumsum(RemWoodyHa),
         TotalRemWoodyHa = sum(RemWoodyHa),
         CumPropRemWoodyHa = CumRemWoodyHa / TotalRemWoodyHa,
         CumNumSUs = cumsum(Num_SUs),
         TotalNumSUs = sum(Num_SUs),
         CumPropNumSUs = CumNumSUs / TotalNumSUs) %>%
  ungroup()

ggplot(Defor_sum_df4, aes(x = RiskBin, y = CumPropRemWoodyHa, color = ClearType))+
  geom_line()+
  labs(x = "Deforestation risk", y = "Cumulative proportion of removable woody biomass")+
  theme_pubr()+
  coord_cartesian(xlim = c(0,0.25))

ggplot(Defor_sum_df2, aes(x = RiskBin, y = RemWoodyHa, color = ClearType))+
  # geom_point()+
  # geom_smooth(method = "loess", se = FALSE)+
  geom_line()+
  # scale y axis to log10 scale
  # scale_y_log10()+
  scale_y_continuous(labels = scales::comma)+
  labs(x = "Deforestation risk", y = "Removable woody biomass (ha)")+
  theme_pubr()+
  coord_cartesian(xlim = c(0,0.25))

ggplot(Defor_sum_df2)+
 geom_boxplot(aes(x = ClearType, y = RemWoodyHa))+
  scale_y_continuous(labels = scales::comma)+
  labs(x = "Clearing type", y = "Removable woody biomass (ha)")+
  theme_pubr()

ggplot(Defor_sum_df3, aes(x = PredAll, y = RemWoodyHa, colour = ClearType))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "loess", se = FALSE)+
  labs(x = "Deforestation risk", y = "Removable woody biomass (ha)") +
  scale_y_continuous(labels = scales::comma)+
  theme_pubr()+
  coord_cartesian(xlim = c(0,0.25))

summary(Defor_sum_df3$PredAll)

summary(Defor_sum_df3[Defor_sum_df3$ClearType == "Agriculture",]$PredAll)
summary(Defor_sum_df3[Defor_sum_df3$ClearType == "Forestry",]$PredAll)
summary(Defor_sum_df3[Defor_sum_df3$ClearType == "Infrastructure",]$PredAll)

## High-quality koala habitat clearing risk ----
Khab_risk_Ag <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Ag.qs"))
Khab_risk_Fo <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_Fo.qs"))
Khab_risk_In <- qread(file.path(OUTPUT_DIR, "predictions/Khab_risk_In.qs"))

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


# Compare Deforestation risk and Koala habitat deforestation risk ----
KhabRisk_Diff <- SUs_Pred_SF %>%
  #standardise each risk to 0-1 scale
  mutate(across(c(starts_with("Pred_"), starts_with("KhabRisk_")), ~(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))),
         dKhabRisk_All = Pred_All - KhabRisk_All,
         dKhabRisk_Ag = Pred_Ag - KhabRisk_Ag,
         dKhabRisk_Fo = Pred_Fo - KhabRisk_Fo,
         dKhabRisk_In = Pred_In - KhabRisk_In) %>% 
  dplyr::select(dKhabRisk_All, dKhabRisk_Ag, dKhabRisk_Fo, dKhabRisk_In)

summary(KhabRisk_Diff)

## PLot output maps ----
#### Base Layers
KMR_shp <- st_read(file.path(INPUT_DIR, "spatial_units/biodiversity_nsw_koala_modelling_regions_v1p1/NSW_Koala_Modelling_Regions_v1.1.shp"))

ABS_urb <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_UCL_shape/UCL_2016_AUST.shp") %>% 
  st_transform(st_crs(Pred_Ag))

STE <- st_read("D:/Data/NSW_Deforestation/risk-model-covariates/Input/2016_STE_shape/STE_2016_AUST.shp") %>% 
  st_transform(st_crs(Pred_Ag))

NSW_urb_sel_pt <- ABS_urb %>%
  filter(UCL_NAME16 %in% c("Coffs Harbour",
                           "Nowra - Bomaderry", "Bega",
                           "Armidale", "Tamworth",
                           "Narrabri", "Dubbo",
                           "Orange", "Wagga Wagga",
                           "Griffith", "Deniliquin",
                           "Brewarrina (L)", 
                           "Broken Hill", "Ivanhoe (L)", "Cobar")) %>% 
  mutate(UCL_NAME16 = case_when(UCL_NAME16 == "Brewarrina (L)" ~ "Brewarrina",
                                UCL_NAME16 == "Ivanhoe (L)" ~ "Ivanhoe",
                                .default = UCL_NAME16),
         UCL_NAME16 = str_wrap(UCL_NAME16, width = 8)) %>%
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
                                UCL_NAME16 == "Collarenebri (L)" ~ "Collarenebri",
                                .default = UCL_NAME16))

NSW_urb_pt %>% filter(UCL_NAME16 %in% c("Oberon"))

#### File names
dAll_risk_with_Insets_FPath  <- file.path(OUTPUT_DIR, "figures/dPred_All_map1.png")
dAg_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/dPred_Ag_map1.png")
dFo_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/dPred_Fo_map1.png")
dIn_risk_with_Insets_FPath <- file.path(OUTPUT_DIR, "figures/dPred_In_map1.png")

# Inset_BL_All <- data.frame(x = c(9376800, 9820110, 9650498), y = c(4726000, 4903177, 4507888))

Inset_BL_Ag <- data.frame(x = c(9395500, 9376800, 9381722), y = c(4836000, 4726000, 4560187))
Inset_BL_In <- data.frame(x = c(9650498, 9769700, 9599470), y = c(4507888, 4593377, 4352856))
Inset_BL_Fo <- data.frame(x = c(9810110, 9503821, 9327460), y = c(4913177, 4379838, 4196286))
Inset_BL_All <- bind_rows(Inset_BL_Ag[1,], Inset_BL_Fo[1,], Inset_BL_In[1,])


# dAll_risk_with_Insets <- PLOTMAP_risk_with_6Insets(
#   DATA_ALL = SUs_Pred_SF, FILL_ALL = Pred_All ,
#   DATA_Agri = SUs_Pred_SF %>% drop_na(Pred_Ag), FILL_Agri = Pred_Ag ,
#   DATA_Fort = SUs_Pred_SF %>% drop_na(Pred_Fo), FILL_Fort = Pred_Fo ,
#   DATA_Infr = SUs_Pred_SF %>% drop_na(Pred_In), FILL_Infr = Pred_In , 
#   COLOUR = hcl.colors(8, palette = "Purples 3" ,rev = TRUE),
#   LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000, 
#   URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"), 
#   FilenamePath_PNG = dAll_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

dAll_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = KhabRisk_Diff, FILL = dKhabRisk_All , COLOUR = hcl.colors(8, palette = "Greens 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_All, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Walgett"), URB_PT_SUB2 = c("Bonalbo", "Casino"), URB_PT_SUB3 = c("Singleton", "Newcastle"), 
  FilenamePath_PNG = dAll_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

dAg_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = KhabRisk_Diff %>% drop_na(dKhabRisk_Ag), FILL = dKhabRisk_Ag , COLOUR = hcl.colors(8, palette = "Greens 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Ag, Inset_dim = 120000, 
  URB_PT_SUB1 = "Walgett", URB_PT_SUB2 = c("Coonamble"),  URB_PT_SUB3 = "Dubbo", 
  FilenamePath_PNG = dAg_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

dFo_risk_with_Insets <- PLOTMAP_risk_with_3Insets(
  DATA = KhabRisk_Diff %>% drop_na(dKhabRisk_Fo), FILL = dKhabRisk_Fo , COLOUR = hcl.colors(8, palette = "Greens 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_Fo, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Bonalbo", "Casino"), URB_PT_SUB2 = "Oberon",  URB_PT_SUB3 =  c("Holbrook", "Tumbarumba"),
  FilenamePath_PNG = dFo_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)

dIn_risk_with_Insets<- PLOTMAP_risk_with_3Insets(
  DATA = KhabRisk_Diff %>% drop_na(dKhabRisk_In), FILL = dKhabRisk_In , COLOUR = hcl.colors(8, palette = "Greens 3" ,rev = TRUE),
  LEGEND_Title = "Deforestation\nrisk", Inset_BL = Inset_BL_In, Inset_dim = 100000, 
  URB_PT_SUB1 = c("Port Macquarie", "Taree", "Forster - Tuncurry"), URB_PT_SUB2 = c("Singleton", "Newcastle", "Jilliby"),  URB_PT_SUB3 =  c("Blue Mountains", "Sydney", "Galston", "The Oaks"),
  FilenamePath_PNG = dIn_risk_with_Insets_FPath, PNG_width = 11, PNG_height = 11, PNG_dpi = 600)


Inset_BL <- Inset_BL
ggplot(data = NULL)+
  geom_rect(aes(xmin = Inset_BL[1,1], xmax = Inset_BL[1,1]+100000, ymin = Inset_BL[1,2], ymax = Inset_BL[1,2]+100000), fill = "transparent", colour = "black", linewidth = 0.6)

