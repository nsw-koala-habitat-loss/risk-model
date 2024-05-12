# THIS CODE READS IN THE REQUIRED DATA AND ORGANISES THE DATA FOR MODEL FITTING

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

# get woody cover in 2011
Woody <- terra::rast("input/woody_cover/woody_nsw.tif") %>% round()

# load woody loss rasters and create raster stack

# get all tif files
WLFiles <- list.files("./input/clearing_data", pattern = "\\.tif$", full.names = TRUE)
WLStack <- terra::rast(WLFiles) %>% round()

# slpit into the different types of clearing
Ag_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 1, 0, 0))) %>% max(na.rm = TRUE) %>% crop(Woody)
In_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 0, 1, 0))) %>% max(na.rm = TRUE) %>% crop(Woody)
Fo_Loss <- WLStack %>% classify(cbind(c(1, 2, 3, 4), c(0, 0, 0, 1))) %>% max(na.rm = TRUE) %>% crop(Woody)

# create raster stack of woody extent and loss
Stack <- terra::rast(list(aloss = Ag_Loss, iloss = In_Loss, floss = Fo_Loss, woody = Woody))

# save raster stack
saveRDS(Stack, file = "output/raster_stacks/woodyextloss.rds")

# proposed covariates based on workshops

#Agriculture:
#Land use
#Combined Drought Indicator
#Property size
#Property value
#Distance to nearest SUA (significant urban area)
#Ecological condition
#Precipitation
#Temperature
#Slope
#Soil fertility

#Infrastructure:
#Land use
#Property size
#Property value
#Distance to nearest SUA
#Population density
#Ecological condition
#Slope

#Forestry:
#Land use
#Combined Drought Indicator
#Forest tenure
#Property size
#Property value
#Household income
#Ecological condition
#Fire
#Slope

# load covariates
# PROBABLY NEEDS TO CHANGE - FEI TO LOOK AT

# combined drought indicator
Drought <- terra::rast("input/covariates/drought.tif")
# elevation
Elev <- terra::rast("input/covariates/elev.tif")
# fire history (2011-2019)
Fire <- terra::rast("input/covariates/fire.tif")
# forest code
ForCode <- terra::rast("input/covariates/forest_code.tif")
# forest tenure
ForTen <- terra::rast("input/covariates/forest_tenure.tif")
# forest tenure type
ForTenType <- terra::rast("input/covariates/forest_tenure_type.tif")
# 2016 median household income weekly
Income <- terra::rast("input/covariates/income.tif")
# land use (NSW landuse 2013)
LandUse <- terra::rast("input/covariates/landuse.tif")
# population density
PopDen <- terra::rast("input/covariates/pop_den.tif")
# precipitation
Precip <- terra::rast("input/covariates/prec.tif")
# remoteness structure (ABS)
Remote <- terra::rast("input/covariates/remoteness.tif")
# slope
Slope <- terra::rast("input/covariates/slope.tif")
# soil fertility
SoilFert <- terra::rast("input/covariates/soil_fert.tif")
# soil nitrogen
SoilNit <- terra::rast("input/covariates/soil_nitrogen.tif")
# soil type
SoilType <- terra::rast("input/covariates/soil_type.tif")
# temperature
Temp <- terra::rast("input/covariates/temp.tif")
# # Property size
# PropSize <- terra::rast("input/covariates/prop_size.tif")
# # Property value
# PropVal <- terra::rast("input/covariates/prop_value.tif")

# reclassify discrete covariates as appropriate

# FEI TO WRITE CODE TO DO THIS - E.G., SEE LINES 35-37

# create raster stack of continuous covariates
StackCovsC <- terra::rast(list(Elev = Elev, Income = Income, PopDen = PopDen, Precip = Precip, Slope = Slope, SoilNit = SoilNit, Temp = Temp))

# save raster stack
saveRDS(StackCovsC, file = "output/raster_stacks/cont_covs.rds")

# create raster stack of discrete covariates
StackCovsD <- terra::rast(list(Drought = Drought, Fire = Fire, ForCode = ForCode, ForTen = ForTen, ForTenType = ForTenType, LandUse = LandUse, Remote = Remote, SoilFert = SoilFert, SoilType = SoilType))

# save raster stack
saveRDS(StackCovsD, file = "output/raster_stacks/disc_covs.rds")

# load spatial property units for each KMR
SUs <- list(CC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Coast"),
          CST = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Southern_Tablelands"),
          DRP = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Darling_Riverine_Plains"),
          FW = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Far_West"),
          NC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "North_Coast"),
          NT = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northern_Tablelands"),
          NS = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northwest_Slopes"),
          R = st_read("input/spatial_units/lots_kmrs.gdb", layer = "Riverina"),
          SC = st_read("input/spatial_units/lots_kmrs.gdb", layer = "South_Coast"))

# add area to the spatial units attribute table
for (i in names(SUs)) {
  SUs[[i]] <- SUs[[i]] %>% mutate(Area = Shape_Area / 10000)
}

# save processes spatial units
saveRDS(SUs, file = "output/spatial_units/sus.rds")

# load SA1s spatial layer
# note that these are the SA1s for the whole study area, not for each KMR
SA1s_All <- st_read("input/spatial_units/sa1s.gdb", layer = "sa1s")

# dissaggregate SA1s into individual KMRs
SA1s <- list(CC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$CC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            CST = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$CST$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            DRP = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$DRP$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            FW = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$FW$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            NC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            NT = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NT$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            NS = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$NS$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            R = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$R$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),],
            SC = SA1s_All[which(!is.na(left_join(SA1s_All, as_tibble(unique(SUs$SC$SA1)), join_by(SA1_CODE21 == value), keep = TRUE)$value)),])

# save processed SA1s
saveRDS(SA1s, file = "output/spatial_units/sa1s.rds")

# crop woody extent and loss rasters by spatial units for each KMR
CropRast <- map(.x = SUs, .f = get_crop, Raster = Stack)

# calculate the number of cells of woody extent and loss in each spatial unit for each KMR
ZStats_Woody <- map2(.x = SUs, .y = CropRast, .f = get_zonal, Stat = "sum")

# round to the nearest integer
ZStats_Woody <- map(.x = ZStats_Woody, .f = round)

# save data
saveRDS(ZStats_Woody, file = "output/data/ZStats_Woody.rds")

# crop continuous covariate rasters by spatial units for each KMR
CropRastCovsC <- map(.x = SUs, .f = get_crop, Raster = StackCovsC)

# crop discrete covariate rasters by spatial units for each KMR
CropRastCovsD <- map(.x = SUs, .f = get_crop, Raster = StackCovsD)

# get continuous covariate values
ZStats_CovsC <- map2(.x = SUs, .y = CropRastCovsC, .f = get_zonal, Stat = "mean")

# remove "mean" label from continuous covariates names
for (i in names(ZStats_CovsC)) {
  names(ZStats_CovsC[[i]]) <- names(ZStats_CovsC[[i]]) %>% str_remove("mean.")
}

# save data
saveRDS(ZStats_CovsC, file = "output/data/ZStats_CovsC.rds")

# get discrete covariate values
ZStats_CovsD <- map2(.x = SUs, .y = CropRastCovsD, .f = get_zonal, Stat = "mode") # note here could use Stat = "frac" to get the fraction of each discrete type in each property

# remove "mode" label from discrete covariates names and covert to factors
for (i in names(ZStats_CovsD)) {
  names(ZStats_CovsD[[i]]) <- names(ZStats_CovsD[[i]]) %>% str_remove("mode.")
  ZStats_CovsD[[i]] <- ZStats_CovsD[[i]] %>% mutate_all(~as.factor(.))
}

# save data
saveRDS(ZStats_CovsD, file = "output/data/ZStats_CovsD.rds")


######################################
# For checking normality of the data #
######################################
library(MASS)
library(ggpubr)
library(ggokabeito)
theme_set(theme_pubr())

# Read in all spatial units
SUs <- readRDS("output/spatial_units/sus.rds")

# Read in all continuous covariates
ZStats_CovsC <- readRDS("output/data/ZStats_CovsC.rds")

# Check for negative value in  continuous covariates 
ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC)
map_dbl(ZStats_CovsC_all, min, na.rm = TRUE)

# Include area in Continuous Covariates, remove units with area<0, +1+min(x) to all values
for (i in names(ZStats_CovsC)) {
  ZStats_CovsC[[i]] <- ZStats_CovsC[[i]] %>%
    mutate(Area = SUs[[i]]$Area,
           Elev = Elev+abs(floor(min(ZStats_CovsC_all$Elev, na.rm = TRUE))),
           Income = Income+1,
           PopDen = PopDen+1,
           Slope = Slope+1) %>% 
    filter(Area > 0)
}

# bind all continuous covariates for all spatial units
ZStats_CovsC_all <- ZStats_CovsC
for (i in names(ZStats_CovsC_all)) {
  ZStats_CovsC_all[[i]] <- ZStats_CovsC_all[[i]] %>% mutate(KMR = i)
}
ZStats_CovsC_all <- do.call(rbind, ZStats_CovsC_all)

# Create the histogram plot with facet_wrap
hist_plot <- ggplot(ZStats_CovsC_long, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  labs(title = "Histograms of Variables", x = "Value", y = "Frequency") +
  facet_wrap(~name, scales = "free")

# View the histogram plot
hist_plot

hist_Elev <- ggplot(ZStats_CovsC_all, aes(x = Elev, fill = KMR)) +
  geom_histogram(binwidth = 10) +
  scale_fill_okabe_ito()+
  labs(x = "Elevation", y = "Frequency")+
  theme(legend.position = "none")
# hist_Elev

hist_Income <- ggplot(ZStats_CovsC_all, aes(x = Income, fill = KMR)) +
  geom_histogram(binwidth = 100) +
  scale_fill_okabe_ito()+
  labs(x = "Income", y = "Frequency")+
  theme(legend.position = "none")
# hist_Income

hist_PopDen <- ggplot(ZStats_CovsC_all, aes(x = PopDen, fill = KMR)) +
  geom_histogram(binwidth = 1000) +
  scale_fill_okabe_ito()+
  labs(x = "Population Density", y = "Frequency")+
  theme(legend.position = "none")
# hist_PopDen

hist_Precip <- ggplot(ZStats_CovsC_all, aes(x = Precip, fill = KMR)) +
  geom_histogram(binwidth = 10) +
  scale_fill_okabe_ito()+
  labs(x = "Precipitation", y = "Frequency")+
  theme(legend.position = "none")
hist_Precip

hist_Slope <- ggplot(ZStats_CovsC_all, aes(x = Slope, fill = KMR)) +
  geom_histogram(binwidth = 1) +
  scale_fill_okabe_ito()+
  labs(x = "Slope", y = "Frequency")+
  theme(legend.position = "none")
# hist_Slope

hist_SoilNit <- ggplot(ZStats_CovsC_all, aes(x = SoilNit, fill = KMR)) +
  geom_histogram(binwidth = 0.01) +
  scale_fill_okabe_ito()+
  labs(x = "Soil Nitrogen", y = "Frequency")+
  theme(legend.position = "none")
# hist_SoilNit

hist_Temp <- ggplot(ZStats_CovsC_all, aes(x = Temp, fill = KMR)) +
  geom_histogram(binwidth = 1) +
  scale_fill_okabe_ito()+
  labs(x = "Temperature", y = "Frequency")+
  theme(legend.position = "none")
# hist_Temp

hist_Area <- ggplot(ZStats_CovsC_all, aes(x = Area, fill = KMR)) +
  geom_histogram(binwidth = 10000) +
  scale_fill_okabe_ito()+
  labs(x = "Area", y = "Frequency")+
  theme(legend.position = "none")
# hist_Area

hist_all <- ggarrange(hist_Elev, hist_Income, hist_PopDen, hist_Precip, hist_Slope, hist_SoilNit, hist_Temp, hist_Area, 
                      ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom")
ggsave("hist_ZStats_CovsC.png", hist_all, width = 3000, height = 2000, units = "px", dpi = 300)



ZStats_CovsC_LB <- data.frame(matrix(nrow = length(ZStats_CovsC), ncol = ncol(ZStats_CovsC[[1]])))
colnames(ZStats_CovsC_LB) <- colnames(ZStats_CovsC_CC)
rownames(ZStats_CovsC_LB) <- names(ZStats_CovsC)

for(j in 1:length(ZStats_CovsC)){
  for (i in 1:8) {
    lm_mod <- lm(ZStats_CovsC[[j]][[i]] ~ 1, na.action = na.omit)
    b <- boxcox(lm_mod, plotit = FALSE)
    ZStats_CovsC_LB[j,i] <- b$x[which.max(b$y)]
  }
}

ZStats_CovsC_all_LB <- data.frame(matrix(nrow = 1, ncol = ncol(ZStats_CovsC[[1]])))
colnames(ZStats_CovsC_all_LB) <- colnames(ZStats_CovsC_CC)
rownames(ZStats_CovsC_all_LB) <- "allKMR"

for (i in 1:8) {
    lm_mod <- lm(ZStats_CovsC_all[[i]] ~ 1, na.action = na.omit)
    b <- boxcox(lm_mod, plotit = FALSE)
    ZStats_CovsC_all_LB[1,i] <- b$x[which.max(b$y)]
}

ZStats_CovsC_LB <- rbind(ZStats_CovsC_LB, ZStats_CovsC_all_LB) %>% signif(1)

ZStats_CovsC_LB_cl <- ZStats_CovsC_LB %>% signif(1) %>% 
  mutate_all(function(x) {
    ifelse(x <= -1.7, "1/x^2",
           ifelse(x >= -1.3 & x <= -0.8, "1/x",
                  ifelse(x >= -0.7 & x <= -0.3, "1/sqrt(x)",
                         ifelse(x >= -0.2 & x <= 0.2, "log(x)",
                                ifelse(x >= 0.3 & x <= 0.6, "sqrt(x)",
                                       ifelse(x >= 0.7 & x <= 1.3, "x", 
                                              ifelse(x >= 1.7, "x^2", x)))))))
  })

ZStats_CovsC_LB_tab <- ggtexttable(ZStats_CovsC_LB, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>% 
  tab_add_hline(at.row = 10:11, row.side = "bottom", linewidth = 1) %>% 
  tab_add_vline(at.column = c(1), column.side = "right", from.row = 2)
ggsave("ZStats_CovsC_LB_tab.png", ZStats_CovsC_LB_tab, width = 3000, height = 2000, units = "px", dpi = 300, bg = 'white')

ZStats_CovsC_LB_cl_tab <- ggtexttable(ZStats_CovsC_LB_cl, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>% 
  tab_add_hline(at.row = 10:11, row.side = "bottom", linewidth = 1) %>% 
  tab_add_vline(at.column = c(1), column.side = "right", from.row = 2)
ggsave("ZStats_CovsC_LB_cl_tab.png", ZStats_CovsC_LB_cl_tab, width = 3000, height = 2000, units = "px", dpi = 300, bg = 'white')

