# THIS CODE READS IN THE REQUIRED DATA AND FITS FOREST LOSS RISK MODELS FOR EACH
# KOALA MODELLING REGION (KMR) IN NEW SOUTH WALES

# load libraries
library(sf)

# load spatial units for each KMR
Central_Coast <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Coast")
Central_Southern_Tablelands <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Central_Southern_Tablelands")
Darling_Riverine_Plains <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Darling_Riverine_Plains")
Far_West <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Far_West")
North_Coast <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "North_Coast")
Northern_Tablelands <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northern_Tablelands")
Northwest_Slopes <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Northwest_Slopes")
Riverina <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "Riverina")
South_Coast <- st_read("input/spatial_units/lots_kmrs.gdb", layer = "South_Coast")
