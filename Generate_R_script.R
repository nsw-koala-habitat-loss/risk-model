# Purpose this script is for estimating HPC resource and generating SLURM script for HPC job scheduling

library(qs)
library(tidyverse)



# Step 1: Define the template ----
R_templ <-'{# This script reads in Model Fitting results and makes predictions
# mean for running as a single job on HPC (Test for array job)

# load library
library(tidyverse)
library(INLA)
library(qs)
library(sf)

set.seed(2024)

# Load directory
MODEL_DIR <- file.path("../output/models/")
PRED_DIR <- file.path("../output/predictions/")
R_DIR <- file.path("../R/")

source(file.path(R_DIR, "functions.R"))

MODEL_FPATH <- list.files(MODEL_DIR, pattern = paste0("[KMR]_[CLEARTYPE]\\\\.qs$"), full.names = TRUE)
MODEL <- qread(MODEL_FPATH)
CT <- if(MODEL$ClearType == "1"){"Ag"} else if(MODEL$ClearType == "2"){"In"} else if(MODEL$ClearType == "3"){"Fo"}
cat("\\n\\nRun prediction by sampling 50000 times for\\nKMR: ",  MODEL$KMR, "\\nClear Type: ", CT, "\\n\\n")

# Predictions
Pred <- predict_model3(model = MODEL, N = 50000, RandEff = "SA1ID")

# Save Predictions
output_name <- file.path(PRED_DIR , paste0("Pred" , str_extract(basename(MODEL_FPATH), "_.*")))
cat("\\n\\nSave Predictions to: ", output_name, "\\n")
qsave(Pred, output_name)}'

# write_lines(SLURM_templ, "SLURM/SLURM_template.job")

# Function to replace placeholders and generate .job files
GEN_R <- function(KMR, CLEARTYPE) {
  # Replace placeholders with actual values
  modified_text <- R_templ %>%
    str_replace_all("\\[KMR\\]", KMR) %>%
    str_replace_all("\\[CLEARTYPE\\]", CLEARTYPE) %>% 
    substr(2, nchar(.)-1)
  
  # Write the modified text to a new .job file
  output_name <- file.path("R", paste0("Pred_HPC_", KMR, "_", CLEARTYPE, ".R"))  # Change extension to .job
  write_lines(modified_text, output_name)

  # Optional: Print out the file name to check progress
  cat("Generated file:", output_name, "\n")
}

KMR <- c("CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC")
CLEARTYPE <- c("Ag", "In", "Fo")


walk(KMR, function(k) {
  walk(CLEARTYPE, function(c) {
    # Create R folder if it doesn't exist
    if(!dir.exists("R")) {
      dir.create("R")
    }
    GEN_R(k, c)
  })
})

