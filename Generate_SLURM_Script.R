# Purpose this script is for estimating HPC resource and generating SLURM script for HPC job scheduling

library(qs)
library(fs)
library(tidyverse)
library(tools)

## Define directory
INPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/input")
OUTPUT_DIR <- file.path("D:/Data/NSW_Deforestation/risk-model/output")

# Check the HPC SACCT report ----
## If you have run the job before, you can estimate the resource usage based on the sacct report of the past jobs
## First generate the HPC SACCT report in .CSV format
## sacct -a --starttime now-3days --format JobID,JobName,State,AllocCPUS,REQMEM,TotalCPU,Elapsed,MaxRSS,ExitCode,NNodes,NodeList,NTasks -u $USER -P --delimiter="," > Pd_In_SACCT.csv

# Read in HPC SACCT report in .CSV format
Pd_SACCT <- read.csv("SLURM/Pd_SACCT.csv") %>% 
  filter(State == "COMPLETED")%>%
  # remove "k" from MaxRSS and convert to numeric
  mutate(MaxRSS = round(as.numeric(gsub("K", "", MaxRSS))/1024^2,2),
         TotalCPU = hms(TotalCPU),
         Elapsed = hms(Elapsed),
         JobID = sub("^(\\d{8}).*", "\\1", JobID)) %>% 
  group_by(JobID) %>% 
  summarise(MaxRSS = sum(MaxRSS, na.rm = TRUE)) %>% 
  left_join(read.csv("SLURM/Pd_SACCT.csv") %>% select(JobID, JobName, Elapsed), by = "JobID") %>% 
  select(JobID, JobName, MaxRSS, Elapsed) %>% 
  mutate(matches = str_match(JobName, "^Pd_([^_]+)_([^_]+)$"),
         KMR = matches[, 2],
         CT = matches[, 3],
         Elapsed = as.numeric(hms(Elapsed), units = "seconds")) %>% 
  select(-matches)


# Set the directory
Model_Flist <- list.files(file.path(OUTPUT_DIR, "models"), pattern = glob2rx("Model_*.qs"), full.names = TRUE)

# List files and get their sizes
Model_Finfo <- fs::file_info(Model_Flist) %>% 
  mutate(FNAME = basename(file_path_sans_ext(path)),
         FSIZE_MB = round(as.numeric(size) / 1024^2, 2),
         matches = str_match(FNAME, "^Model_([^_]+)_([^_]+)$"),
         KMR = matches[, 2],
         CT = matches[, 3]) %>%
  dplyr::select(FNAME, FSIZE_MB, KMR, CT)

# Estimate CPU time and memory usage for HPC job scheduling----
## If you have not run the job before, you can estimate the resource usage based on the model file size
### MODEL_HPC <- ...

## if you have run the job before, you can estimate the resource usage based on the sacct report of the past jobs
MODEL_HPC <- Pd_SACCT %>% 
  mutate(MEMORY = round(MaxRSS*1.2, -1),
         TIME = Elapsed * 1.5,
         H = sprintf("%02d", TIME  %/% 3600),
         M = sprintf("%02d", (TIME  %% 3600) %/% 60),
         TIME_HMS = paste(H, M, "00", sep = ":"),
         JOB_NAME = paste0("Pd_", KMR, "_", CT),
         RSCRIPT = paste0("Pred_HPC_", KMR, "_", CT, ".R"),
         SLURM_FN = paste0("Pred_HPC_", KMR, "_", CT, ".job")) %>% 
  select(MEMORY, JOB_NAME, TIME = TIME_HMS, RSCRIPT, SLURM_FN) %>% 
  mutate(across(everything(), as.character))
  
# Step 1: Define the template ----
SLURM_templ <- '{#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=[MEMORY]G
#SBATCH --job-name=[JOB_NAME]
#SBATCH --time=[TIME]
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --batch=epyc3
#SBATCH --account=a_rhodes
#SBATCH --output=HPC_out/[JOB_NAME].out
#SBATCH --error=HPC_out/[JOB_NAME].error

#--- 
echo "Details for this job:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}"
echo "- SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "- SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE}"
echo "- SLURM_JOB_NAME=${SLURM_JOB_NAME}"
echo "- SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}"

echo -e "/n/n"

#--- 
# load the software module
module load r/4.4.0-combo-EPYC3-only

# Check the path or Rscript
which Rscript

# Print the working directory
pwd
ls

#example job using the above variables
srun Rscript $SLURM_SUBMIT_DIR/../R/[RSCRIPT]}'

# write_lines(SLURM_templ, "SLURM/SLURM_template.job")

# Function to replace placeholders and generate .job files
GEN_SLURM <- function(MEMORY, JOB_NAME, TIME, RSCRIPT, SLURM_FN) {
  # Replace placeholders with actual values
  modified_text <- SLURM_templ %>%
    str_replace_all("\\[MEMORY\\]", MEMORY) %>%
    str_replace_all("\\[JOB_NAME\\]", JOB_NAME) %>%
    str_replace_all("\\[TIME\\]", TIME) %>%
    str_replace_all("\\[RSCRIPT\\]", RSCRIPT) %>% 
    substr(2, nchar(.)-1)
  
  # Write the modified text to a new .job file
  output_name <- file.path("SLURM", SLURM_FN)  # Change extension to .job
  write_lines(modified_text, output_name)

  # Optional: Print out the file name to check progress
  cat("Generated file:", output_name, "\n")
}

# Create SLURM folder if it doesn't exist
if(!dir.exists("SLURM")) {
  dir.create("SLURM")
}
pwalk(MODEL_HPC, GEN_SLURM)

# UNIX command to submit jobs to SLURM scheduler
paste0("sbatch ",list.files("SLURM/", pattern = glob2rx("Pred*_Ag.job")))
paste0("sbatch ",list.files("SLURM/", pattern = glob2rx("Pred*_Fo.job")))
paste0("sbatch ",list.files("SLURM/", pattern = glob2rx("Pred*_In.job")))
 
### Final SLURM Scripts need to be modified for better use of resources 