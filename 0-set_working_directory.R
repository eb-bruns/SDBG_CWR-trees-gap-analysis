################################################################################

### 0-set_working_directory.R

### Authors: Shannon M Still & Emily Beckman Bruns
### Funding: Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).

### Creation date: 21 May 2020
### Last updated: 02 December 2022
### R version 4.2.2

### DESCRIPTION:
# This script sets the working environment for the computer on which you are
# working. You can add your own computer and working directories by filling in 
# the blanks (________) below!

################################################################################
# Set working environment depending on your computer
################################################################################

# Use this to check your "nodename"
# Sys.info()[4]

## For Emily Beckman Bruns:
if (Sys.info()[4] == "Africa.local") {
  # set main working directory
    # for IMLS Collections Value project
  #main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"
    # for North American fruit and nut tree CWR project
  main_dir <- "/Users/emily/Library/CloudStorage/GoogleDrive-ebeckman@mortonarb.org/Shared drives/Global Tree Conservation Program/4. GTCP_Projects/Gap Analyses/Conservation Gap Analysis - FRUIT & NUT TREE CWR NORTH AMERICA/Gap_Analysis"
  # set folder structure within main working directory
  gis_dir <- "gis_layers"
  taxa_dir <- "target_taxa"
  exsitu_dir <- "exsitu_data"
  occ_dir <- "occurrence_data"
  # cloned gap analysis repository
  gap_dir <- "/Users/emily/Documents/GitHub/GapAnalysis/R"
  # OPTIONAL: set local working directory, for trialing locally before saving
  #   to main working directory
  local_dir <- "./Desktop/*work/NorthAm-CWR"
  # set location for login information (e.g., for GBIF)
  log_loc <- file.path(local_dir, "IMLS_passwords.txt")
  # prints computer name, to let you know you're in the right spot
  print(paste("Working from the lovely", Sys.info()[4]))

## For additional user or workstation (fill in ________ with your info):
} else if (Sys.info()[4] == "________") {
  # set main working directory
  main_dir <- "________"
  # set folder structure within main working directory
  gis_dir <- "gis_layers"
  taxa_dir <- "target_taxa"
  exsitu_dir <- "exsitu_data"
  occ_dir <- "occurrence_data"
  # OPTIONAL: set local working directory, for trialing locally before saving
  #   to main working directory
  local_dir <- "________"
  # set location for login information (e.g., for GBIF; see 3-get_occurrence_data.R)
  log_loc <- file.path(local_dir, "________.txt")
  # prints computer name, to let you know you're in the right spot
  print(paste("Working from the lovely", Sys.info()[4]))

## can add as many additional "else if" sections as needed to cover other
#   workstations

} else {
  # default, which sets the working driectory as the folder from which you
  #   opened the scripts/project
  setwd(getwd())
  print("You should add your info to the 0-1_set_working_directory.R script so
    this line automatically sets up all the working directories you'll be using!")
}
