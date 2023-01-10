################################################################################

## 5-refine_occurrence_data.R

### Authors: Emily Beckman Bruns & Shannon M Still
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: 09 December 2022
### R version 4.2.2

### DESCRIPTION:
# Flags suspect points by adding a column for each type of flag, where
#   FALSE = flagged. Most of the flagging is done through the
#   'CoordinateCleaner' package, which was created for "geographic cleaning
#   of coordinates from biologic collections."

### INPUTS:
# target_taxa_with_synonyms.csv
# target_taxa_with_native_dist.csv (output from 1-get_taxa_countries.R)
# outputs from 4-compile_occurrence_data.R
# polygons ...

### OUTPUTS:
# folder (spp_edited_points) with CSV of occurrence points for each target
#   taxon (e.g., Quercus_lobata.csv)
# Summary table with one row for each target taxon, listing number of
#   points and number of flagged records in each flag column
#   (flag_summary_by_sp.csv)

################################################################################
# Load libraries
################################################################################

my.packages <- c(
  "tidyverse","rnaturalearth","sp","tools","terra","textclean",
  "CoordinateCleaner","sf"
)
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

# or use 0-set_working_directory.R script:
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")

# set up file structure within your main working directory
data <- "occurrence_data"
standard <- "standardized_occurrence_data"
polygons <- "gis_layers"

################################################################################
# Read in data
################################################################################

# define projection
wgs_proj <- sp::CRS(SRS_string="EPSG:4326")
# get urban areas layer and transform projection to WGS84
# the ne_download function isn't working now... doing manually 9 Jan 2023
#urban.poly <- rnaturalearth::ne_download(scale = "large", type = "urban_areas")
### !! THIS NEEDS TO BE TESTED STILL !!!
urban.poly <- sf:read_st(
  file.path(main_dir,gis_dir,"ne_10m_urban_areas/ne_10m_urban_areas.shp"))
