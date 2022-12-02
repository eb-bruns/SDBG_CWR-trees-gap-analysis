################################################################################

## 3-1_refine_occurrence_points_CWRcopy.R

### Authors: Emily Beckman Bruns & Shannon M. Still
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 26 May 2022
### Last updated: 03 June 2022

### R version 4.1.3

### DESCRIPTION:
  # Flags suspect points by adding a column for each type of flag, where
  #   FALSE = flagged. Most of the flagging is done through the
  #   'CoordinateCleaner' package, which was created for "geographic cleaning
  #   of coordinates from biologic collections."

### DATA IN:
  # output from 3_compile_raw_occurrence_points.R
  # tabular data:
  # - target_taxa_with_syn.csv
  # - globaltreesearch_country_distribution.csv
  # - spatialpolygon data ...
  #

### DATA OUT:
  # spp_edited_points folder with CSV of occurrence points for each target
  #   species (e.g., Quercus_lobata.csv)
  # Summary table with one row for each target species, listing number of
  #   points and number of flagged records in each flag column
  #   (flag_summary_by_sp.csv)

################################################################################
# Load libraries
################################################################################

my.packages <- c(#"raster",
  "sp", "tools","spatialEco",# "rgdal", "geosphere",
  #"readxl", "writexl", "dplyr", "tidyr",
  "tidyverse", #"housingData", "maps",
  "data.table",
  "textclean",
  "CoordinateCleaner", "countrycode",#, "usmap",
  #"RColorBrewer", "leaflet"
  "rnaturalearth","rnaturalearthhires"
)
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/occurrence_points"
#script_dir <- "./Documents/GitHub/OccurrencePoints/scripts"
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

# or use 0-1_set_workingdirectory.R script:
#source("./Documents/GitHub/OccurrencePoints/scripts/0-1_set_workingdirectory.R")
#source("scripts/0-1_set_workingdirectory.R")

################################################################################
##
## YOU ONLY NEED TO DO THIS FIRST PART ONCE -
# could move to another script but keeping together for now.
##
################################################################################
# Get native country-level occurrence data for target taxa
################################################################################

# read in taxa list
taxon_list_orig <- read.csv(file.path(main_dir,
  "target_taxa_with_synonyms.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
# keep only taxa with accepted species name
taxon_list_orig <- taxon_list_orig %>% filter(taxon_name_status=="Accepted")
  nrow(taxon_list_orig) #95

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"taxa_metadata")))
  dir.create(file.path(main_dir,"taxa_metadata"),
  recursive=T)

### GlobalTreeSearch (GTS)

# FIRST, download raw data
  # Go to https://tools.bgci.org/global_tree_search.php
  # Type your target genus name into the "Genus" box
  # Click "Search Plants" then scroll to the bottom and click "download as CSV
  #   file"
  # If you have more than one target genus, repeat the above steps for the
  #   other genera
  # Move all downloads to "inputs/known_distribution" folder
# read in and compile GlobalTreeSearch data
file_list <- list.files(path = file.path(main_dir,"taxa_metadata"),
  pattern = "globaltreesearch_results", full.names = T)
file_dfs <- lapply(file_list, read.csv, colClasses = "character",
  na.strings=c("","NA"),strip.white=T)
gts_list <- data.frame()
  for(file in seq_along(file_dfs)){
    gts_list <- rbind(gts_list, file_dfs[[file]])
  }
head(gts_list)
  # split countries by delimiter
gts_all <- gts_list %>%
  rename(taxon_name = taxon) %>%
  mutate(native_distribution =
    strsplit(as.character(native_distribution), "; ")) %>%
  unnest(native_distribution) %>% mutate(native_distribution =
    str_trim(native_distribution, side="both"))
# write out all GTS countries to check
#spp_countries <- as.data.frame(sort(unique(str_trim(
#  gts_all$native_distribution, side = c("both")))))
#write_xlsx(spp_countries, path=file.path(main_dir,"inputs",
#  "known_distribution","globaltreesearch_countries.xlsx"))

# use countrycode package to translate country codes from the country names
country_set <- as.data.frame(sort(unique(gts_all$native_distribution))) %>%
  add_column(iso3c = countrycode(sort(unique(gts_all$native_distribution)),
      origin="country.name", destination="iso3c")) %>%
  add_column(iso2c = countrycode(sort(unique(gts_all$native_distribution)),
      origin="country.name", destination="iso2c")) %>%
  add_column(iso3n = countrycode(sort(unique(gts_all$native_distribution)),
      origin="country.name", destination="iso3n")) %>%
  add_column(fips = countrycode(sort(unique(gts_all$native_distribution)),
      origin="country.name", destination="fips"))
names(country_set)[1] <- "country_name"
# add country codes to GTS native distribution data
gts_list$gts_native_dist_iso2c <- gts_list$native_distribution
gts_list$gts_native_dist_iso2c <- mgsub(gts_list$gts_native_dist_iso2c,
  array(as.character(country_set$country_name)),
  array(as.character(country_set$iso2c)))
names(gts_list)[4] <- "gts_native_dist"
head(gts_list)
# save the country codes for ISO2, ISO3, and numeric and character
#   codes, FIPS code
#  write_xlsx(country_set, path=file.path(main_dir,"inputs","gis_data",
#      "imls_global_admin_areas.xlsx"))
#gadm <- country_set; rm(country_set)

# add country codes to the taxon list by matching to GTS
  # match to accepted species names
taxon_list <- left_join(taxon_list_orig, gts_list[,c(2,4,5)],
  by=c("taxon_name_acc" = "taxon"))

# see which species have no GTS data
no_match <- taxon_list[which(is.na(taxon_list$gts_native_dist)),]$taxon_name_acc
no_match

### IUCN Red List (RL)

## To use data download manually from IUCN RL website:
  # Go to https://www.iucnredlist.org/search
  # Create an account if you don't have one (click "Login/Register" in top bar)
  # Login to your account
  # Open the "Taxonomy" tab in the left bar
  #   Either search for your target genus or just check "Plantae"
  #   You can limit the search further using the other tabs, if desired,
  #     but further refinement can sometimes exclude assessments you want
  #   When you're ready, on the right click "Download" then "Search Results"
  # You will receive an email when your download is ready
  # Next, go to your account (https://www.iucnredlist.org/account)
  #   Under "Saved downloads" click "Download" for your recent search
  #   Move downloaded folder to "occurrence_points/inputs/known_distribution"

# read in downloaded IUCN RL data for country-level species distribution
countries <- read.csv(file.path(main_dir,"taxa_metadata",
    "redlist_species_data","countries.csv"),
    colClasses = "character",na.strings=c("","NA"),strip.white=T)
# condense output so its one entry per species
countries_c <- countries %>%
  dplyr::filter(presence != "Extinct Post-1500") %>%
  rename(genus_species = scientificName) %>%
  dplyr::arrange(code) %>%
  dplyr::group_by(genus_species,origin) %>%
  dplyr::mutate(
    rl_native_dist_iso2c = paste(code, collapse = '; '),
    rl_native_dist = paste(name, collapse = '; ')) %>%
  dplyr::ungroup() %>%
  dplyr::select(genus_species,origin,rl_native_dist_iso2c,rl_native_dist) %>%
  dplyr::distinct(genus_species,origin,.keep_all=T)

# separate native dist countries from introduced dist countries
rl_native <- countries_c %>% filter(origin == "Native")
rl_introduced <- countries_c %>% filter(origin == "Introduced")
names(rl_introduced)[3] <- "rl_introduced_dist_iso2c"
names(rl_introduced)[4] <- "rl_introduced_dist"
# join both native and introduced together
rl_list <- full_join(rl_native[,c(1,4,3)],rl_introduced[,c(1,4,3)])

# add country codes to the taxon list by matching to RL data
taxon_list <- left_join(taxon_list, rl_list,
  by=c("taxon_name_acc" = "genus_species"))

# see which species have no RL data
no_match <- taxon_list[which(is.na(taxon_list$rl_native_dist)),]$taxon_name_acc
no_match

### Combine GTS and RL results

# keep only added native distribution columns
native_dist <- taxon_list %>%
  dplyr::select(taxon_name_acc,gts_native_dist,
    gts_native_dist_iso2c,rl_native_dist,rl_native_dist_iso2c,
    rl_introduced_dist,rl_introduced_dist_iso2c) %>%
  dplyr::distinct(taxon_name_acc,gts_native_dist,
    gts_native_dist_iso2c,rl_native_dist,rl_native_dist_iso2c,
    rl_introduced_dist,rl_introduced_dist_iso2c)

# see which target taxa have no distribution data matched
native_dist[is.na(native_dist$gts_native_dist) &
            is.na(native_dist$rl_native_dist),]$taxon_name_acc #35
  # *species* without country-level distribution data:
    #"Asimina incana"
    #"Asimina longifolia"
    #"Asimina pygmaea"
    #"Asimina reticulata"
    #"Carya carolinae-septentrionalis"
    #"Carya ovalis"
    #"Corylus californica"
    #"Pistacia texana"
    #"Prunus andersonii"
    #"Prunus fasciculata"
    #"Prunus geniculata"
    #"Prunus havardii"
    #"Prunus minutiflora"
    #"Prunus pumila"
    #"Prunus texana"

# create columns that combine GTS and RL
  # full names
native_dist$all_native_dist <- paste(native_dist$gts_native_dist,native_dist$rl_native_dist,sep="; ")
native_dist$all_native_dist <- str_squish(mgsub(native_dist$all_native_dist,
    c("NA; ","; NA","NA"),""))
native_dist$all_native_dist <- gsub(", ","~ ",native_dist$all_native_dist)
t <- setDT(native_dist)[,list(all_native_dist =
  toString(sort(unique(strsplit(all_native_dist,'; ')[[1]])))), by = taxon_name_acc]
native_dist <- native_dist %>% dplyr::select(-all_native_dist) %>% full_join(t)
native_dist$all_native_dist <- gsub(", ","; ",native_dist$all_native_dist)
native_dist$all_native_dist <- gsub("~ ",", ",native_dist$all_native_dist)
  # iso abb.
native_dist$all_native_dist_iso2 <- paste(native_dist$gts_native_dist_iso2c,
  native_dist$rl_native_dist_iso2c,sep="; ")
native_dist$all_native_dist_iso2 <- str_squish(mgsub(native_dist$all_native_dist_iso2,
    c("NA; ","; NA","NA"),""))
t <- setDT(native_dist)[,list(all_native_dist_iso2 =
  toString(sort(unique(strsplit(all_native_dist_iso2,'; ')[[1]])))), by = taxon_name_acc]
native_dist <- native_dist %>% dplyr::select(-all_native_dist_iso2) %>% full_join(t)
native_dist$all_native_dist_iso2 <- gsub(", ","; ",native_dist$all_native_dist_iso2)
head(native_dist)

# write taxon list with GTS and RL distribution information
write.csv(native_dist, file.path(main_dir,"taxa_metadata",
    "target_taxa_with_native_dist.csv"),row.names=F)
