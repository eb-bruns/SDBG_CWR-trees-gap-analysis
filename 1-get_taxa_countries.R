################################################################################

## 1-get_taxa_countries.R

### Authors: Emily Beckman Bruns & Shannon M Still
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Minor edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: 02 December 2022
### R version 4.2.2

### DESCRIPTION:
# Gets IUCN Red List and BGCI GlobalTreeSearch native countries of occurrence
# for taxa in your target taxa list. These are later used to filter points.

### INPUTS:
# target_taxa_with_synonyms.csv

### OUTPUTS:
#

################################################################################
# Load libraries
################################################################################

my.packages <- c("tidyverse","countrycode","textclean","data.table")
#install.packages(my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"
  
# or use 0-set_working_directory.R script:
source("SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")
  
################################################################################
# Read in target taxa and create folder for outputs
################################################################################

# read in taxa list
taxon_list <- read.csv(file.path(main_dir,"target_taxa_with_synonyms.csv"), 
                       header = T,na.strings=c("","NA"),colClasses="character")
# keep only accepted taxa
taxon_list <- taxon_list %>% 
  filter(taxon_name_status=="Accepted")
nrow(taxon_list) #95

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"taxa_metadata")))
  dir.create(file.path(main_dir,"taxa_metadata"),recursive=T)

################################################################################
# Get BGCI GlobalTreeSearch (GTS) native country data
################################################################################

# FIRST, download raw data
  # Go to https://tools.bgci.org/global_tree_search.php
  # Type your target genus name into the "Genus" box
  # Click "Search" then scroll to the bottom and click "Download as CSV file"
  # If you have more than one target genus, repeat the above steps for the
  #   other genera
  # Move all downloads to your "taxa_metadata" folder

# read in and compile GTS data
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

# can check countries list
spp_countries <- as.data.frame(sort(unique(str_trim(
  gts_all$native_distribution, side = c("both")))))
spp_countries

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

# add country codes to the taxon list by matching to GTS
# match to accepted taxon names
taxon_list <- left_join(taxon_list, gts_list[,c(2,4,5)],
  by=c("taxon_name_accepted" = "taxon"))
head(taxon_list)

# see which taxa have no GTS data
no_match <- taxon_list[which(is.na(
  taxon_list$gts_native_dist)),]$taxon_name_acc
no_match

################################################################################
# Get IUCN Red List (RL) native country data
################################################################################

# FIRST, download raw data
  # Go to https://www.iucnredlist.org/search
  # Create an account if you don't have one (click "Login/Register" in top bar)
  # Log in to your account
  # Open the "Taxonomy" tab in the left bar
  #   Either search for your target genus or just check "Plantae" to download
  #   data for all plant species that have assessments globally.
  #   You can also search by Family/families if preferred.
  #   You can limit the search further using the other tabs, as desired,
  #   but further refinement can sometimes exclude assessments you want.
  #   When you're ready, on the right click "Download" then "Search Results"
  # You will receive an email when your download is ready
  # Next, go to your account (https://www.iucnredlist.org/account)
  #   Under "Saved downloads" click "Download" for your recent search
  #   Move downloaded folder to your "taxa_metadata" folder and rename it
  #   simply "redlist_species_data"

# read in downloaded RL data for country-level species distribution
countries <- read.csv(file.path(main_dir,"taxa_metadata",
    "redlist_species_data","countries.csv"),
    colClasses = "character",na.strings=c("","NA"),strip.white=T)

# condense output so its one entry per species
countries_c <- countries %>%
  dplyr::filter(presence != "Extinct Post-1500") %>%
  rename(taxon_name_accepted = scientificName) %>%
  dplyr::arrange(code) %>%
  dplyr::group_by(taxon_name_accepted,origin) %>%
  dplyr::mutate(
    rl_native_dist_iso2c = paste(code, collapse = '; '),
    rl_native_dist = paste(name, collapse = '; ')) %>%
  dplyr::ungroup() %>%
  dplyr::select(taxon_name_accepted,origin,rl_native_dist_iso2c,
                rl_native_dist) %>%
  dplyr::distinct(taxon_name_accepted,origin,.keep_all=T)

# separate native dist countries from introduced dist countries
rl_native <- countries_c %>% filter(origin == "Native")
rl_introduced <- countries_c %>% filter(origin == "Introduced")
names(rl_introduced)[3] <- "rl_introduced_dist_iso2c"
names(rl_introduced)[4] <- "rl_introduced_dist"
# join both native and introduced together
rl_list <- full_join(rl_native[,c(1,4,3)],rl_introduced[,c(1,4,3)])

# add country codes to the taxon list by matching to RL data
taxon_list <- left_join(taxon_list, rl_list, by="taxon_name_accepted")

# see which species have no RL data
no_match <- taxon_list[which(is.na(taxon_list$rl_native_dist)),]$taxon_name_acc
no_match

################################################################################
# Combine GTS and RL results
################################################################################

# keep only added native distribution columns
native_dist <- taxon_list %>%
  dplyr::select(taxon_name_accepted,gts_native_dist,
    gts_native_dist_iso2c,rl_native_dist,rl_native_dist_iso2c,
    rl_introduced_dist,rl_introduced_dist_iso2c) %>%
  dplyr::distinct(taxon_name_accepted,gts_native_dist,
    gts_native_dist_iso2c,rl_native_dist,rl_native_dist_iso2c,
    rl_introduced_dist,rl_introduced_dist_iso2c)

# see which target taxa have no distribution data matched
native_dist[is.na(native_dist$gts_native_dist) &
            is.na(native_dist$rl_native_dist),]$taxon_name_acc #40
  # *species* without country-level distribution data:
    #"Asimina incana"
    #"Asimina longifolia"
    #"Asimina manasota"
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

# create columns that combine GTS and RL country data
  # full country names
native_dist$all_native_dist <- paste(
  native_dist$gts_native_dist,native_dist$rl_native_dist,sep="; ")
native_dist$all_native_dist <- str_squish(mgsub(native_dist$all_native_dist,
    c("NA; ","; NA","NA"),""))
native_dist$all_native_dist <- gsub(", ","~ ",native_dist$all_native_dist)
t <- setDT(native_dist)[,list(all_native_dist =
  toString(sort(unique(strsplit(all_native_dist,'; ')[[1]])))), 
  by = taxon_name_accepted]
native_dist <- native_dist %>% dplyr::select(-all_native_dist) %>% full_join(t)
native_dist$all_native_dist <- gsub(", ","; ",native_dist$all_native_dist)
native_dist$all_native_dist <- gsub("~ ",", ",native_dist$all_native_dist)
  # ISO country code abbreviations
native_dist$all_native_dist_iso2 <- paste(native_dist$gts_native_dist_iso2c,
  native_dist$rl_native_dist_iso2c,sep="; ")
native_dist$all_native_dist_iso2 <- str_squish(mgsub(
  native_dist$all_native_dist_iso2, c("NA; ","; NA","NA"),""))
t <- setDT(native_dist)[,list(all_native_dist_iso2 =
  toString(sort(unique(strsplit(all_native_dist_iso2,'; ')[[1]])))), 
  by = taxon_name_accepted]
native_dist <- native_dist %>% 
  dplyr::select(-all_native_dist_iso2) %>% 
  full_join(t)
native_dist$all_native_dist_iso2 <- gsub(
  ", ","; ",native_dist$all_native_dist_iso2)

################################################################################
# *Optionally* Add additional native countries manually
################################################################################

# adding United States to all species, since that is our target region;
# can edit or skip depending on your needs
  # rows with no countries:
native_dist[which(native_dist$all_native_dist == ""),
            "all_native_dist"] <- "United States"
native_dist[which(native_dist$all_native_dist == ""),
            "all_native_dist_iso2"] <- "US"
  # rows with countries but missing the US:
native_dist[which(!grepl("United States",native_dist$all_native_dist)),
            "all_native_dist"] <- "United States"
native_dist[which(!grepl("United States",native_dist$all_native_dist)),
            "all_native_dist_iso2"] <- "US"

# check it out
head(native_dist)
unique(native_dist$all_native_dist)

################################################################################
# Write file
################################################################################

write.csv(native_dist, file.path(main_dir,"taxa_metadata",
    "target_taxa_with_native_dist.csv"),row.names=F)
