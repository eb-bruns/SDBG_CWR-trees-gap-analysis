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
  "CoordinateCleaner", "countrycode"#, "usmap",
  #"RColorBrewer", "leaflet"
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
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/In situ - H - records"
poly_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"

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
taxon_list_orig <- read.csv(file.path(main_dir,"inputs",
  "target_taxa_with_syn.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
# keep only taxa with accepted species name
taxon_list_orig <- taxon_list_orig %>% filter(taxon_name_status=="Accepted")
  nrow(taxon_list_orig) #90

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","known_distribution")))
  dir.create(file.path(main_dir,"inputs","known_distribution"),
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
file_list <- list.files(path = file.path(main_dir,"inputs","known_distribution"),
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
countries <- read.csv(file.path(main_dir,"inputs","known_distribution",
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
write.csv(native_dist, file.path(main_dir,"inputs","known_distribution",
    "target_taxa_with_native_dist.csv"),row.names=F)


################################################################################
##
## NOW THE MAIN PART OF THE SCRIPT BEGINS
##
################################################################################
# 1. Read in data
################################################################################

# bring in polygon (load from saved .RData file)
load(file.path(poly_dir, "inputs", "gis_data", "admin_shapefiles.RData"))
# define projection
wgs_proj <- sp::CRS(SRS_string="EPSG:4326")
urban.poly <- spTransform(urban.poly,wgs_proj)

# read in country-level native distribution data
native_dist <- read.csv(file.path(main_dir,"inputs","known_distribution",
  "target_taxa_with_native_dist.csv"), header = T, na.strings = c("","NA"),
  colClasses = "character")

# create new folder for revised points, if not already existing
out.fld.nm <- "taxon_edited_points"
if(dir.exists(file.path(main_dir, "outputs", out.fld.nm)))
  print("directory already created") else dir.create(file.path(main_dir,
    "outputs", out.fld.nm), recursive=TRUE)

################################################################################
# 2. Iterate through species files and flag suspect points
################################################################################

# list of species files to iterate
all.spp.files <- list.files(path=file.path(main_dir, "outputs",
  "taxon_raw_points"), ignore.case=FALSE, full.names=FALSE, recursive=TRUE)
#all.spp.files <- all.spp.files[1:20]
spp_list <- file_path_sans_ext(all.spp.files)
#spp_list <- c("Magnolia_brasiliensis","Magnolia_boliviana","Magnolia_arcabucoana","Magnolia_angustioblonga")#"Magnolia_calimaensis")

# start a table to add summary of results for each species
summary_tbl <- data.frame(taxon_name_acc = "start", total_pts = "start",
  unflagged_pts = "start", .cen = "start", .urb = "start",
  .inst = "start",
  .con = "start", .outl = "start", .gtsnative = "start", .rlnative = "start",
  .yr1950 = "start", .yr1980 = "start",
  .yrna = "start", stringsAsFactors=F)

  # header/column name order and selection
  c.nms <- c("taxon_name_acc", "taxon_name", "scientificName",
    "taxonIdentificationNotes", "database", "all_source_databases", "year",
    "basisOfRecord", "establishmentMeans","decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "geolocationNotes", "localityDescription",
    "county", "stateProvince", "countryCode_standard", "institutionCode",
    "datasetName", "publisher", "rightsHolder", "license", "nativeDatabaseID",
    "references", "informationWithheld", "issue", "taxon_name_status", "UID",
    "country.name", "country.iso_a2", "country.iso_a3", "country.continent",
    ".cen",".urb",
    ".inst",".con",".outl",".gtsnative",".rlnative",
    ".yr1950",".yr1980",".yrna")

# iterate through each species file to flag suspect points
cat("Starting ", "target ", "taxa (", length(spp_list), " total)", ".\n\n",
  sep="")

for (i in 1:length(spp_list)){

  f.nm <- spp_list[i]

  # bring in records
  eo.df <- read.csv(file.path(main_dir, "outputs", "taxon_raw_points",
    paste0(f.nm, ".csv")))

  # create SpatialPointsDataFrame for species
  eo.spdf <- SpatialPointsDataFrame(eo.df[,c("decimalLongitude",
    "decimalLatitude")], eo.df, proj4string = wgs_proj)
  ## add country polygon data to each point based on lat-long location
  eo.post <- point.in.poly(eo.spdf, adm0.poly, sp=TRUE)@data

  ## CHECK POINT LOCATION AGAINST "ACCEPTED" COUNTRY DISTRUBUTION
  ## GlobalTreeSearch
  # species native country distribution list from GTS
  s.nd.gts.l <- unique(unlist(strsplit(native_dist$gts_native_dist_iso2c[
    native_dist$taxon_name_acc==gsub("_"," ", f.nm)], "; ")))
  if(!is.na(s.nd.gts.l)){
  ## flag records where GTS country doesn't match record's coordinate location
  eo.post <- eo.post %>% mutate(.gtsnative=(ifelse(
    country.iso_a2 %in% s.nd.gts.l, TRUE, FALSE)))
  } else {
    eo.post$.gtsnative <- NA
  }
  ## IUCN Red List native
  # species native country distribution list from RL
  s.nd.rln.l <- unique(unlist(strsplit(native_dist$rl_native_dist_iso2c[
    native_dist$taxon_name_acc==gsub("_"," ", f.nm)], "; ")))
  if(!is.na(s.nd.rln.l)){
  ## flag records where RL country doesn't match record's coordinate location
  eo.post <- eo.post %>% mutate(.rlnative=(ifelse(
    country.iso_a2 %in% s.nd.rln.l, TRUE, FALSE)))
  } else {
    eo.post$.rlnative <- NA
  }
  ## SERIES OF VETTED TESTS FROM CoordinateCleaner PACKAGE
  # Geographic Cleaning of Coordinates from Biologic Collections
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152
  #   Cleaning geographic coordinates by multiple empirical tests to flag
  #     potentially erroneous coordinates, addressing issues common in
  #     biological collection databases.
  ## tests included:
  # cc_cen -> Identify Coordinates in Vicinity of Country and Province Centroids
  # cc_inst -> Identify Records in the Vicinity of Biodiversity Institutions
  ## other test not included but could add:
  # cc_iucn -> Identify Records Outside Natural Ranges
  eo.post2 <- clean_coordinates(eo.post,
    lon = "decimalLongitude", lat = "decimalLatitude",
    species = "taxon_name_acc",
    centroids_rad = 500, # radius around capital coords (meters); default=1000
    inst_rad = 100, # radius around biodiversity institutions coord (meters)
    tests = c("centroids","institutions")
  )
  # adding urban area test separately because won't work for just 1 point
  if(nrow(eo.df)<2){
    eo.post2$.urb <- NA
    print("Speices with fewer than 2 records will not be tested.")
  } else {
    eo.post2 <- as.data.frame(eo.post2)
    flag_urb <- cc_urb(eo.post2,
      lon = "decimalLongitude",lat = "decimalLatitude",
      ref = urban.poly, value = "flagged")
    eo.post2$.urb <- flag_urb
  }
  # for some reason the "sea" flag isn't working in the above function...
  #    adding here separately
  # actually, found it flags a lot on islands, etc. Skipping for now.
  #   flag_sea <- cc_sea(eo.post2,
  #     lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  #   eo.post2$.sea <- flag_sea
  # for some reason the outlier section won't work when part of
  #   "clean_coordinates" function above so adding it here
  eo.post2 <- as.data.frame(eo.post2)
  flag_outl <- cc_outl(eo.post2,
    lon = "decimalLongitude",lat = "decimalLatitude",
    species = "taxon_name_acc", method = "quantile",
    mltpl = 4, value = "flagged")
  eo.post2$.outl <- flag_outl

  ## OTHER CHECKS
  ## Given country vs. lat-long country
  # check if given country matches lat-long country (CoordinateCleaner
  #   has something like this but also flags when NA? Didn't love that)
  eo.post2 <- eo.post2 %>% mutate(.con=(ifelse(
    (as.character(country.iso_a3) == as.character(countryCode_standard) &
    !is.na(country.iso_a3) & !is.na(countryCode_standard)) |
    is.na(country.iso_a3) | is.na(countryCode_standard), TRUE, FALSE)))
  ## Year
  eo.post2 <- eo.post2 %>% mutate(.yr1950=(ifelse(
    (as.numeric(year)>1950 | is.na(year)), TRUE, FALSE)))
  eo.post2 <- eo.post2 %>% mutate(.yr1980=(ifelse(
    (as.numeric(year)>1980 | is.na(year)), TRUE, FALSE)))
  eo.post2 <- eo.post2 %>% mutate(.yrna=(ifelse(
    !is.na(year), TRUE, FALSE)))

  # set column order and remove a few unnecessary columns
  eo.post3 <- eo.post2 %>% dplyr::select(all_of(c.nms))
  # df of unflagged points
  unflagged <- eo.post3 %>%
    filter(.cen & .urb &
       .inst & .con & .outl & .yr1950 & .yr1980 & .yrna &
      (.gtsnative | is.na(.gtsnative)) &
      (.rlnative  | is.na(.rlnative)) &
      basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "LIVING_SPECIMEN" &
      establishmentMeans != "INTRODUCED" & establishmentMeans != "MANAGED" &
      establishmentMeans != "INVASIVE")
  # add to summary table
  summary_add <- data.frame(
    taxon_name_acc = spp_list[i],
    total_pts = nrow(eo.post3),
    unflagged_pts = nrow(unflagged),
    .cen = sum(!eo.post3$.cen),
    .urb = sum(!eo.post3$.urb),
    .inst = sum(!eo.post3$.inst),
    .con = sum(!eo.post3$.con),
    .outl = sum(!eo.post3$.outl),
    .gtsnative = sum(!eo.post3$.gtsnative),
    .rlnative = sum(!eo.post3$.rlnative),
    .yr1950 = sum(!eo.post3$.yr1950),
    .yr1980 = sum(!eo.post3$.yr1980),
    .yrna = sum(!eo.post3$.yrna),
    stringsAsFactors=F)
  summary_tbl[i,] <- summary_add

  # WRITE NEW FILE
  write.csv(eo.post3, file.path(main_dir, "outputs", out.fld.nm,
    paste0(f.nm, ".csv")), row.names=FALSE)

  cat("Ending ", f.nm, ", ", i, " of ", length(spp_list), ".\n\n", sep="")
}

# write summary table
write.csv(summary_tbl, file.path(main_dir,"outputs",
  paste0("summary_of_output_points_", Sys.Date(), ".csv")),row.names = F)
