################################################################################

## 5-refine_occurrence_data.R

### Authors: Emily Beckman Bruns & Shannon M Still
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: XX December 2022
### R version 4.2.2

### DESCRIPTION:
# Flags suspect points by adding a column for each type of flag, where
#   FALSE = flagged. Most of the flagging is done through the
#   'CoordinateCleaner' package, which was created for "geographic cleaning
#   of coordinates from biologic collections."

### INPUTS:
# output from 3_compile_raw_occurrence_points.R
# tabular data:
# - target_taxa_with_syn.csv
# - globaltreesearch_country_distribution.csv
# - spatialpolygon data ...
#

### OUTPUTS:
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
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"
  
# or use 0-set_working_directory.R script:
source("SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")
  
################################################################################
# Read in data
################################################################################

# define projection
wgs_proj <- sp::CRS(SRS_string="EPSG:4326")

# get urban areas layer and transform projection
urban.poly <- rnaturalearth::ne_download(scale = "large", type = "urban_areas")
urban.poly <- spTransform(urban.poly,wgs_proj)

# read in country-level native distribution data
native_dist <- read.csv(file.path(main_dir,"taxa_metadata",
  "target_taxa_with_native_dist.csv"), header = T, na.strings = c("","NA"),
  colClasses = "character")

# create new folder for revised points, if not already existing
out.fld.nm <- "taxon_edited_points"
if(dir.exists(file.path(main_dir, "occurrence_points","OUTPUTS_FROM_R", out.fld.nm)))
  print("directory already created") else dir.create(file.path(main_dir,
     "occurrence_points","OUTPUTS_FROM_R", out.fld.nm), recursive=TRUE)

################################################################################
# 2. Iterate through species files and flag suspect points
################################################################################

# list of species files to iterate
all.spp.files <- list.files(path=file.path(main_dir, "occurrence_points","OUTPUTS_FROM_R",
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
  eo.df <- read.csv(file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R","taxon_raw_points",
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
  write.csv(eo.post3, file.path(main_dir, "occurrence_points","OUTPUTS_FROM_R",out.fld.nm,
    paste0(f.nm, ".csv")), row.names=FALSE)

  cat("Ending ", f.nm, ", ", i, " of ", length(spp_list), ".\n\n", sep="")
}

# write summary table
write.csv(summary_tbl, file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R",
  paste0("summary_of_output_points_", Sys.Date(), ".csv")),row.names = F)
