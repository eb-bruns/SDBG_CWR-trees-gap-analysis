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
urban.poly <- sf::read_sf(
  file.path(main_dir,gis_dir,"ne_10m_urban_areas/ne_10m_urban_areas.shp"))
urban.poly <- st_transform(urban.poly,wgs_proj)

# read in country-level native distribution data
native_dist <- read.csv(file.path(main_dir,"taxa_metadata",
  "target_taxa_with_native_dist.csv"), header = T, na.strings = c("","NA"),
  colClasses = "character")

# read in country boundaries shapefile
world_polygons <- vect(
  file.path(main_dir,polygons,
            "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))                             

# create new folder for revised points, if not already existing
if(!dir.exists(file.path(main_dir,data,standard,"taxon_edited_points")))
  dir.create(file.path(main_dir,data,standard,"taxon_edited_points"), 
    recursive=T)

################################################################################
# Iterate through taxon files and flag suspect points
################################################################################

# list of taxon files to iterate through
taxon_files <- list.files(path=file.path(main_dir,data,standard,
                                         "taxon_raw_points"), 
                          ignore.case=FALSE, full.names=FALSE, recursive=TRUE)
taxon_list <- file_path_sans_ext(taxon_files)

# start a table to add summary of results for each species
summary_tbl <- data.frame(
  taxon_name_accepted = "start", 
  total_pts = "start",
  unflagged_pts = "start", 
  selected_pts = "start", 
  .cen = "start", 
  .urb = "start",
  .inst = "start",
  .con = "start", 
  .outl = "start", 
  .nativectry = "start", 
  .yr1950 = "start", 
  .yr1980 = "start",
  .yrna = "start", 
    stringsAsFactors=F)

# select columns and order
  col_names <- c( 
    #data source and unique ID
    "UID","database","all_source_databases",
    #taxon
    "taxon_name_accepted",
    "taxon_name","scientificName","genus","specificEpithet",
    "taxonRank","infraspecificEpithet","taxonIdentificationNotes",
    #event
    "year","basisOfRecord",
    #record-level
    "nativeDatabaseID","institutionCode","datasetName","publisher",
    "rightsHolder","license","references","informationWithheld",
    "issue","recordedBy",
    #occurrence
    "establishmentMeans","individualCount",
    #location
    "decimalLatitude","decimalLongitude",
    "coordinateUncertaintyInMeters","geolocationNotes",
    "localityDescription","locality","verbatimLocality",
    "locationNotes","municipality","higherGeography","county",
    "stateProvince","country","countryCode","countryCode_standard",
    #additional optional data from target taxa list
    "taxon_name_status","iucnredlist_category",
    "natureserve_rank","fruit_nut",
    #flag columns
    "latlong_countryCode",
    ".cen",".urb",".inst",".con",".outl",
    ".nativectry",
    ".yr1950",".yr1980",".yrna"
  )

## iterate through each species file to flag suspect points
cat("Starting ","target ","taxa (", length(taxon_list)," total)",".\n\n",sep="")

for (i in 1:length(taxon_list)){

  taxon_file <- taxon_list[i]
  taxon_nm <- gsub("_", " ", taxon_file)
  taxon_nm <- mgsub(taxon_nm, c(" var "," subsp "), c(" var. "," subsp. "))

  # bring in records
  taxon_df <- read.csv(file.path(main_dir,data,standard,"taxon_raw_points",
    paste0(taxon_file, ".csv")))

  # make taxon points into a spatial object
  taxon_spdf <- vect(
    cbind(taxon_df$decimalLongitude,taxon_df$decimalLatitude),
    atts=taxon_df, crs="EPSG:4326")
  # add country polygon data to each point based on lat-long location
  taxon_now <- intersect(taxon_spdf,world_polygons)
  taxon_now <- as.data.frame(taxon_now)

  ## CHECK POINT LOCATION AGAINST "ACCEPTED" COUNTRY DISTRUBUTION
  #  this previously checked GlobalTreeSearch and IUCN Red List separately, but
  #     I have now combined into one - old code is commented out below
  ctry <- unique(unlist(strsplit(native_dist$all_native_dist_iso2[
    native_dist$taxon_name_accepted==taxon_nm], "; ")))
  if(!is.na(ctry[1])){
  ## flag records where RL country doesn't match record's coordinate location
    taxon_now <- taxon_now %>% mutate(.nativectry=(ifelse(
      ISO %in% ctry, TRUE, FALSE)))
  } else {
    taxon_now$.nativectry <- NA
  }
  ## GlobalTreeSearch: species native country distribution list from GTS
  #gts_ctry <- unique(unlist(strsplit(native_dist$gts_native_dist_iso2c[
  #  native_dist$taxon_name_accepted==taxon_nm], "; ")))
  #if(!is.na(gts_ctry[1])){
  ## flag records where GTS country doesn't match record's coordinate location
  #taxon_now <- taxon_now %>% mutate(.gtsnative=(ifelse(
  #  ISO %in% gts_ctry, TRUE, FALSE)))
  #} else {
  #  taxon_now$.gtsnative <- NA
  #}
  ## IUCN Red List: species native country distribution list from RL
  #rl_ctry <- unique(unlist(strsplit(native_dist$rl_native_dist_iso2c[
  #  native_dist$taxon_name_accepted==taxon_nm], "; ")))
  #if(!is.na(rl_ctry[1])){
  ## flag records where RL country doesn't match record's coordinate location
  #taxon_now <- taxon_now %>% mutate(.rlnative=(ifelse(
  #  ISO %in% rl_ctry, TRUE, FALSE)))
  #} else {
  #  taxon_now$.rlnative <- NA
  #}

  ## SERIES OF VETTED TESTS FROM CoordinateCleaner PACKAGE
  # Geographic Cleaning of Coordinates from Biologic Collections
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152
  #   Cleaning geographic coordinates by multiple empirical tests to flag
  #     potentially erroneous coordinates, addressing issues common in
  #     biological collection databases.
  # tests included:
  ## cc_cen - Identify Coordinates in Vicinity of Country and Province Centroids
  ## cc_inst - Identify Records in the Vicinity of Biodiversity Institutions
  taxon_now <- clean_coordinates(taxon_now,
    lon = "decimalLongitude", lat = "decimalLatitude",
    species = "taxon_name_accepted",
    centroids_rad = 500, # radius around capital coords (meters); default=1000
    inst_rad = 100, # radius around biodiversity institutions coords (meters)
    tests = c("centroids","institutions")
  )
  ## adding urban area test separately because won't work when only 1 point;
  # this is really slow when not using rnaturalearth version..
  if(nrow(taxon_df)<2){
    taxon_now$.urb <- NA
    print("Taxa with fewer than 2 records will not be tested.")
  } else {
    taxon_now <- as.data.frame(taxon_now)
    flag_urb <- cc_urb(taxon_now,
      lon = "decimalLongitude",lat = "decimalLatitude",
      ref = urban.poly, value = "flagged")
    taxon_now$.urb <- flag_urb
  }
  ## for some reason the "sea" flag isn't working in the clean_coordinates 
  #    function above; adding here separately
  # actually, found it flags a lot on islands, etc. Skipping for now.
  #flag_sea <- cc_sea(taxon_now2,
  #   lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  #taxon_now2$.sea <- flag_sea
  ## the outlier section won't work inside the clean_coordinates function either
  taxon_now <- as.data.frame(taxon_now)
  flag_outl <- cc_outl(taxon_now,
    lon = "decimalLongitude",lat = "decimalLatitude",
    species = "taxon_name_accepted", method = "quantile",
    mltpl = 4, value = "flagged")
  taxon_now$.outl <- flag_outl

  ## OTHER CHECKS
  ## Given country vs. lat-long country
  # check if given country matches lat-long country (CoordinateCleaner
  #   has something like this but also flags when NA? Didn't love that)
  taxon_now <- taxon_now %>% mutate(.con=(ifelse(
    (as.character(ISO) == as.character(countryCode_standard) &
    !is.na(ISO) & !is.na(countryCode_standard)) |
    is.na(ISO) | is.na(countryCode_standard), TRUE, FALSE)))
  ## Year
  taxon_now <- taxon_now %>% mutate(.yr1950=(ifelse(
    (as.numeric(year)>1950 | is.na(year)), TRUE, FALSE)))
  taxon_now <- taxon_now %>% mutate(.yr1980=(ifelse(
    (as.numeric(year)>1980 | is.na(year)), TRUE, FALSE)))
  taxon_now <- taxon_now %>% mutate(.yrna=(ifelse(
    !is.na(year), TRUE, FALSE)))

  # set column order and remove a few unnecessary columns
  taxon_now <- taxon_now %>% 
    rename("latlong_countryCode" = "ISO") %>%
    dplyr::select(all_of(col_names))
  # df of completely unflagged points
  total_unflagged <- taxon_now %>%
    filter(.cen & .urb &
             .inst & .con & .outl & .yr1950 & .yr1980 & .yrna &
             (.nativectry | is.na(.nativectry)) &
             basisOfRecord != "FOSSIL_SPECIMEN" & 
             basisOfRecord != "LIVING_SPECIMEN" &
             establishmentMeans != "INTRODUCED" & 
             establishmentMeans != "MANAGED" &
             establishmentMeans != "CULTIVATED" #& establishmentMeans != "INVASIVE"
    )
  # optional df of unflagged points based on selected filters you'd like to use
  select_unflagged <- taxon_now %>%
    filter(.cen & .inst & .outl & 
             #.urb & .con & .yr1950 & .yr1980 & .yrna &
             (.nativectry | is.na(.nativectry)) &
             basisOfRecord != "FOSSIL_SPECIMEN" & 
             basisOfRecord != "LIVING_SPECIMEN" &
             establishmentMeans != "INTRODUCED" & 
             establishmentMeans != "MANAGED" &
             establishmentMeans != "CULTIVATED" #& establishmentMeans != "INVASIVE"
    )
  # add to summary table
  summary_add <- data.frame(
    taxon_name_accepted = taxon_list[i],
    total_pts = nrow(taxon_now),
    unflagged_pts = nrow(total_unflagged),
    selected_pts = nrow(select_unflagged),
    .cen = sum(!taxon_now$.cen),
    .urb = sum(!taxon_now$.urb),
    .inst = sum(!taxon_now$.inst),
    .con = sum(!taxon_now$.con),
    .outl = sum(!taxon_now$.outl),
    .nativectry = sum(!taxon_now$.nativectry),
    .yr1950 = sum(!taxon_now$.yr1950),
    .yr1980 = sum(!taxon_now$.yr1980),
    .yrna = sum(!taxon_now$.yrna),
    stringsAsFactors=F)
  summary_tbl[i,] <- summary_add

  # WRITE NEW FILE
  write.csv(taxon_now, file.path(main_dir,data,standard,"taxon_edited_points",
    paste0(taxon_file, ".csv")), row.names=FALSE)

  cat("Ending ", taxon_nm, ", ", i, " of ", length(taxon_list), ".\n\n", sep="")
}

# write summary table
write.csv(summary_tbl, file.path(main_dir,data,
  paste0("summary_of_output_points_", Sys.Date(), ".csv")),row.names = F)
