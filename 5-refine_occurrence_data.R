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
# output from 4-compile_occurrence_data.R
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
  "tidyverse","rnaturalearth","sp","tools","terra","textclean"
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
  
# set up file structure within your main working directory
data <- "occurrence_data"
standard <- "standardized_occurrence_data"
polygons <- "gis_layers"

################################################################################
# Read in data
################################################################################

# define projection
wgs_proj <- sp::CRS(SRS_string="EPSG:4326")
wgs_proj_terra <- "+proj=longlat +datum=WGS84"

# get urban areas layer and transform projection to WGS84
urban.poly <- rnaturalearth::ne_download(scale = "large", type = "urban_areas")
urban.poly <- spTransform(urban.poly,wgs_proj)

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
  .cen = "start", 
  .urb = "start",
  .inst = "start",
  .con = "start", 
  .outl = "start", 
  .gtsnative = "start", 
  .rlnative = "start",
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
    "taxon_name","scientificName","family","genus","specificEpithet",
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
    #additional optional data
    "taxon_name_status","iucnredlist_category",
    "natureserve_rank","fruit_nut",
    #flag columns
    "country_continent",#"country.name", "country.iso_a2", "country.iso_a3", 
    ".cen",".urb",
    ".inst",".con",".outl",".gtsnative",".rlnative",
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
    paste0(taxon, ".csv")))

  # make taxon points into a spatial object
  taxon_spdf <- vect(cbind(
    taxon_df$decimalLongitude,taxon_df$decimalLatitude),
    atts=taxon_df, crs=wgs_proj_terra)
          #taxon_spdf <- SpatialPointsDataFrame(taxon_df[,c("decimalLongitude",
          #"decimalLatitude")], taxon_df, proj4string = wgs_proj)
  # add country polygon data to each point based on lat-long location
  taxon_now <- intersect(taxon_spdf,world_polygons)
          #taxon_now <- point.in.poly(taxon_spdf, world_polygons, sp=TRUE)@data

  ## CHECK POINT LOCATION AGAINST "ACCEPTED" COUNTRY DISTRUBUTION FROM 2 SOURCES
  ## GlobalTreeSearch
  # species native country distribution list from GTS
  gts_list <- unique(unlist(strsplit(native_dist$gts_native_dist_iso2c[
    native_dist$taxon_name_accepted==taxon_nm], "; ")))
  if(!is.na(gts_list)){
  ## flag records where GTS country doesn't match record's coordinate location
  taxon_now <- taxon_now %>% mutate(.gtsnative=(ifelse(
    country.iso_a2 %in% gts_list, TRUE, FALSE)))
  } else {
    taxon_now$.gtsnative <- NA
  }
  ## IUCN Red List
  # species native country distribution list from RL
  rl_list <- unique(unlist(strsplit(native_dist$rl_native_dist_iso2c[
    native_dist$taxon_name_accepted==taxon_nm], "; ")))
  if(!is.na(rl_list)){
  ## flag records where RL country doesn't match record's coordinate location
  taxon_now <- taxon_now %>% mutate(.rlnative=(ifelse(
    country.iso_a2 %in% s.nd.rln.l, TRUE, FALSE)))
  } else {
    taxon_now$.rlnative <- NA
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
  taxon_now <- clean_coordinates(taxon_now,
    lon = "decimalLongitude", lat = "decimalLatitude",
    species = "taxon_name_accepted",
    centroids_rad = 500, # radius around capital coords (meters); default=1000
    inst_rad = 100, # radius around biodiversity institutions coords (meters)
    tests = c("centroids","institutions")
  )
  # adding urban area test separately because won't work when only 1 point
  if(nrow(taxon_df)<2){
    taxon_now2$.urb <- NA
    print("Speices with fewer than 2 records will not be tested.")
  } else {
    taxon_now <- as.data.frame(taxon_now2)
    flag_urb <- cc_urb(taxon_now2,
      lon = "decimalLongitude",lat = "decimalLatitude",
      ref = urban.poly, value = "flagged")
    taxon_now2$.urb <- flag_urb
  }
  # for some reason the "sea" flag isn't working in the above function...
  #    adding here separately
  # actually, found it flags a lot on islands, etc. Skipping for now.
  #   flag_sea <- cc_sea(taxon_now2,
  #     lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  #   taxon_now2$.sea <- flag_sea
  # for some reason the outlier section won't work when part of
  #   "clean_coordinates" function above so adding it here
  taxon_now2 <- as.data.frame(taxon_now2)
  flag_outl <- cc_outl(taxon_now2,
    lon = "decimalLongitude",lat = "decimalLatitude",
    species = "taxon_name_acc", method = "quantile",
    mltpl = 4, value = "flagged")
  taxon_now2$.outl <- flag_outl

  ## OTHER CHECKS
  ## Given country vs. lat-long country
  # check if given country matches lat-long country (CoordinateCleaner
  #   has something like this but also flags when NA? Didn't love that)
  taxon_now2 <- taxon_now2 %>% mutate(.con=(ifelse(
    (as.character(country.iso_a3) == as.character(countryCode_standard) &
    !is.na(country.iso_a3) & !is.na(countryCode_standard)) |
    is.na(country.iso_a3) | is.na(countryCode_standard), TRUE, FALSE)))
  ## Year
  taxon_now2 <- taxon_now2 %>% mutate(.yr1950=(ifelse(
    (as.numeric(year)>1950 | is.na(year)), TRUE, FALSE)))
  taxon_now2 <- taxon_now2 %>% mutate(.yr1980=(ifelse(
    (as.numeric(year)>1980 | is.na(year)), TRUE, FALSE)))
  taxon_now2 <- taxon_now2 %>% mutate(.yrna=(ifelse(
    !is.na(year), TRUE, FALSE)))

  # set column order and remove a few unnecessary columns
  taxon_now3 <- taxon_now2 %>% dplyr::select(all_of(col_names))
  # df of unflagged points
  unflagged <- taxon_now3 %>%
    filter(.cen & .urb &
       .inst & .con & .outl & .yr1950 & .yr1980 & .yrna &
      (.gtsnative | is.na(.gtsnative)) &
      (.rlnative  | is.na(.rlnative)) &
      basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "LIVING_SPECIMEN" &
      establishmentMeans != "INTRODUCED" & establishmentMeans != "MANAGED" &
      establishmentMeans != "INVASIVE")
  # add to summary table
  summary_add <- data.frame(
    taxon_name_acc = taxon_list[i],
    total_pts = nrow(taxon_now3),
    unflagged_pts = nrow(unflagged),
    .cen = sum(!taxon_now3$.cen),
    .urb = sum(!taxon_now3$.urb),
    .inst = sum(!taxon_now3$.inst),
    .con = sum(!taxon_now3$.con),
    .outl = sum(!taxon_now3$.outl),
    .gtsnative = sum(!taxon_now3$.gtsnative),
    .rlnative = sum(!taxon_now3$.rlnative),
    .yr1950 = sum(!taxon_now3$.yr1950),
    .yr1980 = sum(!taxon_now3$.yr1980),
    .yrna = sum(!taxon_now3$.yrna),
    stringsAsFactors=F)
  summary_tbl[i,] <- summary_add

  # WRITE NEW FILE
  write.csv(taxon_now3, file.path(main_dir, "occurrence_points","OUTPUTS_FROM_R","taxon_edited_points",
    paste0(taxon, ".csv")), row.names=FALSE)

  cat("Ending ", taxon, ", ", i, " of ", length(taxon_list), ".\n\n", sep="")
}

# write summary table
write.csv(summary_tbl, file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R",
  paste0("summary_of_output_points_", Sys.Date(), ".csv")),row.names = F)
