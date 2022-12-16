################################################################################

## 4-compile_occurrence_data.R

### Author: Emily Beckman Bruns
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: 09 December 2022
### R version 4.2.2

### DESCRIPTION:
# This script compiles in situ occurrence point data and ex situ (genebank
# and botanical garden) data. We remove any rows for species not in our
# target taxa list and standardize some key columns. Records without 
# coordinates are separated out and saved, then records with coordinates are 
# saved in taxon-specific files.

### INPUTS:
# target_taxa_with_synonyms.csv
# occurrence data downloaded in 2-get_occurrence_data.R

### OUTPUTS:
#
# folder (taxon_raw_points) with CSV of geolocated (have lat-long) data for 
# each target taxon (e.g., Malus_angustifolia.csv)
#
# CSV of all occurrence points without lat-long but with locality description
#   (need_geolocation.csv)
#
# summary table (occurrence_record_summary.csv) with one row for each
# target taxon, listing the number of records with valid lat-long, and 
# number of records with locality description only
#

################################################################################
# Load libraries
################################################################################

my.packages <- c(
  'tidyverse','CoordinateCleaner','tidyterra','terra','raster','textclean',
  'countrycode','data.table'
)
# install.packages (my.packages) #Turn on to install current versions
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
# Read in and compile occurrence point data
################################################################################

# create foldes for output data
if(!dir.exists(file.path(main_dir,data)))
  dir.create(file.path(main_dir,data), recursive=T)
if(!dir.exists(file.path(main_dir,data,standard,"taxon_raw_points")))
  dir.create(file.path(main_dir,data,standard,"taxon_raw_points"), recursive=T)

# read in raw occurrence data
file_list <- list.files(file.path(main_dir,data,standard), pattern = ".csv", 
                        full.names = T)
file_dfs <- lapply(file_list, read.csv, header = T, na.strings = c("","NA"),
                   colClasses = "character")
length(file_dfs) #7

# stack all datasets using bind_rows, which keeps non-matching columns
#   and fills with NA; 'Reduce' iterates through list and merges with previous.
# this may take a few minutes if you have lots of data
all_data_raw <- Reduce(bind_rows, file_dfs)
nrow(all_data_raw) #2557933
names(all_data_raw) #37
table(all_data_raw$database)
# BIEN    Ex_situ  FIA     GBIF    iDigBio  IUCN_RedList  NorthAm_herbaria
# 871544  9318     922368  481187  38256    52322         182134 

# add unique identifier
nms <- names(all_data_raw)
all_data_raw <- all_data_raw %>% 
  dplyr::mutate(UID=paste0('id', sprintf("%08d",1:nrow(all_data_raw)))) %>% 
  dplyr::select(c('UID', all_of(nms)))
rm(nms, file_dfs, file_list)

all_data <- all_data_raw
rm(all_data_raw)

################################################################################
# Filter by target taxa
################################################################################

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "target_taxa_with_synonyms.csv"),
                       header = T, na.strings = c("","NA"),colClasses = "character")
target_taxa <- unique(taxon_list$taxon_name)
# if needed, add columns that separate out taxon name
taxon_list <- taxon_list %>%
  separate("taxon_name",c("genus","species","infra_rank","infra_name"),
           sep=" ",remove=F,fill="right")

# join data to taxon list
all_data <- all_data %>% dplyr::select(-genus)
all_data <- left_join(all_data,taxon_list)
# join again just by species name if no taxon match
need_match <- all_data[which(is.na(all_data$taxon_name_status)),]
nrow(need_match) #59933
# remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data)-ncol(taxon_list)+1)]
# rename column for matching
need_match <- need_match %>% rename(taxon_name_full = taxon_name)
need_match$taxon_name <- need_match$species_name
# new join
need_match <- left_join(need_match,taxon_list)
# bind together new matches and previously matched
matched <- all_data[which(!is.na(all_data$taxon_name_status)),]
matched$taxon_name_full <- matched$taxon_name
all_data <- rbind(matched,need_match)
table(all_data$taxon_name_status) # Accepted: 2252608   Synonym: 254542
# fill in extra data for synonyms
all_data <- all_data %>% 
  dplyr::select(-iucnredlist_category,-natureserve_rank,-fruit_nut)
taxon_list_add <- taxon_list %>% 
  filter(taxon_name_status == "Accepted") %>%
  dplyr::select(taxon_name_accepted,iucnredlist_category,natureserve_rank,fruit_nut)
all_data <- left_join(all_data,taxon_list_add)

# check names that got excluded.....
still_no_match <- all_data[which(is.na(all_data$taxon_name_status)),]
nrow(still_no_match) #49979
sort(unique(still_no_match$taxon_name))
table(still_no_match$database)
rm(matched,need_match,taxon_list_add,still_no_match)

# keep only rows for target taxa
all_data <- all_data %>% 
  filter(!is.na(taxon_name_status) & !is.na(UID))
nrow(all_data) #2512150

### ! target taxa with no occurrence data:
unique(taxon_list$taxon_name_acc)[
  !(unique(taxon_list$taxon_name_acc) %in% (unique(all_data$taxon_name_acc)))]
# none ! wow !

################################################################################
# Standardize/check some key columns
################################################################################

# create localityDescription column
# "NAs introduced by coercion" warning ok
all_data <- all_data %>%
  unite("localityDescription",
        c(locality,municipality,higherGeography,county,stateProvince,country,
          countryCode,locationNotes,verbatimLocality), remove = F, 
        sep = " | ") %>%
  mutate(decimalLatitude=as.numeric(decimalLatitude),
         decimalLongitude=as.numeric(decimalLongitude))
# get rid of NAs but keep pipes, so you can split back into parts if desired
all_data$localityDescription <- mgsub(all_data$localityDescription,
                                      c("NA "," NA"), "")
# if no locality info at all, make it NA
all_data$localityDescription <- gsub("| | | | | | | |", NA,
                                     all_data$localityDescription, fixed = T)
# check it
head(unique(all_data$localityDescription))

# check year column
# ex situ years were not standardized so most of those get removed here :(
# "NAs introduced by coercion" warning ok
all_data$year <- as.numeric(all_data$year)
# remove values less than 1500 or greater than current year
all_data$year[which(all_data$year < 1500 |
                      all_data$year > as.numeric(format(Sys.time(),"%Y")))] <- NA
sort(unique(all_data$year))

# check basis of record column
unique(all_data$basisOfRecord)
all_data$basisOfRecord[which(is.na(all_data$basisOfRecord))] <- "UNKNOWN"

# check establishment means
unique(all_data$establishmentMeans)
# use this if you need to recode any mistakes:
#all_data <- all_data %>%
#  mutate(establishmentMeans = recode(establishmentMeans,
#    "Introduced" = "INTRODUCED",
#    "Uncertain" = "UNCERTAIN",
#    "Native" = "NATIVE"))
all_data$establishmentMeans[which(is.na(all_data$establishmentMeans))] <-
  "UNCERTAIN"

# check validity of lat and long; remove invalid coordinates
# if coords are both 0, set to NA
zero <- which(all_data$decimalLatitude == 0 & all_data$decimalLongitude == 0)
all_data$decimalLatitude[zero] <- NA; all_data$decimalLongitude[zero] <- NA
rm(zero)
# flag non-numeric and not available coordinates and lat > 90, lat < -90,
# lon > 180, and lon < -180
coord_test <- cc_val(all_data, lon = "decimalLongitude",lat = "decimalLatitude",
                     value = "flagged", verbose = TRUE) #Flagged 317852 records.
# try switching lat and long for invalid points and check validity again
all_data[!coord_test,c("decimalLatitude","decimalLongitude")] <-
  all_data[!coord_test,c("decimalLongitude","decimalLatitude")]
coord_test <- cc_val(all_data, lon = "decimalLongitude",lat = "decimalLatitude",
                     value = "flagged", verbose = TRUE) #Flagged 317852 records.
# mark these as flagged
all_data$flag <- NA
all_data[!coord_test,]$flag <- "Coordinates invalid"
rm(coord_test)
 
################################################################################
# Select records that have coordinates and are on land
################################################################################

# move forward with subset of points that do have lat and long
geo_pts <- all_data %>% dplyr::filter(is.na(flag))
nrow(geo_pts) #2194298
no_geo_pts <- all_data %>% dplyr::filter(!is.na(flag))

# check if points are in water and mark
# you can also skip this if you don't mind the water points - can take 
# a few minutes or more to buffer and intersect the polygons
  # read in world polygons shapefile
world_polygons <- vect(file.path(main_dir,polygons,
   "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))                             
  # add 1km buffer to ecoregion layer, so we keep pts close to land;
  # width is in meters - change as desired; > ~500 threw this error:
    # "IllegalArgumentException: point array must contain 0 or >1 elements"
world_buff <- terra::buffer(world_polygons, width=500)
rm(world_polygons)
  # make occurrence points a spatial object
geo_pts_spatial <- vect(
  cbind(geo_pts$decimalLongitude,geo_pts$decimalLatitude),
  atts=geo_pts, crs="EPSG:4326")
  # intersect points with ecoregions layer
land_pts <- terra::intersect(geo_pts_spatial,world_buff)
on_land <- as.data.frame(land_pts); nrow(on_land) #2185947
land_id <- unique(on_land$UID)
  # flag points not on land (in water)
geo_pts[which(!(geo_pts$UID %in% land_id)),"flag"] <- "Coordinates in water"
all_data <- rbind(no_geo_pts,geo_pts)
table(all_data$flag)
  # Coordinates in water  Coordinates invalid 
  # 10249                 317852

# separate out points with locality description only, including:
# lat-long is invalid or in water
locality_pts <- all_data %>% 
  dplyr::filter(!is.na(localityDescription) & !is.na(flag))
nrow(locality_pts) #326017
  # save the file, for reference
write.csv(locality_pts, file.path(main_dir,data,
  paste0("need_geolocation_", Sys.Date(), ".csv")),
  row.names = F)

# create final subset of geolocated points which are on land
geo_pts <- all_data %>% filter(is.na(flag))
nrow(geo_pts) #2184049
table(geo_pts$database)
# BIEN    Ex_situ   FIA      GBIF     iDigBio   IUCN_RedList  NorthAm_herbaria
# 808255  1429      917467   344725   38286     22429         51458

################################################################################
# Standardize country code column for checking against lat-long later
################################################################################

# country name to 3 letter ISO code
  # fix some issues first (can add anything that is not matched unambiguously)
geo_pts$country <- mgsub(geo_pts$country,
    c("Bolívia","Brasil","EE. UU.","ESTADOS UNIDOS DE AMERICA",
      "México","MÉXICO","Repubblica Italiana","U. S. A.","United Statese",
      "America","Atats-Unis","CAN","ESP","MA(C)xico","MEX",
      "MX","PER","Unknown","Cultivated"),
    c("Bolivia","Brazil","United States","United States",
      "Mexico","Mexico","Italy","United States","United States",
      "United States","United States","Canada","Spain","Mexico","Mexico",
      "Mexico","Peru",NA,NA))
  # find codes for names
country_codes1 <- as.data.frame(sort(unique(geo_pts$country))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$country)),
      origin="country.name", destination="iso2c")) %>%
  rename("country" = "sort(unique(geo_pts$country))",
         "countryCode_standard" = "iso3c")
  # add to data
geo_pts <- left_join(geo_pts,country_codes1)

# country code to standard 2 letter ISO code
  # fix some issues (can add anything that is not matched unambiguously)
geo_pts$countryCode <- mgsub(geo_pts$countryCode,
    c("XK","ZZ"),c("RS",NA))
geo_pts$countryCode <- str_to_upper(geo_pts$countryCode)
  # check codes 
country_codes2 <- as.data.frame(sort(unique(geo_pts$countryCode))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$countryCode)),
      origin="iso2c", destination="iso2c")) %>%
  rename("countryCode" = "sort(unique(geo_pts$countryCode))",
         "countryCode_standard2" = "iso3c")
# add to data
geo_pts <- left_join(geo_pts,country_codes2)
geo_pts <- geo_pts %>% 
  unite("countryCode_standard",
        c("countryCode_standard","countryCode_standard2"),
        sep=";",remove=T,na.rm=T) %>%
  separate("countryCode_standard","countryCode_standard",sep=";",extra="drop")
geo_pts$countryCode_standard[which(geo_pts$countryCode_standard == "")] <- NA
sort(table(geo_pts$countryCode_standard))

################################################################################
# Remove spatial duplicates
################################################################################

## this section could potentially be moved to another script

## OTHER WAYS OF REMOVING DUPLICATES ARE ALSO POSSIBLE AND COULD MAKE MORE
##    SENSE FOR A SPECIFIC WAY OF USING THE POINTS, including
##    by grid cell, distance between points (e.g., randomly via spThin), etc...
## The segment below removes spatial duplicates based on rounded latitude
##    and longitude. This is a simple fix that doesn't involve spatial 
##    calculations. One additional positive of this method is that you can
##    choose priority datasets to keep points from; this can also help with 
##    citations since data from some databases are harder to cite than others.

# create rounded latitude and longitude columns for removing duplicates;
#   number of digits can be changed based on how dense you want data
geo_pts$lat_round <- round(geo_pts$decimalLatitude,digits=2)
geo_pts$long_round <- round(geo_pts$decimalLongitude,digits=2)

# create subset of all ex situ points, to add back in at end, if desired
ex_situ <- geo_pts[which(geo_pts$database=="Ex_situ"),]
ex_situ <- ex_situ %>% 
  dplyr::select(-flag) %>%
  mutate(coordinateUncertaintyInMeters=as.numeric(coordinateUncertaintyInMeters))

# sort before removing duplicates;
# whatever you sort to the top will be kept when there is a duplicate further down;
# you can turn any of these steps on/off, or add others
  # sort by basis of record
unique(geo_pts$basisOfRecord)
geo_pts$basisOfRecord <- factor(geo_pts$basisOfRecord,
  levels = c("PRESERVED_SPECIMEN","MATERIAL_SAMPLE","MATERIAL_CITATION",
             "OBSERVATION","HUMAN_OBSERVATION","OCCURRENCE","PHYSICAL_SPECIMEN",
             "MACHINE_OBSERVATION","FOSSIL_SPECIMEN","LIVING_SPECIMEN",
             "UNKNOWN"))
geo_pts <- geo_pts %>% arrange(basisOfRecord)
  # sort by establishment means
unique(geo_pts$establishmentMeans)
geo_pts$establishmentMeans <- factor(geo_pts$establishmentMeans,
  levels = c("NATIVE","WILD","UNCERTAIN","INTRODUCED","MANAGED","CULTIVATED",
             "DEAD")) #"INVASIVE"
geo_pts <- geo_pts %>% arrange(establishmentMeans)
  # sort by coordinate uncertainty
geo_pts$coordinateUncertaintyInMeters <-
  as.numeric(geo_pts$coordinateUncertaintyInMeters)
geo_pts <- geo_pts %>% arrange(geo_pts$coordinateUncertaintyInMeters)
  # sort by year
geo_pts <- geo_pts %>% arrange(desc(year))
  # sort by source database
unique(geo_pts$database)
geo_pts$database <- factor(geo_pts$database,
  levels = c("GBIF","NorthAm_herbaria","iDigBio",
             "IUCN_RedList","FIA","BIEN","Ex_situ"))
geo_pts <- geo_pts %>% arrange(database)

# remove duplicates
# can create "all_source_databases" column, to capture
#    databases from which duplicates were removed
# can take a while to remove duplicates if there are lots a rows
geo_pts2 <- geo_pts %>%
  group_by(taxon_name_accepted,lat_round,long_round) %>%
  mutate(all_source_databases = paste(unique(database), collapse = ', ')) %>%
  distinct(taxon_name_accepted,lat_round,long_round,.keep_all=T) %>%
  ungroup() %>%
  dplyr::select(-flag)
nrow(geo_pts2) #366757

# add ex situ data back in
geo_pts2 <- full_join(geo_pts2,ex_situ)
nrow(geo_pts2) #367559

## set header/column name order
keep_col <- c( #data source and unique ID
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
  "natureserve_rank","fruit_nut")
# set column order and remove a few unnecessary columns
geo_pts2 <- geo_pts2[,keep_col]

# take a look
head(as.data.frame(geo_pts2))
nrow(geo_pts2) #367559
table(geo_pts2$database)
# BIEN    Ex_situ   FIA      GBIF     iDigBio   IUCN_RedList  NorthAm_herbaria
# 11471   1429      160120   173202   2392      4903          14042
rm(geo_pts)

# write file if you'd like
#write.csv(geo_pts2, file.path(main_dir,data,standard,
#  paste0("Occurrences_Compiled_", Sys.Date(), ".csv")),row.names = F)

################################################################################
# 6. Look at results
################################################################################

# summarize results for each target taxon
  # lat-long records
count_geo <- geo_pts2 %>% count(taxon_name_accepted)
names(count_geo)[2] <- "num_latlong_records"
  # locality-only records (invalid coords or in water)
count_locality <- locality_pts %>% count(taxon_name_accepted)
names(count_locality)[2] <- "num_locality_records"
  # make table
files <- list(count_geo,count_locality)
summary <- setorder(Reduce(full_join, files),num_latlong_records,na.last=F)
head(summary)
  # write file
write.csv(summary, file.path(main_dir,data,
  paste0("occurrence_record_summary_", Sys.Date(), ".csv")),row.names = F)
#as.data.frame(summary)

################################################################################
# Split by species to save
################################################################################

# first, if you're working at the taxon level, add infrataxon records to their 
# parent species too
add_again <- geo_pts2 %>% filter(grepl("var\\.|subsp\\.", taxon_name_accepted))
unique(add_again$taxon_name_accepted)
add_again$taxon_name_accepted <- gsub(" var\\.*\\s.+", "",
                                      add_again$taxon_name_accepted)
unique(add_again$taxon_name_accepted)
geo_pts2 <- rbind(geo_pts2,add_again)
table(geo_pts2$database)
# BIEN    Ex_situ   FIA      GBIF     iDigBio   IUCN_RedList  NorthAm_herbaria
# 12824   1530      160120   177601   2755      4947          17201

# split records to create one CSV for each target taxon
sp_split <- split(geo_pts2, as.factor(geo_pts2$taxon_name_accepted))
names(sp_split) <- gsub(" ","_",names(sp_split))
names(sp_split) <- gsub("\\.","",names(sp_split))

# write files
lapply(seq_along(sp_split), function(i) write.csv(sp_split[[i]],
  file.path(main_dir,data,standard,"taxon_raw_points",
  paste0(names(sp_split)[[i]],".csv")),row.names = F))
