################################################################################

## 4-compile_occurrence_data.R

### Author: Emily Beckman Bruns
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: XX December 2022
### R version 4.2.2

### DESCRIPTION:
# This script compiles in situ occurrence point data and ex situ (genebank
# and botanical garden) data. We remove any rows for species not in our
# target taxa list, standardize some key columns, and write two CSVs of
# records: one with all records and one only with records that are geolocated
# (include a latitude and longitude).

### INPUTS:
# target_taxa_with_synonyms.csv
#

### OUTPUTS:
#
# folder (taxon_raw_records) with CSV of raw data for each target
# taxon (e.g., Malus_angustifolia.csv)
#
# folder (taxon_geo_points) with CSV of geolocated data for each target
# taxon (e.g., Malus_angustifolia.csv)
#
# CSV of all occurrence points without lat-long but with locality description
#   (need_geolocation.csv)
#
# summary table (occurrence_record_summary.csv) with one row for each
# target taxon, listing total number of records, number of records with valid
# lat-long, and number of records with locality description only
#

################################################################################
# Load libraries
################################################################################

my.packages <- c('plyr', 'tidyverse', 'data.table', 'textclean',
  'CoordinateCleaner', 'maps', 'rnaturalearth', 'rnaturalearthdata',
  'countrycode', 'raster', 'terra'
)
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

# or use 0-set_working_directory.R script:
source("SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")

################################################################################
################################################################################
# 1. Standardize occurrence data columns from different sources
################################################################################


# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "target_taxa_with_synonyms.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
target_sp <- unique(taxon_list$genus_species)

# exploring multicore...
#https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
#library(doParallel)
#numCores <- detectCores()
#numCores #8
#registerDoParallel(numCores)


################################################################################
################################################################################
# 2. Compile occurrence data
################################################################################

# create folder for output data
if(!dir.exists(file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R")))
  dir.create(file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R"), recursive=T)

# read in raw datasets
file_list <- list.files(file.path(main_dir,"occurrence_points",
  "standardized_occurrence_data"), pattern = ".csv", full.names = T)
file_dfs <- lapply(file_list, read.csv, header = T, na.strings = c("","NA"),
  colClasses = "character")
length(file_dfs) #7

# stack all datasets using rbind.fill, which keeps non-matching columns
#   and fills with NA; 'Reduce' iterates through list and merges with previous.
# this may take a few minutes if you have lots of data
all_data_raw <- Reduce(rbind.fill, file_dfs)
  nrow(all_data_raw) #582673
  names(all_data_raw) #37
  table(all_data_raw$database)
#     Ex_situ    GBIF           iDigBio      IUCN_RedList  US_Herbaria
#      10753       349328        38256        22721       161615

#add unique identifier
  nms <- names(all_data_raw)
all_data_raw <- all_data_raw %>% mutate(UID=paste0('id', sprintf("%08d",
  1:nrow(all_data_raw)))) %>% select(c('UID', all_of(nms)))
#rm(nms, file_dfs, file_list)
# all_data$UID <- seq.int(nrow(all_data))

## write out file to review if needed
#write.csv(all_data_raw, file.path(main_dir, "outputs",
#  paste0("all_data_cleaning_", Sys.Date(), ".csv")), row.names=FALSE)

################################################################################
# 3. Filter by target taxa
################################################################################

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "target_taxa_with_synonyms.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
taxon_list <- taxon_list %>%
  # if needed, add columns that separate out taxon name
  separate("taxon_name",c("genus","species","infra_rank","infra_name"),
    sep=" ",remove=F,fill="right") %>%
  # select necessary columns
  select(taxon_name,genus,species,infra_rank,
    infra_name,taxon_name_status,taxon_name_acc,species_name_acc)

# full join to taxon list
all_data_raw <- all_data_raw %>% select(-genus)
all_data_raw <- left_join(all_data_raw,taxon_list)
# join again just by species name if no taxon match
need_match <- all_data_raw[which(is.na(all_data_raw$taxon_name_status)),]
  nrow(need_match) #23281
  # remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data_raw)-ncol(taxon_list)+1)]
  # rename column for matching
need_match <- need_match %>% rename(taxon_name_full = taxon_name)
need_match$taxon_name <- need_match$species_name
  # new join
need_match <- left_join(need_match,taxon_list)
  # bind together new matches and previously matched
matched <- all_data_raw[which(!is.na(all_data_raw$taxon_name_status)),]
matched$taxon_name_full <- matched$taxon_name
all_data <- rbind(matched,need_match)
  table(all_data$taxon_name_status) # Accepted: 572711   Synonym: 9931

# check names that got excluded.....
still_no_match <- all_data[which(is.na(all_data$taxon_name_status)),]
  nrow(still_no_match) #4947760
table(still_no_match$database)
sort(table(still_no_match$taxon_name))
## write out file to review if needed
#write.csv(still_no_match, file.path(main_dir, "outputs",
#  paste0("no_taxon_match_", Sys.Date(), ".csv")),
#  row.names=FALSE)

# keep only rows for target taxa
all_data <- all_data[which(!is.na(all_data$taxon_name_status)),]
  nrow(all_data) #582642

### ! target taxa with no occurrence data:
unique(taxon_list$taxon_name_acc)[
  !(unique(taxon_list$taxon_name_acc) %in% (unique(all_data$taxon_name_acc)))]
# Prunus x orthosepala

#save(all_data, all_data_raw, file="all_data_to_clean.RData")
#  rm(still_no_match, matched, need_match, all_data_raw)

################################################################################
# 4. Standardize some key columns
################################################################################

#load("all_data_to_clean.RData")
## this section could potentially be moved to separate script

# create localityDescription column
all_data <- all_data %>%
  unite("localityDescription",
    c(locality,municipality,higherGeography,county,stateProvince,country,
      countryCode,locationNotes,verbatimLocality), remove = F, sep = " | ") %>%
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
all_data <- all_data %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "Introduced" = "INTRODUCED",
    "Uncertain" = "UNKNOWN",
    "Native" = "NATIVE"))
all_data$establishmentMeans[which(is.na(all_data$establishmentMeans))] <-
  "UNKNOWN"

# check validity of lat and long
  # if coords are both 0, set to NA
zero <- which(all_data$decimalLatitude == 0 & all_data$decimalLongitude == 0)
all_data$decimalLatitude[zero] <- NA; all_data$decimalLongitude[zero] <- NA
  # flag non-numeric and not available coordinates and lat > 90, lat < -90,
  # lon > 180, and lon < -180
coord_test <- cc_val(all_data, lon = "decimalLongitude",lat = "decimalLatitude",
  value = "flagged", verbose = TRUE) #Flagged 121093 records.
  # try switching lat and long for invalid points and check validity again
all_data[!coord_test,c("decimalLatitude","decimalLongitude")] <-
  all_data[!coord_test,c("decimalLongitude","decimalLatitude")]
coord_test <- cc_val(all_data, lon = "decimalLongitude",lat = "decimalLatitude",
  value = "flagged", verbose = TRUE) #Flagged 121093 records.
  ## mark these as flagged
all_data$flag <- NA
all_data[!coord_test,]$flag <- paste0("Coordinates invalid")
  # make invalid lat-long NA
#all_data[!coord_test,c("decimalLatitude","decimalLongitude")] <- c(NA,NA)

## set header/column name order
h.nms <- c("UID", "taxon_name_acc", "taxon_name", "scientificName",
  "taxonIdentificationNotes", "database", "year", "basisOfRecord",
  "establishmentMeans","decimalLatitude", "decimalLongitude",
  "coordinateUncertaintyInMeters", "geolocationNotes", "localityDescription",
  "county", "stateProvince", "country", "countryCode","institutionCode",
  "datasetName", "publisher", "rightsHolder", "license", "nativeDatabaseID",
  "references", "informationWithheld", "issue", "taxon_name_status", "flag")
# set column order and remove a few unnecessary columns
all_data <- all_data %>% select(all_of(h.nms))

# separate out points with locality description only (no lat-long)
locality_pts <- all_data %>% filter(!is.na(localityDescription) &
    !is.na(flag)) %>%
  arrange(desc(year)) %>%
  distinct(taxon_name_acc,localityDescription,.keep_all=T)
nrow(locality_pts) #82147

# move forward with subset of points that do have lat and long
geo_pts <- all_data %>%
  filter(is.na(flag)) #%>%
  #select(-localityDescription)
  nrow(geo_pts) #461549

# check if points are in water, mark, and separate out
world_polygons <- ne_countries(type = 'countries', scale = 'large')
# add buffer; 0.01 dd = ~ 0.4 to 1 km depending on location
world_buff <- buffer(world_polygons, width=0.04, dissolve=F)
  ## another option is data(buffland)
# check if in water and mark, then separate out
geo_pts[is.na(map.where(world_buff, geo_pts$decimalLongitude,
  geo_pts$decimalLatitude)),]$flag <- paste("Coordinates in water",sep="; ")
water_pts <- geo_pts %>% filter(grepl("water",flag))
  nrow(water_pts) #1524 --> there are ~7,000 now... I think this is due to the layer
                  # keeping in to see what's happening
table(water_pts$database)
#     Ex_situ         GBIF      iDigBio IUCN_RedList  US_Herbaria
#          88         1072          150           68          146

# add water points to locality points if they have locality data
#locality_pts_add <- water_pts %>% filter(!is.na(localityDescription)) %>%
#  arrange(desc(year)) %>%
#  distinct(taxon_name_acc,localityDescription,.keep_all=T)
#nrow(locality_pts_add) #1046
#locality_pts <- rbind(locality_pts,locality_pts_add)
# write file of locality-only points
#table(locality_pts$database)
#    Ex_situ        GBIF     iDigBio US_Herbaria
#       4193         771          78       73489
#write.csv(locality_pts, file.path(main_dir,"outputs",
#  paste0("need_geolocation_", Sys.Date(), ".csv")),
#  row.names = F)

# create final subset of geolocated points which are on land
#geo_pts <- geo_pts %>% filter(!grepl("Coordinates in water",flag))
#nrow(geo_pts) #450866
#table(geo_pts$database)
#     Ex_situ         GBIF      iDigBio IUCN_RedList  US_Herbaria
#        1357       340888        37339        22246        49036
# can write a file just to look it over
#write.csv(geo_pts, file.path(main_dir,"outputs","are_geolocated.csv"),
#  row.names = F)

# standardize country code column for checking against lat-long later
  # country name to 3 letter ISO code
    # fix some issues first
geo_pts$country <- mgsub(geo_pts$country,
    c("Bolívia","Brasil","EE. UU.","ESTADOS UNIDOS DE AMERICA",
      "México","MÉXICO","Repubblica Italiana","U. S. A.","United Statese"),
    c("Bolivia","Brazil","United States","United States",
      "Mexico","Mexico","Italy","United States","United States"))
country_set <- as.data.frame(sort(unique(geo_pts$country))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$country)),
      origin="country.name", destination="iso3c"))
names(country_set) <- c("country","iso3c")
  country_set[which(is.na(country_set$iso3c)),]
  # country code to 3 letter ISO code
geo_pts$countryCode <- str_to_upper(geo_pts$countryCode)
#geo_pts$countryCode <- mgsub(geo_pts$countryCode,
#    c("XK","ZZ"),c("SRB",NA))
country_set2 <- as.data.frame(sort(unique(geo_pts$countryCode))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$countryCode)),
      origin="iso2c", destination="iso3c"))
names(country_set2) <- c("countryCode","iso3c_2")
  country_set2[which(is.na(country_set2$iso3c)),]
country_set3 <- country_set2[which(is.na(country_set2$iso3c)),]
country_set3$iso3c_2 <- country_set3$countryCode
names(country_set3) <- c("countryCode","iso3c_3")
  country_set3[which(is.na(country_set3$iso3c)),]
  # add country codes to data
geo_pts <- join(geo_pts,country_set)
geo_pts <- join(geo_pts,country_set2)
geo_pts <- join(geo_pts,country_set3)
# errors here are ok! just ignore
geo_pts[which(geo_pts$iso3c == geo_pts$iso3c_2),]$iso3c_2 <- NA
geo_pts[which(geo_pts$iso3c == geo_pts$iso3c_3),]$iso3c_3 <- NA
geo_pts <- tidyr::unite(geo_pts,"countryCode_standard",
  c("iso3c","iso3c_2","iso3c_3"),sep=";",remove=T,na.rm=T)
sort(unique(geo_pts$countryCode_standard))
geo_pts$countryCode_standard[which(geo_pts$countryCode_standard == "")] <- NA

################################################################################
# 5. Remove duplicates
################################################################################

## this section could potentially be moved to another script as well
## OTHER WAYS OF REMOVING DUPLICATES ARE ALSO POSSIBLE AND COULD MAKE MORE
##    SENSE FOR A SPECIFIC WAY OF USING THE POINTS, including
##    by grid cell, distance between points, etc...
## The segement below removes spatial duplicates based on rounded lattitude
##    and longitude. This is a simple fix that doesn't involved spatial data
##    or complex spatial calculations.

# create rounded latitude and longitude columns for removing duplicates
#   number of digits can be changed based on how dense you want data
geo_pts$lat_round <- round(geo_pts$decimalLatitude,digits=1)
geo_pts$long_round <- round(geo_pts$decimalLongitude,digits=1)

# create subset of all ex situ points, to add back in at end, if desired
ex_situ <- geo_pts[which(geo_pts$database=="Ex_situ"),]

# sort before removing duplicates; can turn any of these on/off, or add others
  # sort by basis of record
geo_pts$basisOfRecord <- factor(geo_pts$basisOfRecord,
  levels = c("PRESERVED_SPECIMEN","MATERIAL_SAMPLE","OBSERVATION",
  "HUMAN_OBSERVATION","OCCURRENCE","MACHINE_OBSERVATION","LITERATURE",
  "FOSSIL_SPECIMEN","LIVING_SPECIMEN","UNKNOWN"))
geo_pts <- geo_pts %>% arrange(basisOfRecord)
  # sort by establishment means
geo_pts$establishmentMeans <- factor(geo_pts$establishmentMeans,
  levels = c("NATIVE","UNKNOWN","INTRODUCED","MANAGED"))#,"CUT","INVASIVE","DEAD"))
geo_pts <- geo_pts %>% arrange(establishmentMeans)
  # sort by coordiante uncertainty
geo_pts$coordinateUncertaintyInMeters <-
  as.numeric(geo_pts$coordinateUncertaintyInMeters)
geo_pts <- geo_pts %>% arrange(geo_pts$coordinateUncertaintyInMeters)
  # sort by year
geo_pts <- geo_pts %>% arrange(desc(year))
  # sort by dataset
geo_pts$database <- factor(geo_pts$database,
  levels = c(#"FIA",
    "GBIF","US_Herbaria","iDigBio",#"BISON","BIEN",
    "IUCN_RedList","Ex_situ"))
geo_pts <- geo_pts %>% arrange(database)

# remove duplicates
# can create "all_source_databases" column, to capture
#    databases from which duplicates were removed
# can take a while to remove duplicates if there are lots a rows
geo_pts2 <- geo_pts %>%
  group_by(taxon_name_acc,lat_round,long_round) %>%
  mutate(all_source_databases = paste(unique(database), collapse = ',')) %>%
  distinct(taxon_name_acc,lat_round,long_round,.keep_all=T) %>%
  ungroup() %>%
  select(-flag)

# add ex situ data back in
geo_pts2 <- geo_pts2 %>%
  filter(!grepl("Ex_situ",all_source_databases))
  # error here is ok if you don't have ex situ data!
ex_situ$all_source_databases <- "Ex_situ"
ex_situ_add <- ex_situ %>% arrange(UID) %>%
  select(basisOfRecord,establishmentMeans)
dups <- unique(geo_pts2$UID)
ex_situ <- ex_situ %>%
  select(-flag) %>%
  filter(!(UID %in% dups))
geo_pts2 <- rbind(geo_pts2,ex_situ)
geo_pts2 <- geo_pts2 %>% arrange(UID)
geo_pts2$basisOfRecord <- as.character(geo_pts2$basisOfRecord)
geo_pts2$establishmentMeans <- as.character(geo_pts2$establishmentMeans)
geo_pts2[which(grepl("Ex_situ",geo_pts2$all_source_databases)),8:9] <-
  ex_situ_add

## set header/column name order
h.nms2 <- c("taxon_name_acc", "taxon_name", "scientificName",
  "taxonIdentificationNotes", "database", "all_source_databases", "year",
  "basisOfRecord", "establishmentMeans","decimalLatitude", "decimalLongitude",
  "coordinateUncertaintyInMeters", "geolocationNotes", "localityDescription",
  "county", "stateProvince", "countryCode_standard", "institutionCode",
  "datasetName", "publisher", "rightsHolder", "license", "nativeDatabaseID",
  "references", "informationWithheld", "issue", "taxon_name_status", "UID")
# set column order and remove a few unnecessary columns
geo_pts2 <- geo_pts2 %>% select(all_of(h.nms2))

# take a look
head(geo_pts2)
nrow(geo_pts2) #91005
table(geo_pts2$all_source_databases)
table(geo_pts2$database)
#        GBIF  US_Herbaria      iDigBio IUCN_RedList      Ex_situ
#       78168         8274          969         2237         1357

write.csv(geo_pts2, file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R",
  paste0("Occurrences_Compiled_", Sys.Date(), ".csv")),row.names = F)

################################################################################
# 6. Look at results
################################################################################

# summarize results for each target taxon
  # lat-long records
count_geo <- geo_pts2 %>% count(taxon_name_acc)
names(count_geo)[2] <- "num_latlong_records"
  # water records
#count_water <- water_pts %>% count(taxon_name_acc)
#names(count_water)[2] <- "num_water_records"
  # locality-only records
count_locality <- locality_pts %>% count(taxon_name_acc)
names(count_locality)[2] <- "num_locality_records"
  # make table of all categories
files <- list(count_geo,count_locality)#count_water,
summary <- setorder(Reduce(full_join, files),num_latlong_records,na.last=F)
head(summary)
  # write file
write.csv(summary, file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R",
  paste0("occurrence_record_summary_", Sys.Date(), ".csv")),row.names = F)
as.data.frame(summary)

## can save data out to a file so don't have to rerun
#save(all_data, taxon_list, s, geo_pts2, need_match,
#  file=file.path(main_dir,"outputs","EO_data.RData"))
#  rm(all_data_raw, file_dfs, geo_pts, locality_pts, matched, need_match,
#  source_standard); head(all_data)

################################################################################
# 7. Split by species
################################################################################

#load(file.path(main_dir,"outputs","EO_data.RData"))

# split lat-long points to create one CSV for each target taxon
sp_split <- split(geo_pts2, as.factor(geo_pts2$taxon_name_acc))
names(sp_split) <- gsub(" ","_",names(sp_split))

# write files
if(!dir.exists(file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R","taxon_raw_points")))
  dir.create(file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R","taxon_raw_points"),
  recursive=T)
lapply(seq_along(sp_split), function(i) write.csv(sp_split[[i]],
  file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R","taxon_raw_points",
  paste0(names(sp_split)[[i]], ".csv")),row.names = F))

#unlink("all_data_to_clean.RData")
#unlink(file.path(main_dir, "outputs", paste0("all_data_cleaning_",
#Sys.Date(), ".csv")))

######
###### NOW RUN 3-1_refine_occurrence_points.R !!
######
