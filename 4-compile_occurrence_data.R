################################################################################

## 3-0_compile_occurrence_points_CWRcopy.R

### Author: Emily Beckman Bruns
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant #MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 26 May 2022
### Last updated: 10 November 2022

### R version 4.2.2

### DESCRIPTION:
  #
  # This script compiles in situ occurrence point data and ex situ (genebank
  # and botanical garden) data. We remove any rows for species not in our
  # target taxa list, standardize some key columns, and write two CSVs of
  # records: one with all records and one only with records that are geolocated
  # (include a latitude and longitude).
  #

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

# set manually to wherever you want
#   mine is a Google Drive folder
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

################################################################################
################################################################################
# 1. Standardize occurrence data columns from different sources
################################################################################

# create folders for output data
if(!dir.exists(file.path(main_dir,"occurrence_data")))
  dir.create(file.path(main_dir,"occurrence_data"), recursive=T)
if(!dir.exists(file.path(main_dir,"occurrence_data","standardized_occurrence_data")))
  dir.create(file.path(main_dir,"occurrence_data","standardized_occurrence_data"), recursive=T)

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

###############
### A) Global Biodiversity Information Facility (GBIF)
###############

# FIRST, download raw data at https://www.gbif.org/occurrence/search
#   Query: Scientific name = added each genus individually;
#          Location = Including coordinates
#   Citation for all data together:
#     GBIF.org (17 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.ngyqg7
#   Citations for data when broken into 3 (alphabetically 5 & 5, then synonym genera):
#     GBIF.org (17 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.yhe3xt
#     GBIF.org (17 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.rqqhg3
#     GBIF.org (17 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.p5r8rr
# Rename the downloaded folder "GBIF" and place in the "raw_occurrence_data" folder

# read in data (this is how you do it all at once)
#gbif_raw <- fread(file.path(main_dir,"occurrence_points","raw_occurrence_data",
#  "GBIF","occurrence.txt"), quote="", na.strings="")
#nrow(gbif_raw) #5214388

# read in raw data and loop through each file:
file_list <- list.files(path = file.path(main_dir,"occurrence_points",
  "raw_occurrence_data","GBIF"), pattern = "occurrence", full.names = T)
length(file_list) #3

for(i in 1:length(file_list)){

  # read in data
  gbif_raw <- fread(file_list[[i]], quote="", na.strings="")
  print(nrow(gbif_raw))

  ### standardize column names

  # remove genus-level records
  gbif_raw <- gbif_raw %>% filter(taxonRank != "GENUS")
  print(nrow(gbif_raw))

  # keep only necessary columns
  gbif_raw <- gbif_raw %>% select(
      # taxon name
    "scientificName",
    "family","genus","specificEpithet","taxonRank","infraspecificEpithet",
      # taxon IDs
    #"taxonID","speciesKey","taxonKey",
      # taxon identification notes (GROUP)
    "identificationRemarks","identificationVerificationStatus","identifiedBy",
      "taxonRemarks",
      # lat-long
    "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
      # record details
    "basisOfRecord","year","gbifID","references",
      #"identifier","occurrenceID","recordNumber",
      # locality description
    "locality","verbatimLocality","county","municipality","stateProvince",
      "higherGeography","countryCode",
      # location notes (GROUP)
    "associatedTaxa","eventRemarks","fieldNotes","habitat","locationRemarks",
      "occurrenceRemarks","occurrenceStatus",
      # geolocation details (GROUP)
    "georeferencedBy","georeferencedDate",
      "georeferenceProtocol","georeferenceRemarks","georeferenceSources",
      "georeferenceVerificationStatus",
      # data source details
    "datasetName","publisher","recordedBy","institutionCode",
      "rightsHolder","license",#"collectionCode"
      # other caveats
    "establishmentMeans","informationWithheld","issue"
    #"dataGeneralizations","hasGeospatialIssues"
  )
  gbif_raw$database <- "GBIF"

  # rename columns
  gbif_raw <- gbif_raw %>% rename(nativeDatabaseID = gbifID)

  # combine a few similar columns
  gbif_raw <- gbif_raw %>% unite("taxonIdentificationNotes",
      identificationRemarks:taxonRemarks,na.rm=T,remove=T,sep=" | ")
    gbif_raw$taxonIdentificationNotes <-
      gsub("^$",NA,gbif_raw$taxonIdentificationNotes)
  gbif_raw <- gbif_raw %>% unite("locationNotes",
    associatedTaxa:occurrenceStatus,na.rm=T,remove=T,sep=" | ")
    gbif_raw$locationNotes <- gsub("^$",NA,gbif_raw$locationNotes)
  gbif_raw <- gbif_raw %>% unite("geolocationNotes",
    georeferencedBy:georeferenceVerificationStatus,na.rm=T,remove=T,sep=" | ")
    gbif_raw$geolocationNotes <- gsub("^$",NA,gbif_raw$geolocationNotes)
  #str(gbif_raw)

  # create taxon_name column
  gbif_raw$taxon_name <- NA
  unique(gbif_raw$taxonRank)
  subsp <- gbif_raw %>% filter(taxonRank == "SUBSPECIES")
    subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
      subsp$infraspecificEpithet)
  var <- gbif_raw %>% filter(taxonRank == "VARIETY")
    var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
      var$infraspecificEpithet)
  form <- gbif_raw %>% filter(taxonRank == "FORM")
    form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
      form$infraspecificEpithet)
  spp <- gbif_raw %>% filter(taxonRank == "SPECIES")
    spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
  gbif_raw <- Reduce(rbind,list(subsp,var,form,spp))
    rm(subsp,var,form,spp)
  #str(gbif_raw)

  # fix hybrid names
  #unique(gbif_raw$taxon_name)
  gbif_raw$taxon_name <- mgsub(gbif_raw$taxon_name,
    c("Asimina xnashii","Prunus xorthosepala","Carya xlecontei",
      "Carya xludoviciana","Castanea xneglecta"),
    c("Asimina x nashii","Prunus x orthosepala","Carya x lecontei",
      "Carya x ludoviciana","Castanea x neglecta"))

  # create species_name column
  gbif_raw$species_name <- NA
  gbif_raw$species_name <- sapply(gbif_raw$taxon_name, function(x)
    unlist(strsplit(x," var. | subsp. | f. "))[1])
  #sort(unique(gbif_raw$species_name))

  # keep only target species
  gbif_raw <- gbif_raw %>%
    filter(species_name %in% target_sp)
  print(nrow(gbif_raw))

  # write file
  write.csv(gbif_raw, file.path(main_dir,"occurrence_points","standardized_occurrence_data",
    paste0("gbif",i,".csv")),row.names=FALSE)
  rm(gbif_raw)

}

###############
# B) Integrated Digitized Biocollections (iDigBio)
###############

# FIRST, download raw data at https://www.idigbio.org/portal/search
# Check the "Must have map point" checkbox
# Click "Add a field" dropdown on the left and select "Genus"; type your
#   target genus name into the "Genus" box
# Click the "Download" tab, type in your email, and click the download button
#   (down arrow within circle)
# If you have more than one target genus, repeat the above steps for the
#   other genera
# Your downloads will pop up in the "Downloads" section;
#   "Click To Download" for each
# Move all the folders you downloaded into a "iDigBio" folder
# Pull the "occurrence_raw.csv" file out into the
#   "iDigBio" folder and add appropriate genus name to the file name
# Place iDigBio folder in the "raw_occurrence_data" folder

# read in raw occurrence points
file_list <- list.files(path = file.path(main_dir,"occurrence_points",
  "raw_occurrence_data","iDigBio"), pattern = "occurrence_raw", full.names = T)
file_dfs <- lapply(file_list, read.csv, colClasses = "character",
  na.strings=c("","NA"),strip.white=T,fileEncoding="UTF-8")
length(file_dfs) #21

# stack datasets to create one dataframe
idigbio_raw <- data.frame()
for(file in seq_along(file_dfs)){
  idigbio_raw <- rbind(idigbio_raw, file_dfs[[file]])
};
  nrow(idigbio_raw) #148158
  rm(file_list,file_dfs)
# replace prefixes in column names
colnames(idigbio_raw) <- gsub(x = colnames(idigbio_raw),
  pattern = "dwc\\.", replacement = "")

### standardize column names

# combine a few similar columns
idigbio_raw <- idigbio_raw %>% unite("taxonIdentificationNotes",
    c(identificationID:identifiedBy,taxonRemarks),na.rm=T,remove=T,sep=" | ")
  idigbio_raw$taxonIdentificationNotes <-
    gsub("^$",NA,idigbio_raw$taxonIdentificationNotes)
idigbio_raw <- idigbio_raw %>% unite("locationNotes",
  c(associatedTaxa,eventRemarks,fieldNotes,habitat,locationRemarks,
    occurrenceRemarks,occurrenceStatus),na.rm=T,remove=T,sep=" | ")
  idigbio_raw$locationNotes <- gsub("^$",NA,idigbio_raw$locationNotes)
idigbio_raw <- idigbio_raw %>% unite("geolocationNotes",
  georeferenceProtocol:georeferencedDate,na.rm=T,remove=T,sep=" | ")
  idigbio_raw$geolocationNotes <- gsub("^$",NA,idigbio_raw$geolocationNotes)

# split date collected column to just get year
idigbio_raw <- idigbio_raw %>% separate("eventDate","year",sep="-",remove=T)
idigbio_raw$year <- gsub("[[:digit:]]+/[[:digit:]]+/","",idigbio_raw$year)
idigbio_raw$year <- gsub("^[1-9][0-9][0-9][0-9]/","",idigbio_raw$year)
idigbio_raw$year <- gsub("^[0-9][0-9]/","",idigbio_raw$year)
idigbio_raw$year <- gsub("/","",idigbio_raw$year)
idigbio_raw$year <- as.numeric(idigbio_raw$year)
# keep only years in reasonable timeframe (1500-current year)
idigbio_raw$year[which(idigbio_raw$year < 1500 |
  idigbio_raw$year > as.numeric(format(Sys.time(),"%Y")))] <- NA
  sort(unique(idigbio_raw$year))
nrow(idigbio_raw) #148158

# keep only necessary columns
idigbio_raw <- idigbio_raw %>% select(
    "scientificName","taxonIdentificationNotes",
    "family","genus","specificEpithet","taxonRank","infraspecificEpithet",
    "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
    "basisOfRecord","establishmentMeans","year",
    "locality","verbatimLocality","locationNotes","county","municipality",
    "higherGeography","stateProvince","country","countryCode",
    "institutionCode","datasetName","rightsHolder","dcterms.license","geolocationNotes",
    "recordedBy","coreid","dcterms.references","informationWithheld") %>%
  rename(license = dcterms.license, references = dcterms.references)
idigbio_raw$database <- "iDigBio"

# recode standard column: establishment means
idigbio_raw$establishmentMeans <- str_to_upper(idigbio_raw$establishmentMeans)
  unique(idigbio_raw$establishmentMeans)
idigbio_raw <- idigbio_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "CULTIVATED" = "MANAGED",
    "LAWN DECORATION" = "MANAGED",
    "CAPTIVE" = "MANAGED",
    "PLANTED" = "MANAGED",
    "MANAGED OR ESCAPED" = "MANAGED",
    "ESCAPED" = "INTRODUCED",
    "WILD." = "NATIVE",
    "CULTIVATED." = "MANAGED",
    "ALIEN" = "INVASIVE",
    "SHADE TREE PLANTED" = "MANAGED",
    .default = "UNKNOWN",
    .missing = "UNKNOWN"))

# create taxon_name column
sort(unique(idigbio_raw$taxonRank))
genus <- c("gen.","genus","Genus")
subspecies <- c("subsp.","subspecies","Subspecies")
variety <- c("var.","Varietas","variety","Variety")
forma <- c("f.","fm.","form","Form","form.","forma","Forma","Subform")
species <- c("espécie","sp.","specie","species","Species","spp.")
hybrid <- c("×","x","X","?")
#aff <- c("aff.")
idigbio_raw <- idigbio_raw %>% filter(!(taxonRank %in% genus))
subsp <- idigbio_raw %>% filter(taxonRank %in% subspecies)
  subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
    subsp$infraspecificEpithet)
var <- idigbio_raw %>% filter(taxonRank %in% variety)
  var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
    var$infraspecificEpithet)
form <- idigbio_raw %>% filter(taxonRank %in% forma)
  form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
    form$infraspecificEpithet)
spp <- idigbio_raw %>% filter(taxonRank %in% species)
  spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
h <- idigbio_raw %>% filter(taxonRank %in% hybrid)
  h$taxon_name <- paste0(h$genus," x ",h$specificEpithet)
#a <- idigbio_raw %>% filter(taxonRank %in% aff)
#  a$taxon_name <- paste(a$genus,"aff.",a$specificEpithet)
idigbio_raw <- Reduce(rbind.fill,list(subsp,var,form,spp,h))#a
  rm(subsp,var,form,spp,h)#a
  str(idigbio_raw)
# replace other hybrid signifiers
idigbio_raw$taxon_name <- gsub(" ×"," x ",idigbio_raw$taxon_name)
idigbio_raw$taxon_name <- gsub(" X "," x ",idigbio_raw$taxon_name)
idigbio_raw$taxon_name <- gsub(" \\?"," x ",idigbio_raw$taxon_name)
idigbio_raw$taxon_name <- str_squish(idigbio_raw$taxon_name)
sort(unique(idigbio_raw$taxon_name))

# create species_name column
idigbio_raw$species_name <- NA
idigbio_raw$species_name <- sapply(idigbio_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(idigbio_raw$species_name))

# recode standard columns
  # basis of record
unique(idigbio_raw$basisOfRecord)
idigbio_raw <- idigbio_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "PreservedSpecimen" = "PRESERVED_SPECIMEN",
    "Preserved Specimen" = "PRESERVED_SPECIMEN",
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "Preserved specimen" = "PRESERVED_SPECIMEN",
    "Preservedspecimen" = "PRESERVED_SPECIMEN",
    "Physical specimen" = "PRESERVED_SPECIMEN",
    "physicalspecimen" = "PRESERVED_SPECIMEN",
    "Physicalspecimen" = "PRESERVED_SPECIMEN",
    "PhysicalSpecimen" = "PRESERVED_SPECIMEN",
    "MaterialSample" = "MATERIAL_SAMPLE",
    "Occurrence" = "OCCURRENCE",
    "machineobservation" = "MACHINE_OBSERVATION",
    "HumanObservation" = "HUMAN_OBSERVATION",
    "FossilSpecimen" = "FOSSIL_SPECIMEN",
    .missing = "UNKNOWN"))
unique(idigbio_raw$basisOfRecord)

# keep only target species
idigbio_raw <- idigbio_raw %>%
  filter(species_name %in% target_sp)
nrow(idigbio_raw)

# write file
write.csv(idigbio_raw, file.path(main_dir, "occurrence_points",
  "standardized_occurrence_data", "idigbio.csv"), row.names=FALSE)
rm(idigbio_raw)

###############
# C) IUCN Red List
###############

# FIRST, download raw data:
  # If you don't have an IUCN RL account yet, create one, and sign in
  # Go to https://www.iucnredlist.org/search
  # Open the "Taxonomy" tab on the left, and type in your target genus name
  #   and check the box next to the genus name when it comes up;
  #   or, alternatively, if you are just looking for a few
  #   taxa you can search for them individually
  # You should be able to add each genus/taxon to your search so only
  #   one file needs to be exported
  # In the far-left bar, scroll down and, if desired, check
  #   "Subspecies and varieties"
  # Click the grey "Download" button and select "Range data - Points (CSV)";
  #   then fill in the popup window information
  # Go to https://www.iucnredlist.org/account to find your query
  # Click "Download" next to your query
  # Rename the folder you downloaded to "IUCN Red List" and place in
  #   the "raw_occurrence_data" folder

# read in raw occurrence points
redlist_raw <- read.csv(file.path(main_dir,"occurrence_points","raw_occurrence_data",
  "IUCN_Red_List", "points_data.csv"), colClasses = "character", na.strings=c("", "NA"),
  strip.white=T, fileEncoding="UTF-8")
nrow(redlist_raw) #70955

### standardize column names

# create taxon_name column
subsp <- redlist_raw %>% filter(!is.na(subspecies))
  subsp$taxon_name <- paste(subsp$binomial,"subsp.",subsp$subspecies)
spp <- redlist_raw %>% filter(is.na(subspecies))
  spp$taxon_name <- spp$binomial
redlist_raw <- rbind(subsp,spp)
sort(unique(redlist_raw$taxon_name))

# create rightsHolder column (to use for citations)
redlist_raw$rightsHolder <- paste(redlist_raw$citation,redlist_raw$year)

# standardize other columns
redlist_raw <- redlist_raw %>%
# keep only necessary columns
  select(taxon_name,binomial,tax_comm,event_year,
    basisofrec,origin,latitude,longitude,dist_comm,source,
    compiler,rightsHolder,presence,subspecies) %>%
# recode standard columns
  mutate(origin = recode(origin,
    "1" = "NATIVE",
    "2" = "REINTRODUCED",
    "3" = "INTRODUCED",
    "4" = "INTRODUCED", #VAGRANT on the RL
    "5" = "UNKNOWN", #ORIGIN_UNCERTAIN on the RL
    "6" = "ASSISTED_COLONISATION")) %>%
  mutate(presence = recode(presence,
    "1" = "EXTANT",
    "2" = "EXTANT",
    "3" = "POSSIBLY_EXTANT",
    "4" = "POSSIBLY_EXTINCT",
    "5" = "EXTINCT",
    "6" = "PRESENCE_UNCERTAIN")) %>%
  mutate(basisofrec = recode(basisofrec,
    "HumanObservation" = "HUMAN_OBSERVATION",
    "PreservedSpecimen" = "PRESERVED_SPECIMEN",
    "Literature" = "LITERATURE",
    "Expert" = "HUMAN_OBSERVATION",
    "FossilSpecimen" = "FOSSIL_SPECIMEN",
    "LivingSpecimen" = "LIVING_SPECIMEN",
    "Unknown" = "UNKNOWN",
    "Liturature" = "LITERATURE",
    "literature" = "LITERATURE",
    "Human_Observance" = "HUMAN_OBSERVATION",
    "Human_observation" = "HUMAN_OBSERVATION",
    "specimen" = "PRESERVED_SPECIMEN",
    "Online data" = "UNKNOWN",
    "Preserved_specimen" = "PRESERVED_SPECIMEN",
    .missing = "UNKNOWN")) %>%
# remove "EXTINCT" rows
  filter(presence != "EXTINCT") %>%
# rename columns
  rename(species_name = binomial,
                taxonIdentificationNotes = tax_comm,
                year = event_year,
                basisOfRecord = basisofrec,
                decimalLatitude = latitude,
                decimalLongitude = longitude,
                locality = dist_comm,
                datasetName = source,
                recordedBy = compiler,
                issue = presence,
                establishmentMeans = origin,
                taxonRank = subspecies)
redlist_raw$database <- "IUCN_RedList"
# check everything has been standardized
unique(redlist_raw$basisOfRecord)

# check data
head(redlist_raw)
nrow(redlist_raw) #70947

# keep only target species
redlist_raw <- redlist_raw %>%
  filter(species_name %in% target_sp)
nrow(redlist_raw)

# write file
write.csv(redlist_raw, file.path(main_dir,"occurrence_points","standardized_occurrence_data",
  "redlist.csv"),row.names=FALSE)
rm(redlist_raw)

###############
# D) U.S. Herbaria Consortia (SERNEC, SEINet, etc.)
###############

# FIRST, download raw data:
  # Go to http://sernecportal.org/portal/collections/harvestparams.php
  # Type your target genus/genera name(s) into the "scientific name" box
  #   and click "List Display"; or, alternatively, if you are just looking
  # for a few taxa you can search for and download them individually
  # Click the Download Specimen Data button (arrow pointing down into a box),
  #   in the top right corner
  # In the pop-up window, select the "Darwin Core" radio button,
  #   uncheck everything in the "Data Extensions" section, and
  #   select the "UTF-8 (unicode)" radio button
  #   leave other fields as-is
  # Click "Download Data"
  # Move the folder you downloaded into a "Herbaria Consortia" folder
  ## QUERY for CWR: Scientific Name = Asimina,Carya,Castanea,Corylus,Diospyros,Juglans,Malus,Persea,Pistacia,Prunus,Amygdalus,Annona,Celastrus,Cerasus,Emplectocladus,Fagus,Fagus-castanea,Hicoria,Hicorius,Laurocerasus,Laurus,Orchidocarpum,Padus,Pyrus,Tamala

# read in raw occurrence points
herb_raw <- read.csv(file.path(main_dir,"occurrence_points","raw_occurrence_data",
  "Herbaria_Consortia","occurrences.csv"), colClasses = "character",
  na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
nrow(herb_raw) #201998

### standardize column names

# create taxon_name column
  # this method is not perfect; the taxonRank isn't always categoried correctly
unique(herb_raw$taxonRank)
subsp <- herb_raw %>% filter(taxonRank == "Subspecies")
subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
  subsp$infraspecificEpithet)
var <- herb_raw %>% filter(taxonRank == "Variety")
var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
  var$infraspecificEpithet)
form <- herb_raw %>% filter(taxonRank == "Form")
form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
  form$infraspecificEpithet)
spp <- herb_raw %>% filter(
  is.na(taxonRank) | taxonRank == "Species" | taxonRank == "Subform")
  spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
herb_raw <- Reduce(rbind.fill,list(subsp,var,form,spp))
herb_raw$taxon_name[which(is.na(herb_raw$taxon_name))] <-
  herb_raw$scientificName[which(is.na(herb_raw$taxon_name))]
herb_raw$taxon_name <- gsub("Ã\u0097","",herb_raw$taxon_name)
herb_raw$taxon_name <- gsub("Ã«","e",herb_raw$taxon_name)
  #hybrid signifier replace
herb_raw$taxon_name <- gsub(" × "," x ",herb_raw$taxon_name)
herb_raw$taxon_name <- str_squish(herb_raw$taxon_name)
sort(unique(herb_raw$taxon_name))
  # some records got removed if they are genus-level
nrow(herb_raw) #201961

# keep only necessary columns
herb_raw <- herb_raw %>% select(
  "taxon_name",
  "family","genus","specificEpithet","taxonRank","infraspecificEpithet",
    "scientificName",
  #"taxonID",
  "identificationRemarks","identifiedBy","taxonRemarks",
  "decimalLatitude","decimalLongitude",
  "coordinateUncertaintyInMeters",
  "basisOfRecord","year","id","references",#"occurrenceID","recordNumber",
  "locality","county","municipality","stateProvince","country",
  "associatedTaxa","habitat","locationRemarks","occurrenceRemarks",
  "georeferencedBy","georeferenceProtocol","georeferenceRemarks",
    "georeferenceSources","georeferenceVerificationStatus",
  "institutionCode","rightsHolder","recordedBy",
  #"collectionCode",
  "establishmentMeans","informationWithheld")
herb_raw$database <- "US_Herbaria"
herb_raw$datasetName <- herb_raw$institutionCode

# rename columns
herb_raw <- herb_raw %>% rename(nativeDatabaseID = id)

# combine a few similar columns
herb_raw <- herb_raw %>% unite("taxonIdentificationNotes",
  identificationRemarks:taxonRemarks,na.rm=T,remove=T,sep=" | ")
  herb_raw$taxonIdentificationNotes <-
    gsub("^$",NA,herb_raw$taxonIdentificationNotes)
herb_raw <- herb_raw %>% unite("locationNotes",
  associatedTaxa:occurrenceRemarks,na.rm=T,remove=T,sep=" | ")
  herb_raw$locationNotes <- gsub("^$",NA,herb_raw$locationNotes)
herb_raw <- herb_raw %>% unite("geolocationNotes",
  georeferencedBy:georeferenceVerificationStatus,na.rm=T,remove=T,sep=" | ")
  herb_raw$geolocationNotes <- gsub("^$",NA,herb_raw$geolocationNotes)

# create species_name column
herb_raw$species_name <- NA
herb_raw$species_name <- sapply(herb_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(herb_raw$species_name))

# recode standard columns
  # basis of record
herb_raw <- herb_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "PreservedSpecimen" = "PRESERVED_SPECIMEN",
    "Preserved specimen" = "PRESERVED_SPECIMEN",
    "Preserved Specimen" = "PRESERVED_SPECIMEN",
    "preserved specimen" = "PRESERVED_SPECIMEN",
    "Espécimen preservado" = "PRESERVED_SPECIMEN",
    "Physicalspecimen" = "PRESERVED_SPECIMEN",
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "preservedSpecimen" = "PRESERVED_SPECIMEN",
    "Observation" = "OBSERVATION",
    "LivingSpecimen" = "LIVING_SPECIMEN",
    "Ejemplar herborizado" = "PRESERVED_SPECIMEN",
    "HumanObservation" = "HUMAN_OBSERVATION",
    .default = "UNKNOWN"))
  # year
sort(unique(herb_raw$year))
herb_raw <- herb_raw %>%
  mutate(year = recode(year,
    "9999" = "0",
    "18914" = "1891",
    "19418" = "1941"))
herb_raw$year <- as.integer(herb_raw$year)
herb_raw$year[which(herb_raw$year < 1500)] <- NA
sort(unique(herb_raw$year))
  # establishment means
sort(unique(herb_raw$establishmentMeans))
herb_raw <- herb_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "clonal" = "UNKNOWN",
    "Native" = "NATIVE",
    "Native." = "NATIVE",
    "native" = "NATIVE",
    "Wild." = "NATIVE",
    "wild" = "NATIVE",
    "wild collection" = "NATIVE",
    "Naturalized." = "INTRODUCED",
    "Alien" = "INTRODUCED",
    "Exotic" = "INTRODUCED",
    "Escape from cultivation." = "INTRODUCED",
    "Escaped from cultivation" = "INTRODUCED",
    "established non-native" = "INTRODUCED",
    "introduced" = "INTRODUCED",
    "Introduced" = "INTRODUCED",
    "Introduced; Volunteer" = "INTRODUCED",
    "Introduced." = "INTRODUCED",
    "Native/naturalizing" = "INTRODUCED",
    "Naturalized?" = "INTRODUCED",
    "Non-native" = "INTRODUCED",
    "Non-native invasive" = "INTRODUCED",
    "Non-native." = "INTRODUCED",
    "Nonnative" = "INTRODUCED",
    "Volunteer" = "INTRODUCED",
    "Uncertain" = "UNKNOWN",
    "wild caught" = "UNKNOWN",
    .default = "MANAGED"))

# check data
head(herb_raw)

# keep only target species
herb_raw <- herb_raw %>%
  filter(species_name %in% target_sp)
nrow(herb_raw)

# write file
write.csv(herb_raw, file.path(main_dir, "occurrence_points",
  "standardized_occurrence_data", "herbaria_consortia.csv"), row.names=FALSE)
rm(herb_raw)

###############
# E) Ex situ points
###############

# read in ex situ data, both dead (past) and current
exsitu_dead <- read.csv(list.files(path = file.path(main_dir,"occurrence_points",
  "raw_occurrence_data","Ex-situ"), pattern = "ExSitu_Dead", full.names = T),
  colClasses = "character", na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
nrow(exsitu_dead) #936
exsitu_postgeo <- read.csv(list.files(path = file.path(main_dir,"occurrence_points",
  "raw_occurrence_data","Ex-situ"), pattern = "ExSitu_Compiled_Post-Geolocation",
  full.names = T), colClasses = "character", na.strings=c("", "NA"), strip.white=T,
  fileEncoding="UTF-8")
nrow(exsitu_postgeo) #9817

# join the datasets together
exsitu_raw <- rbind.fill(exsitu_dead,exsitu_postgeo)
head(exsitu_raw)
nrow(exsitu_raw) #10753

# add centroid for county-level desciptions flagged during geolocating
  # read in shapefile of U.S. EPA Level IV ecoregions
    #https://catalog.data.gov/dataset/tiger-line-shapefile-current-nation-u-s-counties-and-equivalent-entities
#us_counties <- vect(file.path(main_dir,"gis_layers",
#  "tl_2021_us_county/tl_2021_us_county.shp"))
#county_cent <- centroids(us_counties,inside=TRUE)
#county_cent <- county_cent %>% mutate(
#  rename(NAME = county,
         ###THE STATE NAMES AREN"T HERE>>> THEY"RE JUST CODES :'''(((((( ))))))))
  # get records that need a centroid added and be sure county and state are filled in
#need_cent <- exsitu_raw %>%
#  filter(gps_det == "County-level") %>%
#  select(-lat_dd,long_dd)
  # remove the words 'county' and 'parish'
#need_cent$county <- gsub(" County","",need_cent$county)
#need_cent$county <- gsub(" Country","",need_cent$county)
#need_cent$county <- gsub(" Parish","",need_cent$county)
#  need_cent$county
  # add county and state when its not filled already:
#fix_cty <- need_cent %>%
#  filter(is.na(county) | is.na(state)) %>%
#  select(UID,all_locality,county,state)
#  fix_cty # see what is needed and edit table below:
#add <- data.frame(county = c("Clay","Scotland","Ontario","Oscoda","Oscoda",
#                             "Jackson","Jackson","Jackson","Jackson","Jackson",
#                             "Clark","Clark"),
#                  state = c("Arkansas","North Carolina","Canada","Michigan","Michigan",
#                            "Minnesota","Minnesota","Minnesota","Minnesota","Minnesota",
#                            "Arkansas","Arkansas"))
#fix_cty <- fix_cty %>% select(UID)
#fix_cty <- cbind(fix_cty,add)
#need_cent[need_cent$UID %in% fix_cty$UID,]$county <- fix_cty$county
#need_cent[need_cent$UID %in% fix_cty$UID,]$state <- fix_cty$state
  # add centroid data
#need_cent <-


### standardize column names

# keep only necessary columns
exsitu_raw <- exsitu_raw %>%
  select(
    "taxon_name_acc","species_name_acc","taxon_name","taxon_verif",
    "coll_year","prov_type","num_indiv","lat_dd","long_dd",
    "uncertainty","flag","all_locality","locality","germ_type",
    "garden_loc","state","inst_short","UID","gps_det") %>%
# rename columns
  rename(
    "taxon_name" = "taxon_name_acc",
    "species_name" = "species_name_acc",
    "scientificName" = "taxon_name",
    "taxonIdentificationNotes" = "taxon_verif",
    "year" = "coll_year",
    "basisOfRecord" = "prov_type",
    "establishmentMeans" = "num_indiv",
    "decimalLatitude" = "lat_dd",
    "decimalLongitude" = "long_dd",
    "coordinateUncertaintyInMeters" = "uncertainty",
    "issue" = "flag",
    "locality" = "all_locality",
    "verbatimLocality" = "locality",
    "stateProvince" = "state",
    "datasetName" = "inst_short",
    "nativeDatabaseID" = "UID",
    "geolocationNotes" = "gps_det") %>%
  # fix random error in inst_short
  mutate(datasetName = recode(datasetName,
    "UCalifornia BGerkeley" = "UCaliforniaBGBerkeley"))
exsitu_raw$database <- "Ex_situ"
exsitu_raw$rightsHolder <- exsitu_raw$datasetName

# combine a few similar columns
exsitu_raw <- exsitu_raw %>% unite("locationNotes",
  germ_type:garden_loc,na.rm=T,remove=T,sep=" | ")
  exsitu_raw$taxonIdentificationNotes <-
    gsub("^$",NA,exsitu_raw$locationNotes)

# create species_name column
exsitu_raw$species_name <- NA
exsitu_raw$species_name <- sapply(exsitu_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(exsitu_raw$species_name))

# check data
head(exsitu_raw)

# write file
write.csv(exsitu_raw, file.path(main_dir,"occurrence_points",
  "standardized_occurrence_data","ex_situ.csv"),row.names=FALSE)
rm(exsitu_raw)


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
