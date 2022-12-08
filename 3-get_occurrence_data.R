################################################################################

### 3-get_occurrence_data.R

### Authors: Emily Beckman Bruns & Shannon M Still
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Creation date: 05 February 2020
### Last full check and update: 17 November 2022
### R version 4.2.2

### DESCRIPTION:
# This script provides instructions and code chunks for downloading and
# standardizing wild occurrence points from:
  # GLOBAL DATABASES (though all likely have U.S. bias?)
    # A) Global Biodiversity Information Facility (GBIF)
    # B) Integrated Digitized Biocollections (iDigBio)
    # C) IUCN Red List of Threatened Species
    # D) North American herbaria consortia (SEINet, etc.)
    # E) Botanical Information and Ecology Network (BIEN)
  # U.S. NATIONAL DATABASES
    # F) Forest Inventory and Analysis (FIA) Program of USDA Forest Service
    # NOTE: BISON data used to be downloaded as well, but as of 2021
    # it has now been subsumed under GBIF.us, which is part of GBIF
  # EX SITU ACCESSION-LEVEL WILD COLLECTION LOCATION
    # G) If you ran script 2-compile_exsitu_data.R, this section can be 
    # used to prep wild collection location lat-longs to be used used within the
    # in situ occurrence data as well.
## NOTE #1: Not all data from these sources are reliable. The aim of this
#  script is to get all easily-downloadable occurrence data, which
#  can then be sorted and vetted for the user's specific purposes.
## NOTE #2: You can add other occurrence point data (e.g., expert comment,
#  NatureServe, floras, USDA PLANTS, BONAP, IUCN Red List, private
#  sources, etc.) by standardizing column names and formatting to match the
#  schema in the "Occurrence data column schema" tab of 
#  Gap-analysis-workflow_metadata, then save as a CSV in 
#  your standardized occurrence data folder
## NOTE #3: Each database has it's own citation guidelines; please review
#  them prior to using the data. Information has been compiled in the
#  "Occurrence data citation guidance" tab of Gap-analysis-workflow_metadata

### INPUTS:
# (optional) target_taxa_with_synonyms.csv
# see example in Gap-analysis-workflow_metadata workbook
  # columns:
    # 1. "taxon_name" --> genus, species, infra rank, and infra name, all
    #    separated by one space each; hybrid symbol should be " x ", rather
    #    than "_" or "✕", and go between genus name and species epithet
    # 2. "taxon_name_accepted" --> the accepted name you're using
     # 3. "taxon_name_status" --> "Accepted" or "Synonym"
     # 4+. (optional) other metadata you want to keep with target taxa

### OUTPUTS [all are optional]:
# gbif.csv
# idigbio.csv
# redlist.csv
# herbaria.csv
# bien.csv
# fia.csv
# exsitu.csv
# [additional files if added manually via instructions in NOTE #2 above]

################################################################################
# Load libraries
################################################################################

my.packages <- c('plyr', 'tidyverse', 'data.table', 'textclean',
                 'rgbif', 'ridigbio', 'BIEN')
# install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually...
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"
    
# ...or use 0-set_working_directory.R script:
source("SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")
    
# set up file structure within your main working directory
data <- "occurrence_data"
raw <- "raw_occurrence_data"
standard <- "standardized_occurrence_data"

################################################################################
# Load or create target taxa list
################################################################################

## IF YOU HAVE CSV OF TARGET TAXA AND SYNONYMS:
# read in taxa list
taxon_list <- read.csv(file.path(main_dir,"target_taxa_with_synonyms.csv"),
  header = T, colClasses="character")
head(taxon_list)
nrow(taxon_list) #287
# list of target taxon names
taxon_names <- sort(taxon_list$taxon_name)

## IF YOU HAVE JUST ONE OR A FEW TARGET TAXA, CREATE A LIST BY HAND:
#taxon_names <- c("Quercus havardii")

################################################################################
# Download & basic standardization of occurrence data from each target database
################################################################################

## NOTE #1: you can pick and choose any/all sections below (A-F), depending on
#  which databases you'd like to use
## NOTE #2: some sections have two options for downloading data: manually via
#  the website, or automacially using the API; choose whichever works for you

# create folder for all occurrence data
if(!dir.exists(file.path(main_dir,data)))
  dir.create(file.path(main_dir,data), recursive=T)
# create folder for raw data we will download
if(!dir.exists(file.path(main_dir,data,raw)))
  dir.create(file.path(main_dir,data,raw), recursive=T)
# create folder for data after basic standardization
if(!dir.exists(file.path(main_dir,data,standard)))
  dir.create(file.path(main_dir,data,standard), recursive=T)



###
###############
###############################################
### A) Global Biodiversity Information Facility (GBIF)
###    https://www.gbif.org
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,raw,"GBIF")))
  dir.create(file.path(main_dir,data,raw,"GBIF"),
  recursive=T)

###
### OPTION 1: automatic download via API
### (can go down to option 2 -manual download- if this isn't working)
###

# load GBIF account user information
  # if you don't have account yet, go to https://www.gbif.org then click
  # "Login" in top right corner, then click "Register"
# either read in a text file with username, password, and email (one on each
#   line) or manually fill in each below:
login <- read_lines(log_loc)
  user  <- login[1] #"user"
  pwd   <- login[2] #"password"
  email <- login[3] #"email"
  rm(login)
# get GBIF taxon keys for all taxa in target list
keys <- sapply(taxon_names,function(x) name_backbone(name=x)$speciesKey,
  simplify = "array"); keys
# remove duplicate and NULL keys
keys_nodup <- keys[!duplicated(keys) & keys != "NULL"]
# create data frame of keys and matching taxon_name
gbif_codes <- map_df(keys_nodup,~as.data.frame(.x),.id="taxon_name")
names(gbif_codes)[2] <- "speciesKey"
# create vector of keys as input into gbif download
gbif_taxon_keys <- vector(mode="numeric")
for(i in 1:length(keys_nodup)){
  gbif_taxon_keys <- c(gbif_taxon_keys,keys_nodup[[i]][1])
}; sort(gbif_taxon_keys)
  rm(i)
# download GBIF data (Darwin Core Archive format)
gbif_download <- occ_download(
                pred_in("taxonKey", gbif_taxon_keys),
                  ## you can also add additional filters:
                #pred_in("basisOfRecord", c("PRESERVED_SPECIMEN",
                #    "HUMAN_OBSERVATION","FOSSIL_SPECIMEN","OBSERVATION",
                #    "UNKNOWN","MACHINE_OBSERVATION","MATERIAL_SAMPLE",
                #    "LITERATURE")),
                #pred("hasCoordinate", TRUE),
                #pred("hasGeospatialIssue", FALSE),
                  ## "SIMPLE_CSV" format is another option
                format = "DWCA",
                user=user,pwd=pwd,
                email=email)
  rm(user, pwd, email)
# load gbif data just downloaded
  # download and unzip before reading in
download_key <- gbif_download
  # must wait for download to complete before continuing;
  # it may take a while (up to 3 hours) if you have a large taxa list;
  # function below will pause script until the download is ready;
  # you can also login on the gbif website and go to your profile to see the
  # progress of your download
occ_download_wait(download_key, status_ping=10, quiet=TRUE)
  # get download when its ready then unzip and read in
occ_download_get(key=download_key[1],
  path=file.path(main_dir,data,raw,"GBIF"),
  overwrite=TRUE) #Download file size: 184.11 MB
unzip(zipfile=paste0(
  file.path(main_dir,data,raw,"GBIF",download_key[1]),".zip"),
  files="occurrence.txt",
  exdir=file.path(main_dir,data,raw,"GBIF"))
  # create file list for next step (standardizing)
file_list <- list.files(path = file.path(main_dir,data,
  raw,"GBIF"), pattern = "occurrence", full.names = T)
length(file_list) #1

# NOTE !!! NOW SCROLL DOWN TO "STANDARDIZE THE DATA" SECTION !!!

###
### OPTION 2: manual download from website
###

# First, download raw data:
# Go to https://www.gbif.org/occurrence/search
# Note that you need an account and to be logged in.
# You can query by genus or scientific name, and use other filters as desired.
# Remember that you need to cite these data using the DOI provided with the
# download!
# Place all downloaded data (single or multiple) in the folder you just created:
# occurrence_data/raw_occurrence_data/GBIF

# get path(s) to raw data; we will then loop through each file to standardize it
file_list <- list.files(path = file.path(main_dir,data,raw,"GBIF"),
  pattern = "occurrence", full.names = T)
length(file_list) #3

###
### STANDARDIZE THE DATA
### Run this section no matter which option your chose for getting the data !!
###

# loop through file(s)
for(i in 1:length(file_list)){
  # read in data
  gbif_raw <- fread(file_list[[i]], quote="", na.strings="")
  print(nrow(gbif_raw))
  # remove genus-level records
  gbif_raw <- gbif_raw %>% filter(taxonRank != "GENUS")
  print(nrow(gbif_raw))
  # keep only necessary columns
  gbif_raw <- gbif_raw %>% select(
      # Taxon
    "scientificName","family","genus","specificEpithet","taxonRank",
    "infraspecificEpithet",
        # *concatenate into taxonIdentificationNotes
    "identificationRemarks","identificationVerificationStatus","identifiedBy",
    "taxonRemarks",
      # Event
    "year",
      # Record-level
    "basisOfRecord","gbifID","datasetName","publisher","institutionCode",
    "rightsHolder","license","references","informationWithheld","issue",
      # Occurrence
    "recordedBy","establishmentMeans",
      # Location
    "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
    "county","municipality","stateProvince","higherGeography","countryCode",
        # *concatenate into golocationNotes
    "georeferencedBy","georeferencedDate","georeferenceProtocol",
    "georeferenceRemarks","georeferenceSources",
    "georeferenceVerificationStatus",
        # *concatenate into locationNotes
    "locality","verbatimLocality","associatedTaxa","eventRemarks","fieldNotes",
    "habitat","locationRemarks","occurrenceRemarks","occurrenceStatus"
  ) %>% rename(nativeDatabaseID = gbifID)
  # add database column
  gbif_raw$database <- "GBIF"
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
  # fix hybrid names as needed (I don't think any were found)
  sort(unique(gbif_raw$taxon_name))
  gbif_raw$taxon_name <- mgsub(gbif_raw$taxon_name,
    c("Asimina xnashii","Prunus xorthosepala","Carya xlecontei",
      "Carya xludoviciana","Castanea xneglecta"),
    c("Asimina x nashii","Prunus x orthosepala","Carya x lecontei",
      "Carya x ludoviciana","Castanea x neglecta"))
  # create species_name column
  gbif_raw$species_name <- NA
  gbif_raw$species_name <- sapply(gbif_raw$taxon_name, function(x)
    unlist(strsplit(x," var. | subsp. | f. "))[1])
  # keep only target species
  gbif_raw <- gbif_raw %>%
    filter(species_name %in% taxon_names)
  print(nrow(gbif_raw))
  # check a few standards and recode if needed
    # establishmentMeans
  sort(unique(gbif_raw$establishmentMeans)) # check and add below as needed
  gbif_raw <- gbif_raw %>%
    mutate(establishmentMeans = recode(establishmentMeans,
      "Uncertain" = "UNCERTAIN",
      "Native" = "NATIVE",
      "Introduced" = "INTRODUCED"))
    # basisOfRecord
  sort(unique(gbif_raw$basisOfRecord))
  # write file
  write.csv(gbif_raw, file.path(main_dir,data,standard,
    paste0("gbif",i,".csv")),row.names=FALSE)
  rm(gbif_raw)
}



###
###############
###############################################
# B) Integrated Digitized Biocollections (iDigBio)
#    https://www.idigbio.org/portal/search
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,raw,"iDigBio")))
  dir.create(file.path(main_dir,data,raw,"iDigBio"),
  recursive=T)

###
### OPTION 1: automatic download via API
### (can go down to option 2 -manual download- if this isn't working;
###  also, all fields are not downloaded via the API... you can use
###  manual download method if you want *everything* that's available)
###

# download iDigBio data for target taxa
  # we have to go taxon by taxon; function can only return 100,000 records at
  # once and e.g. Quercus has more than that so can't download by genus
idigbio_raw <- data.frame()
for(i in 1:length(taxon_names)){
  output_new <- idig_search_records(rq=list(scientificname=taxon_names[[i]]),
    fields="all")
  idigbio_raw <- rbind(idigbio_raw,output_new)
  print(paste(round(i/length(taxon_names)*100,digits=1),"% complete",sep=""))
}
  rm(output_new, i)
nrow(idigbio_raw) #128374
# remove rows that are lists
idigbio_raw <- idigbio_raw %>% dplyr::select(everything(),-commonnames,-flags,
  -mediarecords,-recordids)
# write raw file
write.csv(idigbio_raw,file.path(main_dir,data,
  raw,"iDigBio","idigbio_R_download.csv"),row.names=FALSE)
# standarize a few fields:
# split date collected column to just get year
idigbio_raw <- idigbio_raw %>% separate("eventdate","year",sep="-",remove=T)
idigbio_raw$year <- gsub("[[:digit:]]+/[[:digit:]]+/","",idigbio_raw$year)
  sort(unique(idigbio_raw$year))
# keep only necessary columns & rename to fit standard
idigbio_raw <- idigbio_raw %>%
  select("scientificname","family","genus","specificepithet","taxonrank",
         "infraspecificepithet","year","basisofrecord","uuid","institutioncode",
         "occurrenceid","collectioncode","individualcount","geopoint.lon",
         "geopoint.lat","coordinateuncertainty","locality","verbatimlocality",
         "municipality","county","stateprovince","country","countrycode") %>%
  rename("scientificName" = "scientificname",
         "specificEpithet" = "specificepithet",
         "taxonRank" = "taxonrank",
         "infraspecificEpithet" = "infraspecificepithet",
         "decimalLongitude" = "geopoint.lon",
         "decimalLatitude" = "geopoint.lat",
         "coordinateUncertaintyInMeters" = "coordinateuncertainty",
         "basisOfRecord" = "basisofrecord",
         "nativeDatabaseID" = "uuid",
         "references" = "occurrenceid",
         "verbatimLocality" = "verbatimlocality",
         "stateProvince" = "stateprovince",
         "countryCode" = "countrycode",
         "institutionCode" = "institutioncode",
         "individualCount" = "individualcount")
idigbio_raw$datasetName <- idigbio_raw$institutionCode
idigbio_raw$rightsHolder <- idigbio_raw$institutionCode
idigbio_raw$establishmentMeans <- "UNKNOWN"
# create taxon_name column
idigbio_raw$taxon_name <- str_to_sentence(idigbio_raw$scientificName)
sort(unique(idigbio_raw$taxon_name))
# fix name capitalization
idigbio_raw$family <- str_to_title(idigbio_raw$family)
idigbio_raw$genus <- str_to_title(idigbio_raw$genus)

# NOTE !!! NOW SCROLL DOWN TO "STANDARDIZE THE DATA" SECTION !!!

###
### OPTION 2: manual download from website
###

# First, download raw data:
# Go to https://www.idigbio.org/portal/search
# Click "Add a field" dropdown on the left and select "Genus"
# Type your target genus name into the "Genus" box
# Click the "Download" tab, type in your email, and
# click the download button (down arrow within circle).
# If you have more than one target genus, repeat the above steps for the
# other genera.
# Your downloads will pop up in the "Downloads" section
# "Click To Download" for each
# Move all the folders you downloaded into the "iDigBio" folder
# in your working directory.
# Pull the "occurrence_raw.csv" file out into the "iDigBio" folder.

# read in raw occurrence points
file_list <- list.files(path = file.path(main_dir,data,raw,"iDigBio"),
  pattern = "occurrence", full.names = T)
file_dfs <- lapply(file_list, read.csv, colClasses = "character",
  na.strings=c("","NA"),strip.white=T,fileEncoding="UTF-8")
length(file_dfs) #21
# stack datasets to create one dataframe
idigbio_raw <- data.frame()
for(file in seq_along(file_dfs)){
  idigbio_raw <- rbind(idigbio_raw, file_dfs[[file]])
}; nrow(idigbio_raw) #148158
# replace prefixes in column names
colnames(idigbio_raw) <- gsub(x = colnames(idigbio_raw),
  pattern = "dwc\\.", replacement = "")
# create taxon_name column
sort(unique(idigbio_raw$taxonRank))
subspecies <- c("ssp.","subsp.","subspecies","Subspecies")
variety <- c("var.","Varietas","variety","Variety")
forma <- c("f.","fm.","form","Form","form.","forma","Forma","Subform")
species <- c("espécie","sp.","specie","species","Species","spp.")
hybrid <- c("×","x","X")
aff <- c("aff.")
  # note that it's ok for these to error "replacement has 1 row, data has 0"
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
  h$taxon_name <- paste(h$genus,"x",h$specificEpithet)
a <- idigbio_raw %>% filter(taxonRank %in% aff)
  a$taxon_name <- paste(a$genus,"aff.",a$specificEpithet)
idigbio_raw <- Reduce(rbind.fill,list(subsp,var,form,spp,h,a))
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
  sort(unique(idigbio_raw$year))
# keep only necessary columns & rename to fit standard
idigbio_raw <- idigbio_raw %>%
  select(
    "taxon_name","scientificName","taxonIdentificationNotes",
    "family","genus","specificEpithet","taxonRank","infraspecificEpithet",
    "year","basisOfRecord","coreid","institutionCode","datasetName",
    "rightsHolder","dcterms.license","dcterms.references","informationWithheld",
    "recordedBy","establishmentMeans","individualCount","decimalLongitude",
    "decimalLatitude","coordinateUncertaintyInMeters",
    "geolocationNotes","locality","verbatimLocality","locationNotes",
    "county","municipality","higherGeography","stateProvince","country",
    "countryCode") %>%
  rename("license" = "dcterms.license",
         "references" = "dcterms.references",
         "nativeDatabaseID" = "coreid")

###
### STANDARDIZE THE DATA
### Run this section no matter which option your chose for getting the data !!
###

# add database column
idigbio_raw$database <- "iDigBio"
# create species_name column
idigbio_raw$species_name <- NA
idigbio_raw$species_name <- sapply(idigbio_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
# keep only target species
idigbio_raw <- idigbio_raw %>%
  filter(species_name %in% taxon_names)
# check a few standards and recode if needed
  # establishmentMeans
idigbio_raw$establishmentMeans <- str_to_upper(idigbio_raw$establishmentMeans)
unique(idigbio_raw$establishmentMeans) # check and add below as needed
idigbio_raw <- idigbio_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "WILD: NATIVE/NATURALIZED" = "NATIVE",
    "WILD." = "NATIVE",
    "WILD" = "NATIVE",
    "WILD COLLECTION" = "NATIVE",
    "WILD CAUGHT" = "NATIVE",
    "LAWN DECORATION" = "CULTIVATED",
    "CAPTIVE" = "CULTIVATED",
    "SHADE TREE PLANTED" = "CULTIVATED",
    "CULTIVATED." = "CULTIVATED",
    "MANAGED OR ESCAPED" = "MANAGED",
    "ESCAPED" = "INTRODUCED",
    "ALIEN" = "INTRODUCED",
    "PERSISTING FROM CULTIVATION" = "INTRODUCED",
    .default = "UNCERTAIN",
    .missing = "UNCERTAIN"))
  # basisOfRecord
idigbio_raw$basisOfRecord <- str_to_lower(idigbio_raw$basisOfRecord)
idigbio_raw$basisOfRecord <- gsub(" |_","",idigbio_raw$basisOfRecord)
unique(idigbio_raw$basisOfRecord) # check and add below as needed
idigbio_raw <- idigbio_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "machineobservation" = "MACHINE_OBSERVATION",
    "humanobservation" = "HUMAN_OBSERVATION",
    "fossilspecimen" = "FOSSIL_SPECIMEN",
    "materialsample" = "MATERIAL_SAMPLE",
    "physicalspecimen" = "PHYSICAL_SPECIMEN",
    "occurrence" = "OCCURRENCE",
    .default = "UNKNOWN",
    .missing = "UNKNOWN"))
# write file
write.csv(idigbio_raw, file.path(main_dir,data,standard,"idigbio.csv"),
  row.names=FALSE)
rm(idigbio_raw)



###
###############
###############################################
# C) IUCN Red List of Threatened Species
#    https://www.iucnredlist.org
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,raw,"IUCN_RedList")))
  dir.create(file.path(main_dir,data,raw,"IUCN_RedList"),
  recursive=T)

# There is an API for the IUCN Red List ('rredlist')
#   https://cran.r-project.org/web/packages/rredlist/rredlist.pdf
# but it does not have a function for downloading point data
# so, we must download manually.

# First, download raw data:
# If you don't have an IUCN RL account yet, create one
# Sign in to your account
# Go to https://www.iucnredlist.org/search (or click "Advanced" on home page)
# Open the "Taxonomy" tab on the left, type in your target genus name
# and check the box next to the genus name when it comes up - note that
# sometimes there is a glitch where the box is greyed out; you can navigate
# to the genus you want by going through the tree, starting with Plantae.
# Or, alternatively, if you are just looking for a few
# taxa you can search for them individually.
# You should be able to add each genus/taxon to your search so only
# one file needs to be exported.
# In the far-left bar, scroll down and, if desired, check
# "Subspecies and varieties" (if you want infrataxa).
# Click the grey "Download" button and select "Range data - Points (CSV)"
# then fill in the prompts in the popup window.
# Go to https://www.iucnredlist.org/account to find your query
# Click "Download" next to your query
# Move the folder you downloaded into the "IUCN_RedList" folder in the
# "occurrence_data/raw_occurrence_data" folder
# Pull the "points_data.csv" file out into the "IUCN_RedList" folder

# read in data
redlist_raw <- read.csv(file.path(main_dir,data,raw,
  "IUCN_RedList","points_data.csv"), colClasses = "character",
  na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
nrow(redlist_raw) #52330
# create taxon_name column
subsp <- redlist_raw %>% filter(!is.na(subspecies))
  subsp$taxon_name <- paste(subsp$binomial,"subsp.",subsp$subspecies)
spp <- redlist_raw %>% filter(is.na(subspecies))
  spp$taxon_name <- spp$binomial
redlist_raw <- rbind(subsp,spp)
sort(unique(redlist_raw$taxon_name))
# combine a few similar columns
redlist_raw$rightsHolder <- paste(redlist_raw$citation,redlist_raw$year)
# keep only necessary columns
redlist_raw <- redlist_raw %>%
  select(taxon_name,binomial,tax_comm,event_year,
    basisofrec,origin,latitude,longitude,dist_comm,source,
    compiler,rightsHolder,presence,subspecies,id_no)
# check a few standards and recode if needed
  # establishmentMeans
redlist_raw <- redlist_raw %>%
  mutate(origin = recode(origin,
    "1" = "NATIVE",
    "2" = "REINTRODUCED",
    "3" = "INTRODUCED",
    "4" = "VAGRANT",
    "5" = "UNCERTAIN",
    "6" = "ASSISTED_COLONISATION"))
  # issue
redlist_raw <- redlist_raw %>%
  mutate(presence = recode(presence,
    "1" = "EXTANT",
    "2" = "EXTANT",
    "3" = "POSSIBLY_EXTANT",
    "4" = "POSSIBLY_EXTINCT",
    "5" = "EXTINCT",
    "6" = "PRESENCE_UNCERTAIN")) %>%
    # remove extinct rows
  filter(presence != "EXTINCT")
  # basisOfRecord
redlist_raw$basisofrec <- str_to_lower(redlist_raw$basisofrec)
redlist_raw$basisofrec <- gsub(" |_","",redlist_raw$basisofrec)
unique(redlist_raw$basisofrec) # check and add below as needed
redlist_raw <- redlist_raw %>%
  mutate(basisofrec = recode(basisofrec,
    "humanobservation" = "HUMAN_OBSERVATION",
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "literature" = "LITERATURE",
    "expert" = "HUMAN_OBSERVATION",
    "fossilspecimen" = "FOSSIL_SPECIMEN",
    "livingspecimen" = "LIVING_SPECIMEN",
    "unknown" = "UNKNOWN",
    "liturature" = "LITERATURE",
    "materialsample" = "MATERIAL_SAMPLE",
    "observation" = "OBSERVATION",
    "specimen" = "PRESERVED_SPECIMEN",
    .missing = "UNKNOWN"))
# rename to fit standard
redlist_raw <- redlist_raw %>%
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
         taxonRank = subspecies,
         nativeDatabaseID = id_no,
         references = compiler)
# add database column & publisher columns
redlist_raw$database <- "IUCN_RedList"
redlist_raw$publisher <- "IUCN Red List of Threatened Species"
# write file
write.csv(redlist_raw, file.path(main_dir,data,standard,"redlist.csv"),
  row.names=FALSE)
rm(redlist_raw)



###
###############
###############################################
# D) North American herbaria consortia (SERNEC, SEINet, etc.)
#    https://sernecportal.org/portal/collections/index.php
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,standard,"NorthAm_herbaria")))
  dir.create(file.path(main_dir,data,raw,"NorthAm_herbaria"),
  recursive=T)

# First, download raw data:
# Go to http://sernecportal.org/portal/collections/harvestparams.php
# Type your target genus name into the "scientific name" box and click
# "List Display"; or, alternatively, if you are just looking for a few
# taxa you can search for and download them individually.
# Click the Download Specimen Data button (grey square with arrow pointing
# down into a box) in the top right corner.
# In the pop-up window, select the "Darwin Core" radio button,
# uncheck everything in the "Data Extensions" section, and
# select the "UTF-8 (unicode)" radio button; leave other fields as-is.
# Click "Download Data"
# If you have more than one target genus, repeat the above steps for the
# other genera.
# Move all the folders you downloaded into the "NorthAm_herbaria" folder.

# read in data
# this code pull the "occurrences.csv" from each genus folder, for compilation
file_list <- list.files(file.path(main_dir,data,raw,"NorthAm_herbaria"),
  pattern = "SymbOutput", full.names = T)
file_dfs <- lapply(file_list, function(i){
  read.csv(file.path(i,"occurrences.csv"), colClasses = "character",
   na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")})
sernec_raw <- data.frame()
for(file in seq_along(file_dfs)){
  sernec_raw <- rbind(sernec_raw, file_dfs[[file]])
}; nrow(sernec_raw) #182142
# remove genus-level records
sernec_raw <- sernec_raw %>% filter(taxonRank != "Genus")
# create taxon_name column
# this method is not perfect; the taxonRank isn't always categoried correctly
subsp <- sernec_raw %>% filter(taxonRank == "Subspecies")
subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
  subsp$infraspecificEpithet)
var <- sernec_raw %>% filter(taxonRank == "Variety")
var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
  var$infraspecificEpithet)
form <- sernec_raw %>% filter(taxonRank == "Form")
form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
  form$infraspecificEpithet)
spp <- sernec_raw %>% filter(is.na(taxonRank) | taxonRank == "Species" |
  taxonRank == "Subform")
  spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
sernec_raw <- Reduce(rbind.fill,list(subsp,var,form,spp))
sernec_raw$taxon_name[which(is.na(sernec_raw$taxon_name))] <-
  sernec_raw$scientificName[which(is.na(sernec_raw$taxon_name))]
sort(unique(sernec_raw$taxon_name))
sernec_raw$taxon_name <- gsub("Ã\u0097","",sernec_raw$taxon_name)
sernec_raw$taxon_name <- gsub("Ã«","e",sernec_raw$taxon_name)
sernec_raw$taxon_name <- str_squish(sernec_raw$taxon_name)
# keep only necessary columns & rename to fit standard
sernec_raw <- sernec_raw %>% select(
  "taxon_name","family","genus","specificEpithet","taxonRank",
  "infraspecificEpithet","scientificName","identificationRemarks",
  "identifiedBy","taxonRemarks","decimalLatitude","decimalLongitude",
  "coordinateUncertaintyInMeters","basisOfRecord","year","id",
  "references","locality","county","municipality","stateProvince",
  "country","associatedTaxa","habitat","locationRemarks",
  "occurrenceRemarks","georeferencedBy","georeferenceProtocol",
  "georeferenceRemarks","georeferenceSources",
  "georeferenceVerificationStatus","institutionCode",
  "rightsHolder","recordedBy","individualCount",
  "establishmentMeans","informationWithheld")
sernec_raw <- sernec_raw %>% rename(nativeDatabaseID = id)
# add database column & datasetName column
sernec_raw$database <- "NorthAm_herbaria"
sernec_raw$datasetName <- sernec_raw$institutionCode
# combine a few similar columns
sernec_raw <- sernec_raw %>% unite("taxonIdentificationNotes",
  identificationRemarks:taxonRemarks,na.rm=T,remove=T,sep=" | ")
  sernec_raw$taxonIdentificationNotes <-
    gsub("^$",NA,sernec_raw$taxonIdentificationNotes)
sernec_raw <- sernec_raw %>% unite("locationNotes",
  associatedTaxa:occurrenceRemarks,na.rm=T,remove=T,sep=" | ")
  sernec_raw$locationNotes <- gsub("^$",NA,sernec_raw$locationNotes)
sernec_raw <- sernec_raw %>% unite("geolocationNotes",
  georeferencedBy:georeferenceVerificationStatus,na.rm=T,remove=T,sep=" | ")
  sernec_raw$geolocationNotes <- gsub("^$",NA,sernec_raw$geolocationNotes)
# create species_name column
sernec_raw$species_name <- NA
sernec_raw$species_name <- sapply(sernec_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(sernec_raw$species_name))
# check a few standards and recode if needed
  # basisOfRecord
sernec_raw$basisOfRecord <- str_to_upper(sernec_raw$basisOfRecord)
sernec_raw$basisOfRecord <- gsub(" ","_",sernec_raw$basisOfRecord)
unique(sernec_raw$basisOfRecord) # check and add below as needed
sernec_raw <- sernec_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "PRESERVEDSPECIMEN" = "PRESERVED_SPECIMEN",
    "PHYSICALSPECIMEN" = "PRESERVED_SPECIMEN",
    "ESPÉCIMEN_PRESERVADO" = "PRESERVED_SPECIMEN",
    "EJEMPLAR_HERBORIZADO" = "PRESERVED_SPECIMEN",
    "LIVINGSPECIMEN" = "LIVING_SPECIMEN",
    "HUMANOBSERVATION" = "HUMAN_OBSERVATION",
    .missing = "UNKNOWN"))
  # establishmentMeans
sernec_raw$establishmentMeans <- str_to_upper(sernec_raw$establishmentMeans)
sernec_raw$establishmentMeans <- gsub("\\.","",sernec_raw$establishmentMeans)
sort(unique(sernec_raw$establishmentMeans))
sernec_raw <- sernec_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "NATIVE" = "NATIVE",
    "INTRODUCED" = "INTRODUCED",
    "UNCERTAIN" = "UNCERTAIN",
    "ALIEN" = "INTRODUCED",
    "CLONAL" = "UNKNOWN",
    "WILD" = "NATIVE",
    "NATURALIZED" = "INTRODUCED",
    "ESCAPE FROM CULTIVATION" = "INTRODUCED",
    "ESTABLISHED NON-NATIVE" = "INTRODUCED",
    "INTRODUCED; VOLUNTEER" = "INTRODUCED",
    "NATIVE/NATURALIZING" = "NATIVE",
    "NATURALIZED" = "INTRODUCED",
    "NON-NATIVE" = "INTRODUCED",
    "NONNATIVE" = "INTRODUCED",
    "WILD CAUGHT" = "NATIVE",
    "VOLUNTEER" = "UNCERTAIN",
    "WILD COLLECTION" = "WILD",
    .default = "CULTIVATED"))
# write file
write.csv(sernec_raw, file.path(main_dir,data,standard,"herbaria.csv"),
  row.names=FALSE)
rm(sernec_raw)



###
###############
###############################################
# E) Botanical Information and Ecology Network (BIEN)
#    https://bien.nceas.ucsb.edu/bien/
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,raw,"BIEN")))
  dir.create(file.path(main_dir,data,raw,"BIEN"),
  recursive=T)

# This database can only be accessed via an API (no manual download option)
# For information about the functions in their package, run the next line:
#vignette("BIEN")

# Note that some people have received an error message that they are
# unable to connect to the correct port; for example:
  #   Error in postgresqlNewConnection(drv, ...) :
  #   RS-DBI driver: (could not connect public_bien@vegbiendev.nceas.ucsb.edu:5432 on dbname "public_vegbien": could not connect to server: Connection timed out (0x0000274C/10060)
  #                   Is the server running on host "vegbiendev.nceas.ucsb.edu" (128.111.85.31) and accepting
  #                   TCP/IP connections on port 5432?
# The vignette("BIEN") says this about the issue:
## Database connection issues Some institution and computer programs
## (e.g. some antivirus programs) block the SQL connections that this package
## relies on. While we are exploring ways around this issue, at present the
## simplest method is to use the package on a computer/network that doesn’t
## block SQL connections.

# download occurrence data for target taxa
# this can take a few minutes or more if you have many target taxa
bien_raw <- BIEN_occurrence_species(taxon_names,all.taxonomy=T,native.status=T,
  natives.only=F,observation.type=T,collection.info=T,political.boundaries=T,
  cultivated=T)
nrow(bien_raw) #1514733
write.csv(bien_raw, file.path(main_dir,data,
  raw,"BIEN","bien_R_download.csv"),row.names=FALSE)
# get citation info and write to a file
bien_citation <- BIEN_metadata_citation(dataframe = bien_raw)
write.csv(bien_citation, file.path(main_dir,data,
  raw,"BIEN","citation_info.csv"),row.names=FALSE)

# OPTIONAL: remove rows from GBIF and/or FIA if you are separately downloading
# those datasets
bien_raw <- bien_raw %>% filter(datasource != "FIA" & datasource != "GBIF")
nrow(bien_raw) #871544

# split date collected column to just get year
# first remove extra date_collected column (no idea where this comes from)
bien_raw <- bien_raw[, !duplicated(colnames(bien_raw), fromLast = TRUE)]
bien_raw <- bien_raw %>% separate("date_collected","year",sep="-",remove=T)
sort(unique(bien_raw$year))
# keep only necessary columns & rename to fit standard
bien_raw <- bien_raw %>%
  select("name_matched","verbatim_scientific_name","verbatim_family",
    "identified_by","identification_remarks","date_identified",
    "latitude","longitude","observation_type","year","record_number",
    "locality","county","state_province","country",
    "custodial_institution_codes","dataset","datasource","dataowner",
    "recorded_by","is_cultivated_observation") %>%
  rename("taxon_name" = "name_matched",
         "scientificName" = "verbatim_scientific_name",
         "family" = "verbatim_family",
         "decimalLatitude" = "latitude",
         "decimalLongitude" = "longitude",
         "basisOfRecord" = "observation_type",
         "nativeDatabaseID" = "record_number",
         "stateProvince" = "state_province",
         "institutionCode" = "custodial_institution_codes",
         "datasetName" = "dataset",
         "publisher" = "datasource",
         "rightsHolder" = "dataowner",
         "recordedBy" = "recorded_by")
# add database column
bien_raw$database <- "BIEN"
# combine a few similar columns
bien_raw <- bien_raw %>% unite("taxonIdentificationNotes",
  identified_by:date_identified,na.rm=T, remove=T, sep=" | ")
  bien_raw$taxonIdentificationNotes <-
    gsub("^$", NA, bien_raw$taxonIdentificationNotes)
# fix taxon name notation for forms
bien_raw$taxon_name <- gsub(" fo. "," f. ",bien_raw$taxon_name)
bien_raw$taxon_name <- str_squish(bien_raw$taxon_name)
# create species_name column
bien_raw$species_name <- NA
bien_raw$species_name <- sapply(bien_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(bien_raw$species_name))
# check a few standards and recode if needed
  # basisOfRecord
sort(unique(bien_raw$basisOfRecord))
bien_raw <- bien_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "literature" = "LITERATURE",
    "plot" = "OBSERVATION",
    "specimen" = "PRESERVED_SPECIMEN",
    "trait occurrence" = "OCCURRENCE",
    "occurrence" = "OCCURRENCE",
    "checklist occurrence" = "OCCURRENCE"
  ))
  # establishmentMeans
bien_raw$is_cultivated_observation <- as.character(bien_raw$is_cultivated_observation)
bien_raw <- bien_raw %>%
  mutate(establishmentMeans = recode(is_cultivated_observation,
    "1" = "MANAGED",
    "0" = "UNCERTAIN",
    .missing = "UNCERTAIN")) %>%
  select(-is_cultivated_observation)
# write file
write.csv(bien_raw, file.path(main_dir,data,standard,"bien.csv"),
  row.names=FALSE)
rm(bien_raw)



###
###############
###############################################
# F) USDA Forest Service, Forest Inventory and Analysis (FIA)
#    data: https://apps.fs.usda.gov/fia/datamart/datamart.html
#    metadata doc: https://www.fia.fs.usda.gov/library/database-documentation/current/ver90/FIADB%20User%20Guide%20P2_9-0-1_final.pdf
#    Note that as of late 2022 data for most of Alaska is not available,
#    in addition to District of Columbia.
#    There is also an rFIA package (released 2020) that could be tried- I think
#    currently it's not working due to the recent FIA 2.0 release.
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(main_dir,data,raw,"FIA")))
  dir.create(file.path(main_dir,data,raw,"FIA"),
  recursive=T)

# download and read in supplemental data tables
download.file("https://apps.fs.usda.gov/fia/datamart/CSV/FIADB_REFERENCE.zip",
  file.path(main_dir,data,raw,"FIA","reference_tables"))
  # species codes
unzip(zipfile = file.path(main_dir,data,raw,"FIA","reference_tables"),
      files = "REF_SPECIES.csv",
      exdir = file.path(main_dir,data,raw,"FIA"))
sp_codes <- read.csv(file.path(main_dir,data,raw,"FIA","REF_SPECIES.csv"))
    # create taxon name and select only the columns we need
subsp <- sp_codes %>% filter(SUBSPECIES != "")
subsp$taxon_name <- paste(subsp$GENUS,subsp$SPECIES,"subsp.",subsp$SUBSPECIES)
var <- sp_codes %>% filter(VARIETY != "")
var$taxon_name <- paste(var$GENUS,var$SPECIES,"var.",var$VARIETY)
spp <- sp_codes %>% filter(SUBSPECIES == "" & VARIETY == "")
spp$taxon_name <- paste(spp$GENUS,spp$SPECIES)
sp_codes <- Reduce(rbind.fill,list(subsp,var,spp))
sp_codes <- sp_codes %>% select(SPCD,taxon_name)
    # get FIA species that are in our target taxa list
taxon_fia <- sp_codes[which(sp_codes$taxon_name %in% taxon_names),]
taxon_fia
species_codes <- sort(unique(taxon_fia$SPCD))

# set up list of states and territories to cycle through (57 total)
state_abb <- c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID",
  "IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE",
  "NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN",
  "TX","UT","VT","WV","WA","VA","WI","WY"
    # turn these off if you don't want territories:
  ,"AS","FM","GU","MP","PW","PR","VI"
)

###
### OPTION 1: download FIA data manually (~6.65 GB) and run everything locally.
### this works better if you have enough disk space since it avoids connection
### issues (but, see option 2 below if you don't have disk space)
### NOTE: Option 2 (direct web access) is not working with the new FIA DataMart
### release; files stop downloading at ~50 MB, so most states cannot be
### downloaded in this way (e.g., RI is 10 MB and works; AK is 103 MB and fails)
###

# If you have access to The Morton Arboretum GTCP Shared drive and you don't
# need a more recent version of the data, you can use the files already
# downloaded to the "FIA_data_Nov2022" folder.
# Otherwise, you'll need to download things from scratch :( ...
# First, download the full FIA dataset via this link (6.65 GB !!)
# https://apps.fs.usda.gov/fia/datamart/CSV/CSV_FIADB_ENTIRE.zip
# this can take a few hours or more if you have a slow connection !
# Now unzip the file and move it where you'd like on your computer
# Next, go to the FIA DataMart and download each state tree file manually:
# https://apps.fs.usda.gov/fia/datamart/datamart.html
# To do this, click the "Select States" dropdown arrow and click "Select all"
# Then, scroll through the list and click (download) every TREE file
# (e.g. AL_TREE.csv)... I know this process sucks. You're not alone.
# If not every state shows in the list, you may have to select them one-by-one
# in the dropdown menu.
# Once you've done this for every state, pat yourself on the back and
# move them to a folder called "FIA_TREE"

# If you have access to the Morton Googld Drive, set your working directory
# to this Shared drive folder (change slighly based on your computer)...
fia_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/Global Tree Conservation Program/4. GTCP_Projects/Gap Analyses/FIA_data_Nov2022"
# If you downloaded manually, set the working directory for the location of
# the folder where you placed all the downloaded files...
#fia_dir <- "./Desktop/*work/FIA_data"

# read in metadata tables and select relevant columns
  # plot data (has lat-long information)
fia_plots <- read.csv(file.path(fia_dir,"CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv"))
fia_plots <- fia_plots %>% select(INVYR,STATECD,UNITCD,COUNTYCD,PLOT,LAT,LON)
  # state and county codes and names
county_codes <- read.csv(file.path(fia_dir,"CSV_FIADB_ENTIRE/ENTIRE_COUNTY.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
county_codes <- county_codes %>% select(STATECD,COUNTYCD,COUNTYNM)

# function to extract target species data from each state CSV,
# via a local copy of the data
extract_tree_data_local <- function(state_abb){
  data <- data.frame()
  # read in tree data, which lists all species and the plots they're in;
  # this takes longer for larger states
  cat("Reading in data for",state_abb)
  state_df <- read.csv(file.path(fia_dir,"FIA_TREE",
    paste0(state_abb,"_TREE.csv")))
  # cycle through list of target taxon codes and extract those rows from
  # the state CSV
  for (sp in 1:length(species_codes)){
    target_sp <- state_df[which(state_df$SPCD==species_codes[[sp]]),]
    data <- rbind(data, target_sp)
  }
  # remove state file to make space for reading in the next one
  rm(state_df)
  # take a look at how much data were pulled
  cat("--> ",nrow(data)," observations.<br>")
  # keep only necessary data columns
  data_sm <- data %>% select(SPCD,INVYR,UNITCD,COUNTYCD,PLOT,STATECD,STATUSCD)
  # return reults
  return(data_sm)
  rm(sp)
}

# loop through states and pull data
fia_outputs <- lapply(state_abb, extract_tree_data_local)
length(fia_outputs)
# stack state-by-state data extracted to create one dataframe
fia_raw <- data.frame()
for(file in seq_along(fia_outputs)){
  fia_raw <- rbind(fia_raw, fia_outputs[[file]])
}; nrow(fia_raw)
# join FIA data to supplemental tables
fia_raw <- join(fia_raw,sp_codes)
fia_raw <- join(fia_raw,fia_plots)
fia_raw <- join(fia_raw,county_codes)
# write file of raw data
write.csv(fia_raw,file.path(main_dir,data,raw,"FIA",
  "fia_extracted_local.csv"),row.names=FALSE)

# NOTE !!! NOW SCROLL DOWN TO "STANDARDIZE THE DATA" SECTION !!!

###
### OPTION 2: read each file from the web then remove before moving to the next.
### this helps alleviate disk space concerns but requires stable internet
### connection for a few hours and can be very slow if you don't have a very
### fast connection
### NOTE EBB late 2022: could not download state TREE files :( got this error:
### "download from 'https://apps.fs.usda.gov/fia/datamart/CSV/AL_TREE.csv' failed"
### looks like it downloaded just part of the data but not all
###

# set where you want files to be temporarily downloaded while we extract data
local_dir <- "./Desktop/*work"
# we cant read directly from the web; throws "Transferred a partial file" error

# function to extract target species data from each state CSV,
# via the online FIA DataMart
extract_tree_data_web <- function(state_abb){
  data <- data.frame()
  # read in tree data, which lists all species and the plots they're in;
  # this takes longer for larger states
  cat("Downloading data for",state_abb)
  download.file(
    paste0("https://apps.fs.usda.gov/fia/datamart/CSV/",state_abb,"_TREE.csv"),
    destfile = file.path(local_dir,"fia_state_file.csv"))
  state_df <- read.csv(file.path(local_dir,"fia_state_file.csv"))
  unlink(file.path(local_dir,"fia_state_file.csv"))
  # cycle through list of target taxon codes and extract those rows from
  # the state CSV
  for (sp in 1:length(species_codes)){
    target_sp <- state_df[which(state_df$SPCD==species_codes[[sp]]),]
    data <- rbind(data, target_sp)
  }
  # remove state file to make space for reading in the next one
  rm(state_df)
  # take a look at how much data were pulled
  cat(state_abb,": ",nrow(data)," observations. ")
  # keep only necessary data columns
  data_sm <- data %>% select(SPCD,INVYR,UNITCD,COUNTYCD,PLOT,STATECD,STATUSCD)
  # read in metadata tables and select relevant columns
    # plot data (has lat-long information)
  fia_plots <- read.csv(url(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/",
    state_abb,"_PLOT.csv")))
  fia_plots <- fia_plots %>% select(INVYR,STATECD,UNITCD,COUNTYCD,PLOT,LAT,LON)
    # state and county codes and names
  county_codes <- read.csv(url(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/",
    state_abb,"_COUNTY.csv")))
  county_codes <- county_codes %>% select(STATECD,COUNTYCD,COUNTYNM)
  # join FIA data to supplemental tables
  data_sm <- join(data_sm,fia_plots)
  data_sm <- join(data_sm,county_codes)
  # return reults
  return(data_sm)
  rm(sp)
}

# loop through states and pull data
fia_outputs <- lapply(state_abb, extract_tree_data_web)
length(fia_outputs)
# stack state-by-state data extracted to create one dataframe
fia_raw <- data.frame()
for(file in seq_along(fia_outputs)){
  fia_raw <- rbind(fia_raw, fia_outputs[[file]])
}; nrow(fia_raw)
# join FIA data to species codes
fia_raw <- join(fia_raw,sp_codes)
# write file of raw data
write.csv(fia_raw,file.path(main_dir,data,raw,"FIA",
  "fia_extracted_web.csv"),row.names=FALSE)

###
### STANDARDIZE THE DATA
### Run this section no matter which option your chose for getting the data !!
###

# rename columns to fit standard & create ID column
fia_raw <- fia_raw %>%
  rename("decimalLatitude" = "LAT",
         "decimalLongitude" = "LON",
         "year" = "INVYR",
         "isAlive" = "STATUSCD",
         "county" = "COUNTYNM") %>%
  unite("nativeDatabaseID", c("year","UNITCD","COUNTYCD","PLOT","STATECD"),
      remove=F, sep="-")
 # add a few informational columns
fia_raw$database <- "FIA"
fia_raw$datasetName <- "Forest Inventory and Analysis Database"
fia_raw$rightsHolder <- "USDA Forest Service"
fia_raw$individualCount <- 1
# create species_name column
fia_raw$species_name <- NA
fia_raw$species_name <- sapply(fia_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(fia_raw$species_name))
  # check a few standards and recode if needed
    # establishmentMeans
fia_raw$isAlive <- as.character(fia_raw$isAlive)
fia_raw <- fia_raw %>%
  mutate(establishmentMeans = recode(isAlive,
    "1" = "UNCERTAIN",
    "2" = "DEAD", #CUT
    "3" = "DEAD",
    "0" = "UNCERTAIN"))
fia_raw <- fia_raw %>% select(-isAlive)
    # basisOfRecord
fia_raw$basisOfRecord <- "HUMAN_OBSERVATION"
    # year
fia_raw$year[which(fia_raw$year == 9999)] <- NA
  # keep only necessary columns
fia_raw <- fia_raw %>%
  select(nativeDatabaseID,year,taxon_name,county,decimalLatitude,
    decimalLongitude,database,species_name,establishmentMeans,
    basisOfRecord,datasetName,rightsHolder,individualCount)
# write file
write.csv(fia_raw, file.path(main_dir,data,standard,"fia.csv"),
  row.names=FALSE)
rm(fia_raw)


###
###############
###############################################
# G) Ex situ accession-level data: wild collection locations
# These data are compiled in 2-compile_exsitu_data.R and here we are simply 
# formatting them for use in the occurrence point dataset.
###############################################
###############
###

# read in ex situ data we saved in 2-compile_exsitu_data.R
exsitu_raw1 <- read.csv(file.path(main_dir,data,raw,"Ex-situ",
  "ExSitu_Compiled_Post-Geolocation_2022-12-07.csv"), colClasses = "character",
  na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
exsitu_raw2 <- read.csv(file.path(main_dir,data,raw,"Ex-situ",
  "ExSitu_Dead_2022-12-07.csv"), colClasses = "character",
  na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
exsitu_raw <- rbind.fill(exsitu_raw1,exsitu_raw2)
nrow(exsitu_raw) #9334
str(exsitu_raw)
rm(exsitu_raw1,exsitu_raw2)
# rename columns to fit standard
exsitu_raw <- exsitu_raw %>%
  rename("taxon_name" = "taxon_name_accepted",
         "scientificName" = "taxon_full_name_orig",
         "specificEpithet" = "species",
         "taxonRank" = "infra_rank",
         "infraspecificEpithet" = "infra_name",
         "taxonIdentificationNotes" = "taxon_verif",
         "year" = "coll_year",
         "nativeDatabaseID" = "UID",
         "publisher" = "data_source",
         "rightsHolder" = "inst_short",
         "references" = "acc_num",
         "issue" = "flag",
         "individualCount" = "num_indiv",
         "decimalLatitude" = "lat_dd",
         "decimalLongitude" = "long_dd",
         "coordinateUncertaintyInMeters" = "uncertainty",
         "verbatimLocality" = "all_locality",
         "stateProvince" = "state",
         "establishmentMeans" = "prov_type") %>%
  unite("recordedBy", c("coll_name","coll_num"), remove=T, sep="; ") %>%
  unite("geolocationNotes", c("gps_det","geolocated_by",
                              "gps_notes"), remove=T, sep="; ") %>%
  unite("locationNotes", c("germ_type","garden_loc"), remove=T, sep="; ")
# recode if needed establishmentMeans
sort(unique(exsitu_raw$establishmentMeans))
exsitu_raw <- exsitu_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "H" = "CULTIVATED",
    "H?" = "CULTIVATED",
    "N" = "NATIVE",
    "NG" = "UNCERTAIN",
    "U" = "UNCERTAIN",
    "W" = "NATIVE",
    "Z" = "NATIVE"
  ))
# add a few informational columns
exsitu_raw$database <- "Ex_situ"
exsitu_raw$basisOfRecord <- "HUMAN_OBSERVATION"
exsitu_raw$datasetName <- exsitu_raw$rightsHolder
# keep only necessary columns
exsitu_raw <- exsitu_raw %>%
  select(database,taxon_name,scientificName,genus,specificEpithet,taxonRank,
         infraspecificEpithet,taxonIdentificationNotes,year,basisOfRecord,
         nativeDatabaseID,datasetName,publisher,rightsHolder,references,
         issue,recordedBy,establishmentMeans,individualCount,decimalLatitude,
         decimalLongitude,coordinateUncertaintyInMeters,geolocationNotes,
         locality,verbatimLocality,locationNotes,municipality,county,
         stateProvince,country)
head(exsitu_raw)
# write file
write.csv(exsitu_raw, file.path(main_dir,data,standard,"exsitu.csv"),
          row.names=FALSE)
rm(exsitu_raw)

