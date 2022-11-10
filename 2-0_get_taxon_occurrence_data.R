################################################################################

## 2-0_get_occurrence_points.R
### Authors: Emily Beckman & Shannon Still ### Date: 02/05/2020

### DESCRIPTION:
  # This script provides instructions and code chunks for downloading and
  #   standardizing wild occurrence points from:
    # GLOBAL DATABASES (though all likely have U.S. bias?)
      # Global Biodiversity Information Facility (GBIF)
      # Integrated Digitized Biocollections (iDigBio)
      # IUCN Red List of Threatened Species
      # U.S. Herbarium Consortia (SERNEC, SEINet, etc.)
      # Botanical Information and Ecology Network (BIEN)
    # NATIONAL DATABASES
      # Forest Inventory and Analysis (FIA) Program of the USDA Forest Service
      # Biodiversity Information Serving Our Nation (BISON), USGS
  ## NOTE: Not all data from these sources are reliable. The aim of this
  #        script is to get all easily-downloadable occurrence data, which
  #        can then be sorted and vetted for the user's specific purposes.
  ## NOTE: You can add other occurrence point data (e.g., expert comment,
  #        NatureServe, floras, USDA PLANTS, BONAP, IUCN Red List, private
  #        sources, etc.) by standardizing column names and formatting to match
  #        the schema in the "Renaming Columns" tab:
  #        https://docs.google.com/spreadsheets/d/1dllfDXaZBLvB1AsrY1wDS-sPceKAdOY681bqUbfoQAs/edit?usp=sharing
  #        then save as CSV and place in "inputs/compiled_occurrence" folder.

### DATA IN:
  # (optional) target_taxa_with_syn.csv
    # columns:
      # 1. "taxon_name" (genus, species, infra rank, and infra name, all
      #    separated by one space each; hybrid symbol should be " x ", rather
      #    than "_" or "✕", and go between genus and species)
      # 2. (optional) other data you want to keep with taxa info
  # files needed for FIA data download/standardization:
  #   FIA_AppendixF_TreeSpeciesCodes_2016.csv
  #   US_state_county_FIPS_codes.csv

### DATA OUT:
  # gbif.csv
  # idigbio.csv
  # herbaria.csv
  # bien.csv
  # fia.csv
  # bison.csv

################################################################################
# Load libraries
################################################################################

rm(list=ls())
my.packages <- c('plyr', 'tidyverse', 'spocc', 'rgbif', 'data.table', 'BIEN',
                 'ridigbio', 'batchtools', 'googledrive', 'textclean','rbison',
                 'tools')
#  install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/occurrence_points"
#script_dir <- "./Documents/GitHub/OccurrencePoints/scripts"
#log_loc <- "./Desktop/IMLS_passwords.txt"

# or use 0-1_set_workingdirectory.R script:
 source("./Documents/GitHub/OccurrencePoints/scripts/0-1_set_workingdirectory.R")
#source('scripts/0-1_set_workingdirectory.R')

################################################################################
# Load functions
################################################################################
#source(file.path(script_dir,"0-2_load_IMLS_functions.R"))

# calculates percent of each data frame column that is not NA
#percent.filled <- function(df){
#  for(i in 1:ncol(df)){
#    print(paste(names(df)[i],": ",
#                round((nrow(df)-sum(is.na(df[,i])))/nrow(df),3)*100,"%",sep=""))
#  }
#}


################################################################################
################################################################################
# 1. Load/create target taxa list
################################################################################

## IF YOU HAVE CSV OF TARGET TAXA AND SYNONYMS:
# read in taxa list
taxon_list <- read.csv(file.path(main_dir, "inputs","taxa_list",
  "target_taxa_with_syn.csv"), header = T, colClasses="character")
head(taxon_list)
# as needed, filter out cultivars and blank rows
#taxon_list <- taxon_list %>%
#  filter(taxon_type != "cultivar" & !is.na(taxon_name))
nrow(taxon_list) #805 ; 257
# list of target taxon names
taxon_names <- taxon_list$taxon_name

## IF JUST ONE OR A FEW TAXA, CREATE A LIST BY HAND:
#taxon_names <- c("Quercus havardii")

################################################################################
# Shannon working here 4/21/2020
################################################################################

## We would like to do this if time permitting ##
  ## this section is to get the header column names from a file, then subset
  ## the table by the column names wanted, then set the colnames appropriate
  ## for that data.frame
  ## this can be repeated for each new data source
  ## we may have to move some of these steps around
  # b.h <- read_excel('/Users/sstill/Box/Seed Menus Project/HerbariumRecords/HerbariumRecords_wHeaders/Burke_AllSpp_wHeaders.xlsx',
  #                    col_types = 'text')
  ## this line sets the correct data.frame
    # n.h <- 'burke'
  ## this line selects the correct columns
    # nms  <- nm.set %>% filter(!is.na(!!as.name(n.h))) %>%
    #                    select(!!as.name(n.h)) %>% as.character()
    # # nms1 <- nm.set %>% filter(!is.na(!!as.name(n.h))) %>%
    #                      select(final_nm, !!as.name(n.h))
  # b.h <- b.h %>% select_(nms)
  ## this line sets the correct column names
    # names(b.h) <- nms1$final_nm
  ## these lines add some data to fill in where missing
    # b.h <- b.h %>% mutate(dataSource=paste0('burke', sprintf("%05d",
    #   1:nrow(b.h))), id=as.character(id), coordinateSystem=NA_character_,
    #   zone = NA_character_, x_coord = as.numeric(x_coord),
    #   y_coord = as.numeric(y_coord),
    #   coordinateUncertaintyMeters = as.numeric(coordinateUncertaintyMeters))
    # rm(n.h, nms, nms1)
  # b.h <- b.h %>% mutate(EPSG='undefined', data_type='herbarium voucher')

################################################################################
# 2. Download/compile data from each target database
################################################################################

# NOTE: you can pick and choose any/all sections below, depending on which
#   databases you'd like to use

# create folder for raw data
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence"), recursive=T)
# create folder for compiled data
if(!dir.exists(file.path(main_dir,"inputs","compiled_occurrence")))
  dir.create(file.path(main_dir,"inputs","compiled_occurrence"), recursive=T)

###############
### A) Global Biodiversity Information Facility (GBIF)
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","gbif_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","gbif_raw"),
  recursive=T)

# GBIF account user information
  # if you don't have account yet, go to https://www.gbif.org then click
  #   "Login" in top right corner, then click "Register"
# either read in a text file with username, password, and email (one on each
#   line) or igore that line and manually fill in each below:
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
                #pred_in("basisOfRecord", c("PRESERVED_SPECIMEN",
                #    "HUMAN_OBSERVATION","FOSSIL_SPECIMEN","OBSERVATION",
                #    "UNKNOWN","MACHINE_OBSERVATION","MATERIAL_SAMPLE",
                #    "LITERATURE")),
                #pred("hasCoordinate", TRUE),
                #pred("hasGeospatialIssue", FALSE),
                format = "DWCA", #"SIMPLE_CSV"
                user=user,pwd=pwd,
                email=email)
    rm(user, pwd, email)
# load gbif data just downloaded
  # download and unzip before reading in
download_key <- gbif_download
  # must wait for download to complete before continuing;
  # it may take a while (up to 3 hours) if you have a large taxa list;
  # function below will pause script until the download is ready
occ_download_wait(download_key, status_ping=10, quiet=TRUE)
  # get download when its ready then unzip and read in
    ## ESB: there was a new error here 5/26/21 and I couldn't find any info on it...
    ##    "Error in occ_download_get(key = download_key[1], path = file.path(main_dir,  :
    ##    response content-type != application/octet-stream; q=0.5"
    ## Ended up manually downloading from GBIF website then reading into R (line 214)
occ_download_get(key=download_key[1],path=file.path(main_dir,"inputs",
  "raw_occurrence","gbif_raw"),overwrite=TRUE) #Download file size: 684.49 MB
unzip(zipfile=paste0(file.path(main_dir,"inputs","raw_occurrence","gbif_raw",
  download_key[1]),".zip"), files="occurrence.txt",
  exdir=file.path(main_dir,"inputs","raw_occurrence","gbif_raw"))
  # read in data
gbif_raw <- fread(file.path(main_dir,"inputs","raw_occurrence","gbif_raw",
  "occurrence.txt"), quote="", na.strings="")
nrow(gbif_raw) #3073391

### standardize column names

# create taxon_name column
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

# keep only necessary columns
gbif_raw <- gbif_raw %>% dplyr::select(
    # taxon name
  "taxon_name","scientificName",
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

# recode standards
  # establishmentMeans
sort(unique(gbif_raw$establishmentMeans))
gbif_raw <- gbif_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "Uncertain" = "UNKNOWN",
    "Native" = "NATIVE"))

# fix taxa names
gbif_raw$taxon_name <- mgsub(gbif_raw$taxon_name,
  c("Tilia xeuropaea","Tilia xvulgaris"),
  c("Tilia x europaea","Tilia x vulgaris"))

# create species_name column
gbif_raw$species_name <- NA
gbif_raw$species_name <- sapply(gbif_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(gbif_raw$species_name))

# check data
#percent.filled(gbif_raw)
head(gbif_raw)

# write file
write.csv(gbif_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "gbif.csv"),row.names=FALSE)
#delete files no longer needed (large files)
  #unlink(paste0(file.path(main_dir,"inputs","raw_occurrence","gbif_raw",
  #          "occurrence.txt")))
  # unlink(paste0(file.path(main_dir,"inputs","raw_occurrence","gbif_raw",
  #           download_key[1]),".zip"))
rm(gbif_raw)

###############
# B) Integrated Digitized Biocollections (iDigBio)
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","idigbio_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","idigbio_raw"),
  recursive=T)

# The R interface doesn't get all fields available; use manual
#    method down below (line 341) if you want to be sure

# download iDigBio data for target taxa
  # we have to go taxon by taxon; function can only return 100,000
  #   records at once and Quercus has more than that so can't download by genera
idigbio_raw <- data.frame()
for(i in 1:length(taxon_names)){
  output_new <- idig_search_records(rq=list(scientificname=taxon_names[[i]]),
    fields="all")
  idigbio_raw <- rbind(idigbio_raw,output_new)
  print(paste(round(i/length(taxon_names)*100,digits=1),"% complete",sep=""))
}
  rm(output_new, i)
nrow(idigbio_raw) #150409
# remove rows that are lists
idigbio_raw <- idigbio_raw %>% dplyr::select(everything(),-commonnames,-flags,
  -mediarecords,-recordids)
# write file
write.csv(idigbio_raw,file.path(main_dir,"inputs","raw_occurrence",
  "idigbio_raw","idigbio_R_download.csv"),row.names=FALSE)

# split date collected column to just get year
idigbio_raw <- idigbio_raw %>% separate("eventdate","year",sep="-",remove=T)
idigbio_raw$year <- gsub("[[:digit:]]+/[[:digit:]]+/","",idigbio_raw$year)
  sort(unique(idigbio_raw$year))

# keep only necessary columns
idigbio_raw$taxon_name <- idigbio_raw$scientificname
idigbio_raw <- idigbio_raw %>% dplyr::select(
  "taxon_name","scientificname",
  "family","genus","specificepithet","taxonrank","infraspecificepithet",
  "geopoint.lon","geopoint.lat","coordinateuncertainty",
  "basisofrecord","year","uuid","occurrenceid",
  "locality","verbatimlocality","county","municipality","stateprovince",
    "country","countrycode",
  "institutioncode","collectioncode"
)

# rename columns
setnames(idigbio_raw,
  old = c("scientificname",
          "specificepithet","taxonrank","infraspecificepithet",
          "geopoint.lon","geopoint.lat",
          "coordinateuncertainty",
          "basisofrecord","uuid","occurrenceid",
          "verbatimlocality","stateprovince","countrycode",
          "institutioncode","collectioncode"),
  new = c("scientificName",
          "specificEpithet","taxonRank","infraspecificEpithet",
          "decimalLongitude","decimalLatitude",
          "coordinateUncertaintyInMeters",
          "basisOfRecord","nativeDatabaseID","references",#"recordNumber",
          "verbatimLocality","stateProvince","countryCode",
          "institutionCode","publisher"),
  skip_absent=T)
idigbio_raw$database <- "iDigBio"
idigbio_raw$datasetName <- idigbio_raw$institutionCode
idigbio_raw$rightsHolder <- idigbio_raw$institutionCode

# fix taxa names
idigbio_raw$taxon_name <- str_to_sentence(idigbio_raw$taxon_name)
idigbio_raw$family <- str_to_title(idigbio_raw$family)
idigbio_raw$genus <- str_to_title(idigbio_raw$genus)

############
# MANUAL WAY:
## First, download raw data
## Go to https://www.idigbio.org/portal/search
## Click "Add a field" dropdown on the left and select "Genus"; type your
##   target genus name into the "Genus" box
## Click the "Download" tab, type in your email, and click the download button
##   (down arrow within circle)
## If you have more than one target genus, repeat the above steps for the
##   other genera
## Your downloads will pop up in the "Downloads" section;
##   "Click To Download" for each
## Move all the folders you downloaded into the "idigbio_raw" folder
## Pull the "occurrence_raw.csv" file out into the
##   "idigbio_raw" folder and rename with appropriate genus name
## read in raw occurrence points
#file_list <- list.files(path = file.path(main_dir,"inputs","raw_occurrence",
#  "idigbio_raw"),pattern = "occurrence", full.names = T)
#file_dfs <- lapply(file_list, read.csv, colClasses = "character",
#  na.strings=c("","NA"),strip.white=T,fileEncoding="UTF-8")
#length(file_dfs) #4
# stack datasets to create one dataframe
#idigbio_raw <- data.frame()
#for(file in seq_along(file_dfs)){
#  idigbio_raw <- rbind(idigbio_raw, file_dfs[[file]])
#}; nrow(idigbio_raw) #326242
# replace prefixes in column names
#colnames(idigbio_raw) <- gsub(x = colnames(idigbio_raw),
#  pattern = "dwc\\.", replacement = "")
### standardize column names
## create taxon_name column
#sort(unique(idigbio_raw$taxonRank))
#subspecies <- c("subsp.","subspecies","Subspecies")
#variety <- c("var.","Varietas","variety","Variety")
#forma <- c("f.","fm.","form","Form","form.","forma","Forma","Subform")
#species <- c("espécie","sp.","specie","species","Species","spp.")
#hybrid <- c("×","x","X")
#aff <- c("aff.")
#subsp <- idigbio_raw %>% filter(taxonRank %in% subspecies)
#  subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
#    subsp$infraspecificEpithet)
#var <- idigbio_raw %>% filter(taxonRank %in% variety)
#  var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
#    var$infraspecificEpithet)
#form <- idigbio_raw %>% filter(taxonRank %in% forma)
#  form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
#    form$infraspecificEpithet)
#spp <- idigbio_raw %>% filter(taxonRank %in% species)
#  spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
#h <- idigbio_raw %>% filter(taxonRank %in% hybrid)
#  h$taxon_name <- paste(h$genus,"x",h$specificEpithet)
#a <- idigbio_raw %>% filter(taxonRank %in% aff)
#  a$taxon_name <- paste(a$genus,"aff.",a$specificEpithet)
#idigbio_raw <- Reduce(rbind.fill,list(subsp,var,form,spp,h,a))
## combine a few similar columns
#idigbio_raw <- idigbio_raw %>% unite("taxonIdentificationNotes",
#    c(identificationID:identifiedBy,taxonRemarks),na.rm=T,remove=T,sep=" | ")
#  idigbio_raw$taxonIdentificationNotes <-
#    gsub("^$",NA,idigbio_raw$taxonIdentificationNotes)
#idigbio_raw <- idigbio_raw %>% unite("locationNotes",
#  c(associatedTaxa,eventRemarks,fieldNotes,habitat,locationRemarks,
#    occurrenceRemarks,occurrenceStatus),na.rm=T,remove=T,sep=" | ")
#  idigbio_raw$locationNotes <- gsub("^$",NA,idigbio_raw$locationNotes)
#idigbio_raw <- idigbio_raw %>% unite("geolocationNotes",
#  georeferenceProtocol:georeferencedDate,na.rm=T,remove=T,sep=" | ")
#  idigbio_raw$geolocationNotes <- gsub("^$",NA,idigbio_raw$geolocationNotes)
## split date collected column to just get year
#idigbio_raw <- idigbio_raw %>% separate("eventDate","year",sep="-",remove=T)
#idigbio_raw$year <- gsub("[[:digit:]]+/[[:digit:]]+/","",idigbio_raw$year)
#idigbio_raw$year <- gsub("^[1-9][0-9][0-9][0-9]/","",idigbio_raw$year)
#idigbio_raw$year <- gsub("^[0-9][0-9]/","",idigbio_raw$year)
#idigbio_raw$year <- gsub("/","",idigbio_raw$year)
#idigbio_raw$year <- as.numeric(idigbio_raw$year)
# keep only years in reasonable timeframe (1500-current year)
#idigbio_raw$year[which(idigbio_raw$year < 1500 |
#  idigbio_raw$year > as.numeric(format(Sys.time(),"%Y")))] <- NA
#  sort(unique(idigbio_raw$year))
## keep only necessary columns
#idigbio_raw <- idigbio_raw %>% dplyr::select(
#    "taxon_name","scientificName","taxonIdentificationNotes",
#    "family","genus","specificEpithet","taxonRank","infraspecificEpithet",
#    "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
#    "basisOfRecord","establishmentMeans","year",
#    "locality","verbatimLocality","locationNotes","county","municipality",
#    "higherGeography","stateProvince","country","countryCode",
#    "institutionCode","datasetName","rightsHolder","dcterms.license","geolocationNotes",
#    "recordedBy","coreid","dcterms.references","informationWithheld") %>%
#  rename(license = dcterms.license, references = dcterms.references)
#idigbio_raw$database <- "iDigBio"
## recode standard column: establishment means
#idigbio_raw$establishmentMeans <- str_to_upper(idigbio_raw$establishmentMeans)
#unique(idigbio_raw$establishmentMeans)
#idigbio_raw <- idigbio_raw %>%
#  mutate(establishmentMeans = recode(establishmentMeans,
#    "CULTIVATED" = "MANAGED",
#    "LAWN DECORATION" = "MANAGED",
#    "CAPTIVE" = "MANAGED",
#    "PLANTED" = "MANAGED",
#    "MANAGED OR ESCAPED" = "MANAGED",
#    "ESCAPED" = "INTRODUCED",
#    "WILD." = "NATIVE",
#    "CULTIVATED." = "MANAGED",
#    "ALIEN" = "INVASIVE",
#    "SHADE TREE PLANTED" = "MANAGED",
#    .default = "UNKNOWN",
#    .missing = "UNKNOWN"))
############

# create species_name column
idigbio_raw$species_name <- NA
idigbio_raw$species_name <- sapply(idigbio_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(idigbio_raw$species_name))

# recode standard columns
  # basis of record
idigbio_raw <- idigbio_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "machineobservation" = "MACHINE_OBSERVATION",
    "humanobservation" = "HUMAN_OBSERVATION",
    "fossilspecimen" = "FOSSIL_SPECIMEN",
    .missing = "UNKNOWN"))

# check data
#percent.filled(idigbio_raw)
head(idigbio_raw)

# write file
write.csv(idigbio_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "idigbio.csv"),row.names=FALSE)

rm(idigbio_raw)

###############
# C) IUCN Red List
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","redlist_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","redlist_raw"),
  recursive=T)

# FIRST, download raw data:
  # If you don't have an IUCN RL account yet, create one, and sign in
  # Go to https://www.iucnredlist.org/search
  # Open the "Taxonomy" tab on the left, and type in your target genus name
  #   and check the box next to the genus name when it comes up;
  #   or, alternatively, if you are just looking for a few
  #   taxa you can search for them individually
  # You should be able to add each genus/taxon to your search so only
  #   only file needs to be exported
  # In the far-left bar, scroll down and, if desired, check
  #   "Subspecies and varieties"
  # Click the grey "Download" button and select "Range data - Points (CSV)";
  #   then fill in the popup window information
  # Go to https://www.iucnredlist.org/account to find your query
  # Click "Download" next to your query
  # Move the folder you downloaded into the "redlist_raw" folder
  # Pull the "points_data.csv" file out into the
  #   "redlist_raw" folder

# read in raw occurrence points
redlist_raw <- read.csv(file.path(main_dir,"inputs","raw_occurrence","redlist_raw",
  "points_data.csv"),
  colClasses = "character", na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
nrow(redlist_raw) #161765

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
  dplyr::select(taxon_name,binomial,tax_comm,event_year,
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

# check data
#percent.filled(redlist_raw)
head(redlist_raw)

# write file
write.csv(redlist_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "redlist.csv"),row.names=FALSE)

rm(redlist_raw)

###############
# D) U.S. Herbaria Consortia (SERNEC, SEINet, etc.)
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","sernec_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","sernec_raw"),
  recursive=T)

# FIRST, download raw data:
  # Go to http://sernecportal.org/portal/collections/harvestparams.php
  # Type your target genus name into the "scientific name" box and click
  #   "List Display"; or, alternatively, if you are just looking for a few
  #   taxa you can search for and download them individually
  # Click the Download Specimen Data button (arrow pointing down into a box),
  #   in the top right corner
  # In the pop-up window, select the "Darwin Core" radio button,
  #   uncheck everything in the "Data Extensions" section, and
  #   select the "UTF-8 (unicode)" radio button
  #   leave other fields as-is
  # Click "Download Data"
  # If you have more than one target genus, repeat the above steps for the
  #   other genera
  # Move all the folders you downloaded into the "sernec_raw" folder

# pulls the "occurrences.csv" from each genus folder, for compilation
file_list <- list.files(file.path(main_dir,"inputs","raw_occurrence","sernec_raw"),
  pattern = "SymbOutput", full.names = T)
file_dfs <- lapply(file_list, function(i){
  read.csv(file.path(i,"occurrences.csv"), colClasses = "character",
   na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")})
sernec_raw <- data.frame()
  for(file in seq_along(file_dfs)){
    sernec_raw <- rbind(sernec_raw, file_dfs[[file]])
  }
nrow(sernec_raw) #235775

### standardize column names

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

# keep only necessary columns
sernec_raw <- sernec_raw %>% dplyr::select(
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
sernec_raw$database <- "US_Herbaria"
sernec_raw$datasetName <- sernec_raw$institutionCode

# rename columns
sernec_raw <- sernec_raw %>% rename(nativeDatabaseID = id)

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

# recode standard columns
  # basis of record
sernec_raw <- sernec_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "PreservedSpecimen" = "PRESERVED_SPECIMEN",
    "Preserved specimen" = "PRESERVED_SPECIMEN",
    "Preserved Specimen" = "PRESERVED_SPECIMEN",
    "Physicalspecimen" = "PRESERVED_SPECIMEN",
    "preservedspecimen" = "PRESERVED_SPECIMEN",
    "preservedSpecimen" = "PRESERVED_SPECIMEN",
    "Observation" = "OBSERVATION",
    "LivingSpecimen" = "LIVING_SPECIMEN",
    "Physical specimen" = "PRESERVED_SPECIMEN",
    "Ejemplar herborizado" = "PRESERVED_SPECIMEN",
    "Physical Specimen" = "PRESERVED_SPECIMEN",
    "HumanObservation" = "HUMAN_OBSERVATION",
    .default = "UNKNOWN"))
  # year
sort(unique(sernec_raw$year))
sernec_raw <- sernec_raw %>%
  mutate(year = recode(year,
    "9999" = "0",
    "18914" = "1891",
    "19418" = "1941"))
sernec_raw$year <- as.integer(sernec_raw$year)
sernec_raw$year[which(sernec_raw$year < 1500)] <- NA
sort(unique(sernec_raw$year))
  # establishment means
sort(unique(sernec_raw$establishmentMeans))
sernec_raw <- sernec_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
    "Alien" = "INTRODUCED",
    "clonal" = "UNKNOWN",
    "Native" = "NATIVE",
    "Native." = "NATIVE",
    "native" = "NATIVE",
    "Wild." = "NATIVE",
    "Naturalized." = "INTRODUCED",
    "Uncertain" = "UNKNOWN",
    "wild caught" = "UNKNOWN",
    .default = "MANAGED"))

# check data
#percent.filled(sernec_raw)
head(sernec_raw)

# write file
write.csv(sernec_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "sernec.csv"),row.names=FALSE)

rm(sernec_raw)

###############
# E) Botanical Information and Ecology Network (BIEN)
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","bien_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","bien_raw"),
  recursive=T)

# information about functions in package
#vignette("BIEN")

# download BIEN occurrence data for target taxa
bien_raw <- BIEN_occurrence_species(taxon_names,all.taxonomy=T,native.status=T,
  natives.only=F,observation.type=T,collection.info=T,political.boundaries=T,
  cultivated=T)
nrow(bien_raw) #2210103
write.csv(bien_raw, file.path(main_dir,"inputs","raw_occurrence","bien_raw",
  "bien_R_download.csv"),row.names=FALSE)

# get citation info and write file
bien_citation <- BIEN_metadata_citation(dataframe = bien_raw)
write.csv(bien_citation, file.path(main_dir,"inputs","raw_occurrence","bien_raw",
  "bien_data_citation.csv"),row.names=FALSE)

### standardize column names

# split date collected column to just get year
  # first remove extra date_collected column (no idea where this comes from)
    bien_raw <- bien_raw[, !duplicated(colnames(bien_raw), fromLast = TRUE)] #bien_raw <- bien_raw[,-24]
    bien_raw <- bien_raw %>% separate("date_collected","year",sep="-",remove=T)
  sort(unique(bien_raw$year))

# keep only necessary columns
bien_raw <- bien_raw %>% dplyr::select(
  "name_matched","verbatim_scientific_name","verbatim_family",
  "identified_by","identification_remarks","date_identified",
  "latitude","longitude",
  "observation_type","year","record_number",
  "locality","county","state_province","country",
  "custodial_institution_codes","dataset","datasource","dataowner","recorded_by",
  "is_cultivated_observation")

## On a computer at UC Davis, I get an error message and cannot get data and I
  ##    am unable to connect to the correct port:
  # Here is the error message:
  #   Error in postgresqlNewConnection(drv, ...) :
  #   RS-DBI driver: (could not connect public_bien@vegbiendev.nceas.ucsb.edu:5432 on dbname "public_vegbien": could not connect to server: Connection timed out (0x0000274C/10060)
  #                   Is the server running on host "vegbiendev.nceas.ucsb.edu" (128.111.85.31) and accepting
  #                   TCP/IP connections on port 5432?

## This is the line about the issue on the vignette("BIEN"):
  ##Database connection issues Some institution and computer programs
  ##  (e.g. some antivirus programs) block the SQL connections that this package
  ##  relies on. While we are exploring ways around this issue, at present the
  ##  simplest method is to use the package on a computer/network that doesn’t
  ##  block SQL connections.

## This will not be fixed in the immediate future. If this deos not work in an
  ##  IT setting that bloacks the specific port, then simply have to run from a
  ## computer outside that system.

# rename columns
setnames(bien_raw,
  old = c("name_matched","verbatim_scientific_name",
          "verbatim_family",
          "latitude","longitude",
          "observation_type","record_number",
          "state_province",
          "custodial_institution_codes",
          "dataset","datasource",
          "dataowner","recorded_by"),
  new = c("taxon_name","scientificName",
          "family",
          "decimalLatitude","decimalLongitude",
          "basisOfRecord","nativeDatabaseID",
          "stateProvince",
          "institutionCode",
          "datasetName","publisher",
          "rightsHolder","recordedBy"),
  skip_absent=T)
bien_raw$database <- "BIEN"

# OPTIONAL: remove rows from GBIF or FIA if you are separately downloading those datasets
bien_raw <- bien_raw %>% filter(publisher != "FIA" & publisher != "GBIF")
nrow(bien_raw) #70788

# combine a few similar columns
bien_raw <- bien_raw %>% unite("taxonIdentificationNotes",
  identified_by:date_identified,na.rm=T, remove=T, sep=" | ")
  bien_raw$taxonIdentificationNotes <-
    gsub("^$", NA, bien_raw$taxonIdentificationNotes)

# fix taxa names
bien_raw$taxon_name <- gsub(" fo. "," f. ",bien_raw$taxon_name)
bien_raw$taxon_name <- str_squish(bien_raw$taxon_name)

# create species_name column
bien_raw$species_name <- NA
bien_raw$species_name <- sapply(bien_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(bien_raw$species_name))

# recode standard columns
  # basis of record
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
  # establishment means
bien_raw$is_cultivated_observation <- as.character(
  bien_raw$is_cultivated_observation)
bien_raw <- bien_raw %>%
  mutate(establishmentMeans = recode(is_cultivated_observation,
    "1" = "MANAGED",
    "0" = "UNKNOWN",
    .missing = "UNKNOWN")) %>%
  dplyr::select(-is_cultivated_observation)

# check data
#percent.filled(bien_raw)
head(bien_raw)

# write file
write.csv(bien_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "bien.csv"),row.names=FALSE)

rm(bien_raw)

###############
# F) USDA Forest Service, Forest Inventory and Analysis (FIA)
###############

# create new folders if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","fia_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","fia_raw"),
  recursive=T)
if(!dir.exists(file.path(main_dir,"inputs","fia_tables")))
  dir.create(file.path(main_dir,"inputs","fia_tables"),
  recursive=T)

# download and read in supplemental tables
  # plot data (has lat-long information)
download.file("https://apps.fs.usda.gov/fia/datamart/CSV/PLOT.csv",
  file.path(main_dir,"inputs","fia_tables","PLOT.csv"))
fia_plots <- read.csv(file.path(main_dir,"inputs","fia_tables","PLOT.csv"))
    # remove unnecessary columns from plot data
fia_plots <- fia_plots[,c("INVYR","STATECD","UNITCD","COUNTYCD","PLOT",
  "LAT","LON")]
  # state and county codes and names
county_codes <- read.csv(file.path(main_dir,"inputs","fia_tables",
  "US_state_county_FIPS_codes.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
  # read in FIA species codes
fia_codes <- read.csv(file.path(main_dir,"inputs","fia_tables",
  "FIA_AppendixF_TreeSpeciesCodes_2016.csv"),colClasses="character")
  # join taxa list to FIA species codes
fia_codes <- fia_codes[,c(1,3)]
names(fia_codes) <- c("SPCD","taxon_name")
  glimpse(fia_codes)
taxon_fia <- fia_codes[which(fia_codes$taxon_name %in% taxon_names),]
  # make a list of unique FIA species codes to select from the data
species_codes <- sort(unique(taxon_fia$SPCD))
  # check results
sort(unique(taxon_fia$taxon_name[which(!is.na(taxon_fia$SPCD))]))
length(species_codes) #56 ; 54

# You can either read in a CSV file for each state or download all data in
#   one file ("https://apps.fs.usda.gov/fia/datamart/CSV/TREE.csv"; ~10 GB);
#   you need lots of memory to pull out data from the one big file,
#   so the script below uses the state files and cycles through each

# list of states and territories to cycle through
state_abb <- c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID",
  "IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE",
  "NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN",
  "TX","UT","VT","WV","WA","VA","WI","WY"
    # turn these off if you don't want territories:
  ,"AS","FM","GU","MP","PW","PR","VI"
)
# function to extract target species data from each state CSV from FIA
extract_tree_data <- function(state_abb){
  data <- data.frame()
  # read in tree data, which lists all species and the plots they're in; larger
  #   ones will take time to read in
      ## read from online
  #state_df <- read.csv(url(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/",
  #  state_abb,"_TREE.csv")))
      ## read from folder
  state_df <- read.csv(file.path(main_dir,"inputs","raw_occurrence","fia_raw",
    "state_data_files",paste0(state_abb,"_TREE.csv")))
  # cycle through vector of target species codes and extract those rows from
  #   the state CSV
  for (sp in 1:length(species_codes)){
    target_sp <- state_df[which(state_df$SPCD==species_codes[[sp]]),]
    data <- rbind(data, target_sp)
  }
  # remove state file to make space for reading in next one
  rm(state_df)
  # take a look at how much data were pulled
  cat(state_abb,": ",nrow(data)," observations. ")#,
        #grep(file_name, file_list), " of ", length(file_list), ".")
   #print(paste(nrow(data), basename(file_name)))
  # keep only necessary data columns
  data_sm <- data %>% dplyr::select(
    "SPCD","INVYR","UNITCD","COUNTYCD","PLOT","STATECD","STATUSCD")
  return(data_sm)
  rm(sp)
}
# create list of state files to cycle through
#file_list <- list.files(file.path(main_dir,"inputs","raw_occurrence","fia_raw"),
#  pattern = "TREE", full.names = T)

# loop through states and pull data
# NOTE: you need a good internet connection for this if files are not local
  fia_outputs <- lapply(state_abb, extract_tree_data)
    length(fia_outputs) #57

# stack state-by-state data extracted to create one dataframe
fia_raw <- data.frame()
for(file in seq_along(fia_outputs)){
  fia_raw <- rbind(fia_raw, fia_outputs[[file]])
}
nrow(fia_raw) #3444367
# save(fia_raw, file=file.path(main_dir, "inputs", ""))

# write file of raw data
write.csv(fia_raw,file.path(main_dir,"inputs","raw_occurrence",
  "fia_raw","fia_extracted_raw.csv"),row.names=FALSE)

### standardize column names
# join FIA data to supplemental tables
fia_raw <- join(fia_raw,fia_codes)
fia_raw <- join(fia_raw,county_codes)
fia_raw <- join(fia_raw,fia_plots)
# create ID column
fia_raw <- fia_raw %>% unite("fiaPlotID",
  c("INVYR","UNITCD","COUNTYCD","PLOT","STATECD"),remove=F,sep="-")

# rename columns
setnames(fia_raw,
  old = c("LAT","LON",
          "INVYR",
          "STATUSCD","fiaPlotID"),
  new = c("decimalLatitude","decimalLongitude",
          "year",
          "isAlive","nativeDatabaseID"),
  skip_absent=T)
fia_raw$database <- "FIA"
fia_raw$datasetName <- "Forest Inventory and Analysis Database"
fia_raw$rightsHolder <- "USDA Forest Service"

# create species_name column
fia_raw$species_name <- NA
fia_raw$species_name <- sapply(fia_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(fia_raw$species_name))

# recode standard columns
  # establishment means
fia_raw$isAlive <- as.character(fia_raw$isAlive)
fia_raw <- fia_raw %>%
  mutate(establishmentMeans = recode(isAlive,
    "1" = "UNKNOWN",
    "2" = "CUT",
    "3" = "DEAD",
    "0" = "UNKNOWN"))
fia_raw <- fia_raw %>% dplyr::select(-isAlive)
  # basis of record
fia_raw$basisOfRecord <- "OBSERVATION"
  # year
fia_raw$year[which(fia_raw$year == 9999)] <- NA

# remove extra columns
fia_raw <- fia_raw %>%
  dplyr::select(nativeDatabaseID,year,taxon_name,county,stateProvince,
    decimalLatitude,decimalLongitude,database,species_name,establishmentMeans,
    basisOfRecord,datasetName,rightsHolder)

# check data
#percent.filled(fia_raw)
head(fia_raw)

# write file
write.csv(fia_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "fia.csv"), row.names=FALSE)
rm(fia_raw)

###############
# G) Biodiversity Information Serving Our Nation (BISON), USGS
###############

# create new folder if not already present
if(!dir.exists(file.path(main_dir,"inputs","raw_occurrence","bison_raw")))
  dir.create(file.path(main_dir,"inputs","raw_occurrence","bison_raw"),
  recursive=T)

# download BISON occurrence data for target taxa
#   there is also county distribution data
# if loop gets hung up, generally best to just try again instead of waiting
bison_raw <- data.frame()
us_cty_dist <- data.frame()

for(i in 1:length(taxon_names)){
  occ <- bison(species = taxon_names[i])
  bison_raw <- rbind.fill(bison_raw,occ$points)
  if(length(occ$counties>0)){
    occ$counties$taxon_name <- taxon_names[i]
    us_cty_dist <- rbind.fill(us_cty_dist,occ$counties)
  }
  print(taxon_names[i])
}
nrow(bison_raw) #2961
write.csv(bison_raw, file.path(main_dir,"inputs","raw_occurrence",
  "bison_raw","bison_R_download.csv"),row.names=FALSE)

### standardize column names

# keep only necessary columns
bison_raw <- bison_raw %>% dplyr::select(
  "name",
  "decimalLatitude","decimalLongitude",
  "basis","occurrenceID",
  "provider")

# rename columns
setnames(bison_raw,
  old = c("name",
          "basis","occurrenceID",
          "provider"),
  new = c("taxon_name",
          "basisOfRecord","nativeDatabaseID",
          "datasetName"),
  skip_absent=T)
bison_raw$database <- "BISON"
bison_raw$rightsHolder <- bison_raw$datasetName

# create species_name column
bison_raw$species_name <- NA
bison_raw$species_name <- sapply(bison_raw$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
sort(unique(bison_raw$species_name))

# recode standard columns
  # basis of record
sort(unique(bison_raw$basisOfRecord))
bison_raw <- bison_raw %>%
  mutate(basisOfRecord = recode(basisOfRecord,
    "Living" = "LIVING_SPECIMEN",
    "Observation" = "HUMAN_OBSERVATION",
    "Specimen" = "PRESERVED_SPECIMEN",
    "Fossil" = "FOSSIL_SPECIMEN",
    "Unknown" = "UNKNOWN"))

# check data
#percent.filled(bison_raw)
head(bison_raw)

# write file
write.csv(bison_raw, file.path(main_dir,"inputs","compiled_occurrence",
  "bison.csv"),row.names=FALSE)

rm(bison_raw)

# write file of county distribution
us_cty_dist <- left_join(us_cty_dist,taxon_list)
head(us_cty_dist)
write.csv(us_cty_dist, file.path(main_dir,"inputs","known_distribution",
  "BISON_US_county_distribution.csv"), row.names=FALSE)
