################################################################################

### 2-0_compile_exsitu_data.R
### Authors: Emily Beckman Bruns
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Smaller edits were added with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 13 December 2019
### Last full check and update: 30 November 2022
### R version 4.2.2

### DESCRIPTION:
  # This script takes a folder of CSV files representing accessions data from
  #   different institutions, combines them into one dataset, and standardizes
  #   some important fields.
  #

### INPUTS:
  #
  # 1. exsitu_standard_column_names (folder)
  #    Accessions data from ex situ collections survey; CSV files whose column
  #    names have been standardized by hand using instructions in the
  #    "Processing ex situ data" tab in Gap-analysis-workflow_metadata.xlsx
  #
  # 2. institution_data_table.csv
  #    Table with metadata for institutions who provided accessions data during
  #    the ex situ collections survey; see example in
  #    Gap-analysis-workflow_metadata.xlsx
  #
  # 3. target_taxa_with_synonyms.csv
  #    List of target taxa and synonyms; see example in
  #    Gap-analysis-workflow_metadata.xlsx
  #
  # 4. Accession-level data downloads from international crop genebank
  #    databases; download instructions are provided within the script below
  # (a) Genesys [Global Crop Diversity Trust]:
  #     https://www.genesys-pgr.org/a/overview
  # (b) WIEWS [FAO's World information and early warning system on plant genetic
  #     resources for food and agriculture]:
  #     https://www.fao.org/wiews/data/ex-situ-sdg-251/search/en/?no_cache=1
  #
  # 5. Polygons...

### OUTPUTS:
  # ExSitu_AllDataRaw_GenusFilterOnly_[[Sys.Date]].csv
  # ....

################################################################################
# Load libraries
################################################################################

#rm(list=ls())
my.packages <- c('plyr', 'tidyverse', 'textclean', 'spatialEco',
                 'maps', 'measurements', 'CoordinateCleaner', 'raster',
                 'data.table', 'terra', 'tidyterra')
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate
distinct <- dplyr::distinct

################################################################################
# Set working directory
################################################################################

# either set manually:
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/Global Tree Conservation Program/4. GTCP_Projects/Gap Analyses/ALL ACCESSION-LEVEL DATA FROM SURVEYS: 2017-2022"
# or use 0-1_set_workingdirectory.R script:
#source("./Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-1_set_working_directory.R")

# set folders where you have raw data (in) and want processed data to go (out)
data_in <- "exsitu_standard_column_names"
data_out <- "Quercus data compiled - Nov 2022"

################################################################################
# Load functions
################################################################################

# function to read in ex situ files from different folders/years and stack
read.exsitu.csv <- function(path,submission_year){
  # create list of paths to ex situ accessions CSV files in folder
  file_list <- list.files(path=path,pattern=".csv",full.names=TRUE)
  # read in each csv in path list to create list of dataframes
  file_dfs <- lapply(file_list,read.csv,header=TRUE,fileEncoding="LATIN1",
    strip.white=TRUE,colClasses="character",na.strings=c("","NA"))
  print(paste0("Number of files: ",length(file_dfs)))
    #sapply(file_dfs, nrow) # can look at number of rows in each csv
  for(file in seq_along(file_dfs)){
    df <- file_dfs[[file]]
    # add file name as column, to record home institution for each record
    df$filename <- rep(file_list[file],nrow(df))
    # remove file path portion
    df$filename <- mgsub(df$filename,c(paste0(path,"/"),".csv"),"")
    # add year of submission
    df$submission_year <- submission_year
    # remove extra blank columns
    t <- grepl("^X",names(df))
    if(length(unique(t))>1){
      #print(df$filename[1])
      df <- df[, -grep("^X", names(df))]
    }
    # replace strange characters in column names (arise from saving?)
    names(df) <- gsub(x = names(df), pattern = "ï\\.\\.", replacement = "")
    # add accession number if there isn't one
    if("acc_num" %in% names(df) & nrow(df[which(is.na(df$acc_num)),]) > 0){
      df[which(is.na(df$acc_num)),]$acc_num <- paste0("added",
        sprintf("%04d", 1:nrow(df[which(is.na(df$acc_num)),])))
      #print(nrow(df))
    } else if ("acc_no" %in% names(df) & nrow(df[which(is.na(df$acc_no)),]) > 0){
      df[which(is.na(df$acc_no)),]$acc_no <- paste0("added",
        sprintf("%04d", 1:nrow(df[which(is.na(df$acc_no)),])))
      #print(nrow(df))
    } else if (!("acc_num" %in% names(df)) & !("acc_no" %in% names(df))){
      df$acc_num <- paste0("added", sprintf("%04d", 1:nrow(df)))
      #print(nrow(df))
    } else {
      #print(paste("NO ACC NUM EDITS:",df$filename[1]))
    }
    # replace old df with new df
    file_dfs[[file]] <- df
    #print(head(file_dfs[[file]],n=2))
  }
  # stack all datasets using rbind.fill, which keeps non-matching columns
  #   and fills with NA; 'Reduce' iterates through and merges with previous
  # this may take a few minutes if you have lots of data
  all_data <- Reduce(rbind.fill, file_dfs)
    print(paste0("Number of rows: ",nrow(all_data)))
    print(paste0("Number of columns: ",ncol(all_data)))
  return(all_data)
}

################################################################################
################################################################################
# 1. Read in and stack all ex situ accessions data from survey of collections
################################################################################

### FIRST: After receiving accession data from institutions, you need to process
#   it manually. See here for instructions:
#   https://docs.google.com/spreadsheets/d/1p5HAS7vIE-3CbQcUmwrnuBGv5324dXy-42Iu6LlbX0E/edit?usp=sharing

# read in data from multiple surveys and stack, or just read in from one folder.
# this function also adds columns for 1) the file name [often equivalent to the
# "inst_short" institution nickname] 2) a sumbission year, 3) an accession
# number if one isn't given
### CHANGE BASED ON FOLDER(S) AND YEAR(S) YOU HAVE
raw_2022 <- read.exsitu.csv(file.path(main_dir,data_in,"2022"), "2022") #30
raw_2021 <- read.exsitu.csv(file.path(main_dir,data_in,"2021"), "2021") #9
raw_2020 <- read.exsitu.csv(file.path(main_dir,data_in,"2020"), "2020") #39
raw_2019 <- read.exsitu.csv(file.path(main_dir,data_in,"2019"), "2019") #143
raw_2018 <- read.exsitu.csv(file.path(main_dir,data_in,"2018"), "2018") #143
raw_2017 <- read.exsitu.csv(file.path(main_dir,data_in,"2017"), "2017") #143
  # for 2017 data, if no genus, fill with 'Quercus'
raw_2017[which(!is.na(raw_2017$species) & is.na(raw_2017$genus)),]$genus <- "Quercus"
# stack all data if you had multiple years
to_stack <- list(raw_2022,raw_2021,raw_2020,raw_2019,raw_2018,raw_2017)
all_data <- Reduce(rbind.fill, to_stack)

### IF APPLICABLE:
  # combine multiple inst_short columns
all_data <- tidyr::unite(all_data,"inst_short",
  c("inst_short","inst_short2"),
  sep=";",remove=T,na.rm=T)
  # combine multiple taxon name columns
  # we'll do the rest of the column standardizing later; this is to get the genus
all_data <- tidyr::unite(all_data,"taxon_full_name",
  c("taxon_full_name","full_sp_name","sp_full_name","taxon_name"),
  sep=";",remove=T,na.rm=T)

# fill in inst_short column with filename if none provided
all_data$inst_short[all_data$inst_short==""] <-
  all_data$filename[all_data$inst_short==""]
# remove rows with no inst_short
all_data <- all_data[which(all_data$inst_short!=""),]

# read in institution metadata file, for comparing list to data read in
inst_data <- read.csv(file.path(main_dir,
  "institution_data_table.csv"), stringsAsFactors = F)
inst_data$parent_R_file <- gsub("\\{PARENT\\}","",inst_data$parent_R_file)
## CHECK ALL INSTITUTIONS ARE PRESENT
#sort(unique(all_data$inst_short)) #333
#sort(unique(inst_data$inst_short)) #340
  ### INSTITUTIONS IN THE DATA BUT NOT IN METADATA TABLE
  # this should be "character(0)"
unique(all_data$inst_short)[!(unique(all_data$inst_short) %in% unique(inst_data$inst_short))]
  ### INSTITUTIONS IN METADATA TABLE BUT NOT IN THE DATA:
  # this should just be parent files (have mult institutions),
  # for example, PCNQuercus
unique(inst_data$inst_short)[!(unique(inst_data$inst_short) %in% unique(all_data$inst_short))]

# fill genus column so we can see which datasets have our target genus
  # fill genus column if not already filled
all_data <- all_data %>% separate("taxon_full_name","genus_temp",sep=" ",remove=F)
all_data[which(is.na(all_data$genus)),]$genus <- all_data[which(is.na(all_data$genus)),]$genus_temp
  # standardize capitalization
all_data$genus <- str_to_title(all_data$genus)
  ### FIX GENUS MISSPELLINGS OR ABBREVIATIONS
sort(unique(all_data$genus))
all_data$genus <- mgsub(all_data$genus,
  c("^Q$","^Q\\.$","Querces","Quercues","Querucs","Querus","Qurercus","Cyclobalanopsis"),
  "Quercus",fixed=F)
all_data$taxon_full_name <- mgsub(all_data$taxon_full_name,
  c("Q\\.","Cyclobalanopsis"), "Quercus",fixed=F)

# remove duplicate data - from previous years or networks
### UPDATE THIS SECTION BASED ON YOUR SPECIFIC NEEDS; THIS IS FOR QUERCUS 2022:
nrow(all_data) #259031
all_data$file_inst_year <- paste(all_data$filename,all_data$inst_short,all_data$submission_year)
  # 1) remove datasets that have no Quercus data (so we know they aren't dups)
inst_genera <- all_data %>%
  group_by(filename,inst_short,submission_year) %>%
  mutate(genera = paste(unique(genus),collapse = ", ")) %>%
  ungroup() %>% select(inst_short,filename,submission_year,genera)
yes_Q <- unique(inst_genera[which(grepl("Quercus",inst_genera$genera)),]) %>%
  mutate(file_inst_year = paste(filename,inst_short,submission_year))
all_data <- all_data %>% filter(file_inst_year %in% yes_Q$file_inst_year)
nrow(all_data) #156714
  # 2) remove old datasets if there is newer data from 2018-2020
all_data$inst_and_year <- paste(all_data$inst_short,all_data$submission_year)
find_dups <- all_data %>%
  filter(submission_year != "2021" & submission_year != "2022") %>%
  group_by(inst_short,submission_year) %>%
  count() %>% ungroup() %>%
  arrange(desc(submission_year)); find_dups
remove_dups <- find_dups[which(duplicated(find_dups$inst_short)),]
remove_dups <- paste(remove_dups$inst_short,remove_dups$submission_year); remove_dups
all_data <- all_data %>% filter(!(inst_and_year %in% remove_dups))
nrow(all_data) #127495
  # 3) if there is duplicate network data, keep only the data direct from garden
network_data <- all_data %>% filter(inst_short != filename &
  (submission_year != "2022" & submission_year != "2021"))
own_data <- all_data %>% filter(inst_short == filename &
  (submission_year != "2022" & submission_year != "2021"))
remove <- sort(unique(network_data[which(
  network_data$inst_short %in% own_data$inst_short),]$file_inst_year)); remove
all_data <- all_data %>%
  filter(!(file_inst_year %in% remove))
nrow(all_data) #126189

# check out column names
  #not_all_na <- function(x) any(!is.na(x))
  #all_data <- all_data %>% select(where(not_all_na))
sort(colnames(all_data))
## IF NEEDED: separate column into multiple
#all_data <- all_data %>% separate("specific",
#  c("infra_rank_add","infra_name_add"),sep=" ",remove=T,fill="right")
## IF NEEDED: see which datasets have extraneous columns so you can fix manually
##  in the raw data as desired; update line below
#unique(all_data$filename[all_data$locality.1 !=""])
## IF NEEDED: merge similar columns (you may not need to do this if no schema
##  mistakes were made when manually editing column names)
all_data <- tidyr::unite(all_data,"infra_name",
  c("infra_name","intra_name","specific_name","specific","Infraspecific.Epithet"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"infra_rank",
  c("infra_rank","intra_rank","specific_rank","specific2","Infraspecific.Rank"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"cultivar",
  c("cultivar","Cultivar.Epithet"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"hybrid",
  c("hybrid","Hybrid","hybrid.1"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"acc_num",
  c("acc_num","acc_no","acc_num2"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"coll_year",
  c("coll_year","acq_year","aqu_yr","aqu_year","coll_yr","year","acc_date","Accession.Date"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"notes",
  c("notes","orig_notes","availability","note","Notes","notes_1","notes_2","notes_3","notes_4","num_individuals_sampled","num_samples_per_individual"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"municipality",
  c("municipality","city",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"coll_name",
  c("coll_name","coll_name2",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"coll_num",
  c("coll_num","coll_no",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"uncertainty",
  c("coord_precision","coord_prec",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"country",
  c("country","country_code",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"county",
  c("county","county.1",),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"orig_source",
  c("orig_source","donor","donor_name","land_owner","orig_source.1","source","source2"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"garden_loc",
  c("garden_loc","garden_loc_campus","loc"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"habitat",
  c("habitat","habitat2"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"taxon_verif",
  c("id_by","id_method","id_notes","id_year","name_determ","taxon_det","Taxonomic.Verification"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"orig_lat",
  c("orig_lat","lat"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"lin_num",
  c("lin_num","lin_no","source_acc_num"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"locality",
  c("locality","loaclity","locality_details","locality_notes","locality.1","locality2","locallity","other_region","site_notes"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"orig_long",
  c("orig_long","long"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"state",
  c("state","maj_region"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"county",
  c("county","min_region"),sep=";",remove=T,na.rm=T)
all_data$no_alive[which(!is.na(all_data$no_alive))] <- paste0(all_data$no_alive[which(!is.na(all_data$no_alive))]," alive")
all_data$num_rec[which(!is.na(all_data$num_alive))] <- paste0(all_data$num_alive[which(!is.na(all_data$num_alive))]," alive")
all_data <- tidyr::unite(all_data,"num_indiv",
  c("no_alive","num_alive","num_indiv","no_plants","num_plants"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"condition",
  c("condition","plant_status","status"),sep=";",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"rec_as",
  c("rec_as","rec_material"),sep=";",remove=T,na.rm=T)
## IF NEEDED: rename columns
all_data <- all_data %>% rename(gps_det = coord_det)
## IF NEEDED: remove rows for dead/removed individuals
sort(unique(all_data$condition))
all_data <- all_data %>%
  filter(all_data$condition != "Dead, Flagged" &
         all_data$condition != "Dead" &
         all_data$condition != "Removed" &
         all_data$condition != "Lost" &
         all_data$condition != "Rmvd")
nrow(all_data) #125865
sort(colnames(all_data))

## REMOVE UNUSED COLUMNS (SELECT WHAT YOU WANT TO KEEP)
keep_col <- c("acc_num","assoc_sp","author","coll_name","coll_num",
  "coll_year","country","county","cultivar","dataset_year","filename",
  "garden_loc","genus","germ_type","gps_det","hybrid","infra_name",
  "infra_rank","inst_short","lin_num","locality","municipality","notes",
  "num_indiv","orig_lat","orig_long","orig_source","prov_type","rec_as",
  "species","state","submission_year","taxon_full_name","taxon_verif",
  "uncertainty")
all_data <- all_data[,keep_col]

# remove leading, trailing, and middle (e.g., double space) whitespace,
#   to prevent future errors
all_data <- as.data.frame(lapply(all_data, function(x) str_squish(x)),
  stringsAsFactors=F)
# replace "" cells with NA in whole dataset
all_data[all_data == ""] <- NA
# add data source colunm
all_data$data_source <- "ex_situ_BG_survey"
# first fix some common lat/long character issues before replacing non-ascii
all_data$orig_lat <- mgsub(all_data$orig_lat,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
all_data$orig_long <- mgsub(all_data$orig_long,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
# now replace all non-ascii characters
all_data <- as.data.frame(lapply(all_data,replace_non_ascii),stringsAsFactors=F)

################################################################################
# 2. Compile Genesys data (genebanks)
################################################################################

### FIRST, download data from Genesys [Global Crop Diversity Trust]:
# Go to https://www.genesys-pgr.org/a/overview
# Search for target genus/genera (in the TAXONOMY tab on the left) then click
# "APPLY FILTERS" at the top left.
# Click the "ACCESSIONS" tab in the top bar, then click "DOWNLOAD ZIP" on the
# right.
# Move the downloaded folder to your main working directory (main_dir) and
# rename to "genesys-accessions"
# You can view field metadata here: https://www.genesys-pgr.org/documentation/basics

# load in data
  # collection site description
genPath <- paste0(main_dir,"/genesys-accessions/coll.csv")
  # I think the "Found and resolved improper quoting" warning is ok
gen_col <- data.table::fread(file = genPath, header = TRUE)
  str(gen_col)
  nrow(gen_col) ; length(unique(gen_col$genesysId)) #333
  # metadata about record
genPath <- paste0(main_dir,"/genesys-accessions/core.csv")
gen_core <- data.table::fread(genPath, header = TRUE)
  str(gen_core)
  nrow(gen_core) ; length(unique(gen_core$genesysId)) #428
  # spatial info
  #   for some reason the headers dont read in correctly with fread;
  #   will use read.csv just to get the headers and add them to the fread df
genPath <- paste0(main_dir,"/genesys-accessions/geo.csv")
gen_geo <- data.table::fread(genPath, header = TRUE)
  str(gen_geo)
  colnames(gen_geo)
gen_geo_head <- read.csv(genPath)
  colnames(gen_geo_head)
gen_geo <- gen_geo %>% select(-V8) # remove the last column (nothing in it)
colnames(gen_geo) <- colnames(gen_geo_head)
  str(gen_geo)
  nrow(gen_geo) ; length(unique(gen_geo$genesysId)) #427

# combine all dataframes
genesys <- Reduce(full_join,list(gen_col,gen_core,gen_geo)) #,gen_names
str(genesys); nrow(genesys)
# create taxon name column
## We don't need to do this here since we work with all taxon names below!
  # fix up infra rank/name column a little
#sort(unique(genesys$subtaxa))
#genesys$subtaxa <- mgsub(genesys$subtaxa,
#  c("\\(l.\\) dun.","^s$","~sp","^a$","^i$","^m$"),"",fixed=F)
#genesys$subtaxa <- mgsub(genesys$subtaxa,c("var.","subsp."),c("var. ","subsp. "))
#genesys$subtaxa <- str_squish(genesys$subtaxa)
#sort(unique(genesys$subtaxa))
  # create full name
#genesys$taxon_full_name <- paste(genesys$genus,genesys$species,genesys$subtaxa)
#genesys$taxon_name <- str_squish(genesys$taxon_name)
  # replace hybrid x with +
#genesys$taxon_name <- gsub(" x "," +",genesys$taxon_name)
#sort(unique(genesys$taxon_name))
  # add 'Genesys' to ID
genesys$UID <- paste0("Genesys-",genesys$genesysId)
  # paste some similar columns together
genesys <- tidyr::unite(genesys,"coll_year",
  c("collDate","acqDate"),sep=";",remove=T,na.rm=T)

# rename columns to match ex situ data and select only those we need
genesys_sel <- genesys %>%
  rename(taxon_full_name = fullTaxa,
         coll_num = collNumb,
         locality = collSite,
         notes = uuid,
         inst_short = instCode,
         acc_num = acceNumb,
         infra_rank = subtaxa,
         country = origCty,
         orig_lat = latitude,
         orig_long = longitude,
         gps_det = method,
         germ_type = storage,
         orig_source = collSrc,
         prov_type = sampStat) %>%
  select(UID,taxon_full_name,coll_num,coll_year,locality,notes,
    inst_short,acc_num,infra_rank,country,orig_lat,orig_long,gps_det,
    germ_type,orig_source,prov_type,genus,species,uncertainty)

# filter by target taxa & add taxon data
## We do this later!
#genesys_sel <- left_join(genesys_sel,taxon_list)
#genesys_target <- genesys_sel %>% filter(!is.na(taxon_name_acc))
#nrow(genesys_target)

# add data source colun
genesys_sel$data_source <- "Genesys"

# join to exsitu data
all_data <- rbind.fill(all_data,genesys_sel)
nrow(all_data) #126293

## YOU EITHER NEED TO REMOVE GRIN OR REMOVE US INSTITUTIONS IN GENSYS:
## To remove US gensys:
# US institutions covered by GRIN (I think); we will remove these
# pulled this list from CWR Gap Analysis 2020 (Khoury et al.)...
#USDAcodes <- c("USA003" ,"USA004", "USA005" ,"USA016" ,"USA020",
#"USA022", "USA026", "USA028", "USA029", "USA042" ,"USA047", "USA049",
#"USA074", "USA108", "USA133", "USA148", "USA151", "USA167", "USA176",
#"USA390", "USA955", "USA956", "USA970", "USA971", "USA995")
#all_data <- all_data %>% filter(!(inst_short %in% USDAcodes))
## To remove GRIN:
all_data <- all_data %>% filter(inst_short != "GRIN_NPGS")

nrow(all_data) #126267

################################################################################
# 3. Compile WIEWS data (genebanks)
################################################################################

# 23 Nov 2022: For some reason the data aren't downloading; skipping this
# section for Quercus (Genesys gets most of the genebanks anyways).

### FIRST, download data from WIEWS [FAO's World information and early warning
# system on plant genetic resources for food and agriculture]:
# Go to https://www.fao.org/wiews/data/ex-situ-sdg-251/search/en/?no_cache=1
# In the "Crop Wild Relatives" dropdown, select "Included".
# Type your target genus into the "Genus" box and press the "+" button on the
# right. Add additional target genera like this as needed.
# Scroll down and press the magnifying glass icon on the bottom right.
# Now click the "Download Results" button on the top right.
# Move the downloaded file to your main working directory (main_dir) and
# rename to "wiews-exsitu.csv"

# read in data
wiewsPath <- paste0(main_dir,"/wiews_exsitu.csv")
#wiews <- data.table::fread(wiewsPath, header = TRUE)
wiews <- read.csv(wiewsPath)
str(wiews) #

# filter out Genesys data
wiews <- wiews %>%
  filter(Source.of.information != "Genesys (https://www.genesys-pgr.org)")
nrow(wiews) #

# create taxon name column
## We don't need to do this here since we work with all taxon names below!
  # Spilt name to get at genus and species
#wiews$name <- wiews$Taxon
#wiews <- tidyr::separate(data = wiews, "name",
#  into =c('genus','spec','sub1','sub2','sub3', 'sub4'),sep=' ')
#  # from "wiewsTransform.R" by Dan Carver:
#  #   Function to split full name into taxon/species
#setSpecies <- function(dataFrame){
#  if(!is.na(dataFrame$sub1)){
#    dataFrame$species <- paste(dataFrame$spec,dataFrame$sub1,
#                               dataFrame$sub2, sep="_")
#  }
#  if(is.na(dataFrame$sub1)){
#    dataFrame$species <- dataFrame$spec
#  }
#  return(dataFrame)
#}
#wiews2 <- setSpecies(wiews)
# remove all NA
#df6 <- str_remove_all(df4$species, 'NA') %>%
#  str_remove_all("__")
#df4$species <- df6
#df4 <- subset(x = df4, select = -c(spec,sub1,sub2) )

# rename columns to match ex situ data and select only those we need
wiews_sel <- wiews %>%
  rename(inst_short = Holding.institute.code,
         acc_num = Accession.number,
         taxon_full_name = Taxon,
         genus = Genus,
         species = Species,
         rec_date = Acquisition.date..YYYY.MM.,
         country = Country.of.origin..ISO3.,
         prov_type = Biological.status,
         inst_short2 = Genebank.s..holding.safety.duplications...code,
         orig_lat = Latitude.of.collecting.site..decimal.degrees.format.,
         orig_long = Longitude.of.collecting.site..decimal.degrees.format.,
         orig_source = Collecting.acquisition.source,
         germ_type = Type.of.germplasm.storage) %>%
  select(inst_short,acc_num,taxon_full_name,genus,species,rec_date,country,
         prov_type,inst_short2,orig_lat,orig_long,orig_source,germ_type)

# !!! JUST FOR NOW; need to check exact to filter out !!!
#   US institutions (covered by GRIN?)
# looks like they are all removed by removing Genesys anyways, so we'll skip
#wiews_target <- wiews_sel %>% filter(!grepl("USA",inst_short))
nrow(wiews_sel) #

# add data source colun
wiews_sel$data_source <- "FAO-WIEWS"

# join to exsitu data
all_data <- rbind.fill(all_data,wiews_sel)
nrow(all_data) #

################################################################################
# 4. Save raw output of compiled data for target genera
#     (can be used to look for hybrids/cultivars, which are removed in next step)
################################################################################

# create folder for output data
if(!dir.exists(file.path(main_dir,data_out)))
  dir.create(file.path(main_dir,data_out), recursive=T)

all_data2 <- all_data
nrow(all_data2) #126267

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "target_taxa_with_synonyms.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
str(taxon_list) #739
unique(taxon_list$taxon_region)
  # can select just a specific region if you want
#taxon_list <- taxon_list %>%
#  filter(taxon_region == "Mesoamerican" |
#         taxon_region == "Mesoamerican & US" |
#         taxon_region == "US")
nrow(taxon_list) #728

# preserve original taxon name
all_data2$taxon_full_name_orig <- all_data2$taxon_full_name

# remove rows not in target genus/genera
target_genera <- c("Quercus")
all_data2 <- all_data2 %>% filter(genus %in% target_genera)
nrow(all_data); nrow(all_data2) #126267 ; 52907

### CHECK OUT THE HYBRID COLUMN ###
sort(unique(all_data2$hybrid))
# standardize so only one hybrid symbol is used ( x )
all_data2$hybrid <- mgsub(all_data2$hybrid,
  c(" A ","^A ","^A$"," X ","^X"," _ ","^_ ","^_$","^1$","\\*","^H$","^hyb$","Hybrid"),
  " x ", fixed=F)
# add "x" at beginning then remove duplicates
all_data2$hybrid[!is.na(all_data2$hybrid)] <- paste0("x ",all_data2$hybrid[!is.na(all_data2$hybrid)])
all_data2$hybrid <- str_squish(all_data2$hybrid)
all_data2$hybrid <- str_to_lower(all_data2$hybrid)
all_data2$hybrid <- mgsub(all_data2$hybrid,c("x x","x h$"),"x",fixed=F)
sort(unique(all_data2$hybrid))

# create concatenated taxon_full_name column
all_data2 <- tidyr::unite(all_data2, "taxon_full_name_concat",
  c(genus,hybrid,species,infra_rank,infra_name,cultivar), sep=" ", remove=F,
  na.rm=T)

# when blank, fill taxon_full_name column with concatenated full name
all_data2[is.na(all_data2$taxon_full_name),]$taxon_full_name <-
  all_data2[is.na(all_data2$taxon_full_name),]$taxon_full_name_concat
#unique(all_data2$taxon_full_name)

# standardize common hybrid signifiers in taxon_full_name
all_data2$taxon_full_name <- mgsub(all_data2$taxon_full_name,
  c(" A "," A","_"," X","\\*","×"," \\? ")," x ",fixed=F)
# make sure the author is separated in taxon_full_name
all_data2$taxon_full_name <- gsub("\\("," (",all_data2$taxon_full_name)
all_data2$taxon_full_name <- str_squish(all_data2$taxon_full_name)
#unique(all_data2$taxon_full_name)

# add space after periods in taxon_full_name
all_data2$taxon_full_name <- gsub(".",". ",all_data2$taxon_full_name,fixed=T)
all_data2$taxon_full_name <- str_squish(all_data2$taxon_full_name)

# summary of genera for each institution
gen_summary <- all_data2 %>%
  arrange(genus) %>%
  rename(genera = genus) %>%
  group_by(inst_short) %>%
  mutate(
    genera = paste(unique(genera), collapse = '; ')) %>%
  ungroup() %>%
  distinct(inst_short,genera)
head(as.data.frame(gen_summary))
nrow(gen_summary) # number of institutions

################################################################################
# 5. Further standardize taxon name, then keep data for target taxa only
################################################################################

# fixing some taxon name issues I noticed:
all_data2$taxon_full_name <- mgsub(all_data2$taxon_full_name,
  c("Quercus Quercus"),
  c("Quercus"))
#all_data2$taxon_full_name <- mgsub(all_data2$taxon_full_name,
#  c("Corylus Acolurnoides","Prunus Ayedoensis","Pyrus (Malus) fusca"),
#  c("Corylus x colurnoides","Prunus x yedoensis","Malus fusca"))

## first change hybrid notation temporarily: remove space so it stays together
all_data2$taxon_full_name <- gsub(" x "," +",all_data2$taxon_full_name)
# separate out taxon full name and trim whitespace again
    # this warning is ok: "Expected 9 pieces. Additional pieces discarded..."
all_data2 <- all_data2 %>% separate("taxon_full_name",
  c("genus_new","species_new","extra1","extra2",
    "extra3","extra4","extra5","extra6","extra7"),sep=" ",extra="warn",
    remove=F,fill="right")
all_data2 <- as.data.frame(lapply(all_data2,str_squish),stringsAsFactors=F)
# replace genus_new with genus, since we fixed that up in the previous section
all_data2$genus_new <- all_data2$genus

## REMOVE RECORDS WITHOUT SPECIES NAME

# remove records with no/non-standard specific epithet (mostly cultivars or
#   'sp.' or questionable '?') by looking in species name column
all_data3 <- all_data2 %>%
  filter(!grepl("\"",species_new) &
                !grepl("\'",species_new) &
                !grepl("\\[",species_new) &
                !grepl("\\(",species_new) &
                !grepl("\\.",species_new) &
                !grepl("[A-Z]",species_new) &
                !grepl("[0-9]",species_new) &
                !grepl("\\?",species_new) &
                !is.na(species_new))
nrow(all_data3) #51587
# see records removed; can add anything you want to fix to the
#   "fixing some taxon name issues I noticed" section above:
#sort(unique(anti_join(all_data2,all_data3)$taxon_full_name))

## FIND INFRATAXA

## look for infrataxa key words
# make data in all "extra" columns lower case
sp_col <- grep("^species_new$", colnames(all_data3))
all_data3[,sp_col:(sp_col+5)] <- as.data.frame(sapply(
  all_data3[,sp_col:(sp_col+5)], tolower), stringsAsFactors=F)
# create matrix of all "extra" species name columns, to search for
#   infraspecific key words
search.col <- matrix(cbind(all_data3$extra1,all_data3$extra2,all_data3$extra3,
  all_data3$extra4,all_data3$extra5,all_data3$extra6,all_data3$extra7),
  nrow=nrow(all_data3))
#str(search.col)
# search the "extra" column matrix for matches to infraspecific key words
matches_i <- which(search.col=="variety"|search.col=="var"|search.col=="var."|
                  search.col=="v"|search.col=="v."|search.col=="va"|
                 search.col=="subspecies"|search.col=="subsp"|
                  search.col=="subsp."|search.col=="ssp"|search.col=="ssp."|
                  search.col=="subs."|search.col=="spp."|search.col=="sub."|
                 search.col=="infra"|
                 search.col=="forma"|search.col=="form"|search.col=="fma"|
                  search.col=="fo"|search.col=="fo."|search.col=="f"|
                  search.col=="f.",arr.ind=T)
matches_i[,2] <- matches_i[,2]+sp_col
# create new infra_rank column and fill with "extra" contents that matched
#   infraspecific key words
all_data3$infra_rank_new <- NA
all_data3$infra_rank_new[matches_i] <- all_data3[matches_i]
#unique(all_data3$infra_rank_new) # check results

# create new infra_name column and fill with next column over from "extra"
#   contents that matched infraspecific key word
all_data3$infra_name_new <- NA
matches_i[,2] <- matches_i[,2]+1
all_data3$infra_name_new[matches_i] <- all_data3[matches_i]
#sort(unique(all_data3$infra_name_new))

# standardize infraspecific rank names
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
  grep("^v$|^v.$|^var$|^variety$|^va$",all_data3$infra_rank_new), "var.")
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
  grep("^subspecies$|^subsp$|^ssp$|^ssp.$|^subs.$|^spp.$|^sub.$",
  all_data3$infra_rank_new), "subsp.")
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
 grep("^forma$|^form$|^fma$|^fo$|^fo.$|^f$",all_data3$infra_rank_new), "f.")
unique(all_data3$infra_rank_new)

## CREATE FINAL TAXON FULL NAME FOR FILTERING

# create new taxon full name column
all_data3$taxon_full_name <- NA
  # select rows with infraspecific name and concatenate
yes_infra <- which(!is.na(all_data3$infra_rank_new) &
  !is.na(all_data3$infra_name_new))
all_data3$taxon_full_name[yes_infra] <- paste(all_data3$genus_new[yes_infra],
  all_data3$species_new[yes_infra], all_data3$infra_rank_new[yes_infra],
  all_data3$infra_name_new[yes_infra],sep=" ")
  # select rows without infraspecific name and concatenate
all_data3$taxon_full_name[-yes_infra] <- paste(all_data3$genus_new[-yes_infra],
  all_data3$species_new[-yes_infra],sep=" ")
# switch hybrid symbol back to " x "
all_data3$taxon_full_name <- gsub(" \\+"," x ",all_data3$taxon_full_name)
#sort(unique(all_data3$taxon_full_name))

## FILTER OUT NON-TARGET TAXA

# rename some taxon name columns to preserve originals
all_data3 <- all_data3 %>%
  rename(taxon_name = taxon_full_name,
                genus_orig = genus,
                species_orig = species,
                infra_rank_orig = infra_rank,
                infra_name_orig = infra_name)
all_data3 <- all_data3 %>%
  rename(genus = genus_new,
                species = species_new,
                infra_rank = infra_rank_new,
                infra_name = infra_name_new)

# join dataset to taxa list
  # add genus_species column to taxon list
taxon_list$genus_species <- NA
taxon_list$genus_species <- sapply(taxon_list$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
  # join by taxon name
all_data4 <- left_join(all_data3,taxon_list)
  # if no taxon match, join again just by species name
need_match <- all_data4[which(is.na(all_data4$taxon_name_status)),]
  nrow(need_match) #4149
    # remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data4)-ncol(taxon_list)+1)]
    # remove taxon_name col from taxon data so it doesn't match
taxon_list_sp <- taxon_list %>% select(-taxon_name)
  # create genus_species column
need_match$genus_species <- paste(need_match$genus,need_match$species)
    # new join by genus_species
need_match <- left_join(need_match,taxon_list_sp)
  # bind together new matches and previously matched
matched <- all_data4[which(!is.na(all_data4$taxon_name_status)),]
all_data4 <- rbind(matched,need_match)
  table(all_data4$taxon_name_status) # Accepted: 47975 | Synonym: 1887
  # see how many rows have taxon name match
nrow(all_data4[which(!is.na(all_data4$taxon_name_status)),]) #49862 (was 13762 before)
table(all_data4$taxon_region)
head(all_data4)

### CHECK UNMATCHED SPECIES, TO ADD TO SYNONYM LIST AS NECESSARY ###
check <- all_data4 %>% filter(is.na(taxon_name_status))
check <- data.frame(taxon_name = sort(unique(check$taxon_name)))
head(check); nrow(check) #290
# write file for checking, as desired
  # IF YOU FIND MISSPELLINGS AND/OR ADDITIONAL SYNONYMS, YOU CAN ADD THEM TO
  #   YOUR TARGET TAXA LIST AND GO BACK TO THE LINE WHERE WE "read in target
  #   taxa list" (SECTION 4) AND RUN AGAIN FROM THERE
write.csv(check, file.path(main_dir,data_out,"ExSitu_UnmatchedSpecies.csv"),
  row.names = F)

# keep only matched names
all_data5 <- all_data4 %>% filter(!is.na(taxon_name_status))
nrow(all_data5) #49862

# final part for removing duplicate data - from previous years or networks
### UPDATE THIS SECTION BASED ON YOUR SPECIFIC NEEDS; THIS IS FOR QUERCUS 2022:
  # for 2021 and 2022 data, only remove older data for species provided
all_data6 <- all_data5
dup_2022 <- unique(all_data6[which(all_data6$submission_year=="2022"),
  c("taxon_name","inst_short")]); dup_2022
for(i in 1:nrow(dup_2022)){
  all_data6 <- all_data6 %>%
   filter(!(taxon_name == dup_2022[i,1] &
            inst_short == dup_2022[i,2] &
            submission_year != "2022"))
}; nrow(all_data6) #46666
dup_2021 <- unique(all_data6[which(all_data6$submission_year=="2021"),
  c("taxon_name","inst_short")]); dup_2021
for(i in 1:nrow(dup_2021)){
  all_data6 <- all_data6 %>%
  filter(!(taxon_name == dup_2021[i,1] &
           inst_short == dup_2021[i,2] &
           submission_year != "2021"))
}; nrow(all_data6) #46623 (was 10831 before)

################################################################################
# 6. Standardize important columns
################################################################################

## KEEP ONLY NECESSARY COLUMNS
keep_col <- c(
  # key data
  "UID","inst_short","taxon_name_accepted","acc_num","prov_type","num_indiv",
  # locality
  "orig_lat","orig_long","locality","municipality","county","state","country",
  "gps_det","uncertainty","assoc_sp",
  # source
  "orig_source","lin_num","coll_name","coll_num","coll_year",
  # material info
  "germ_type","garden_loc","rec_as",
  # other metadata
  "notes","filename","submission_year","dataset_year","data_source",
  # taxon name details
  "taxon_full_name_orig","taxon_full_name_concat",
  "genus","hybrid","species","infra_rank","infra_name","cultivar",
  "taxon_name_status","author","taxon_verif",
  # OPTIONAL additional taxon metadata
  "taxon_region","iucnredlist_category"
)
all_data7 <- all_data6[,keep_col]

# add institution metadata to data
str(inst_data)
  # SELECT COLUMNS YOU WANT ADDED TO ALL DATA:
inst_data <- inst_data %>%
  select(inst_short,inst_country,inst_lat,inst_long) %>%
  arrange(inst_lat) %>%
  distinct(inst_short,.keep_all=T) %>%
  arrange(inst_short); inst_data
all_data7 <- left_join(all_data7,inst_data)
str(all_data7)

##
## Provenance type
##

# look at column contents and change below phrases as needed
all_data7$orig_prov_type <- all_data7$prov_type
all_data7$prov_type <- str_to_lower(all_data7$prov_type)
sort(unique(all_data7$prov_type))
## IF NEEDED: transfer contents of one column to another column, if data
#   needs to be preserved but is in the wrong place
# Warnings are fine here
all_data7[which(all_data7$prov_type=="130"|grepl("^130) ",all_data7$prov_type)),]$notes <- "Semi-natural/sown"
all_data7[which(all_data7$prov_type=="410"|grepl("^410) ",all_data7$prov_type)),]$notes <- "Breeding/research material: Breeder's line"
all_data7[which(all_data7$prov_type=="300"|grepl("^300) ",all_data7$prov_type)),]$notes <- "Traditional cultivar/landrace"
all_data7[which(all_data7$prov_type=="400"|grepl("^400) ",all_data7$prov_type)),]$notes <- "Breeding/research material"
all_data7[which(all_data7$prov_type=="500"|grepl("^500) ",all_data7$prov_type)),]$notes <- "Advanced or improved cultivar (conventional breeding methods)"
all_data7[which(grepl("county",all_data7$prov_type)),]$notes <- all_data7[which(grepl("county",all_data7$prov_type)),]$prov_type

# standardize column by searching for keywords and replacing with standard value
  # remove confusing words/phrases
all_data7$prov_type <- mgsub(all_data7$prov_type,
  c("not of known wild origin","could be cultivated","not of wild source",
    "wildsourceunsure","cultivated plant of known "), "")
  # ex wild (Z)
all_data7$prov_type <- ifelse(grepl(paste(
  c("indirect","ex wild","^z$","cultivated from wild","c ex w","cw",
    "g from w plant","g ex w","cultivatedpropagatedfromwildsource",
    "(indirect) wild origin"),
  collapse = "|"), all_data7$prov_type),"Z",all_data7$prov_type)
  # wild (W)
all_data7$prov_type <- ifelse(grepl(paste(
  c("wild","wld","collect","^w$","^w\\*","^\\(w\\)$","wd","w\\?","genetic",
    "100","110","130"),
  collapse = "|"), all_data7$prov_type),"W",all_data7$prov_type)
  # native to site (N)
all_data7$prov_type <- ifelse(grepl(paste(
  c("original to site","spontaneous","^n$","existing on site","native"),
  collapse = "|"), all_data7$prov_type),"N",all_data7$prov_type)
  # unknown (U)
all_data7$prov_type <- ifelse(grepl(paste(
  c("^\\(u\\)$","^u$","^u ","unsure","insufficient data","unknown","\\?","un",
    "need to confirm","breeding","400","410","donated","999","purchased"),
  collapse = "|"), all_data7$prov_type),"U",all_data7$prov_type)
  # cultivated (H)
all_data7$prov_type <- ifelse(grepl(paste(
  c("cultiva","garden","^c$","^g$","^g ","^h$","horticult","landrace","clone",
    "300","500","^g\\."),
  collapse = "|"), all_data7$prov_type),"H",all_data7$prov_type)
# check one last time
sort(unique(all_data7$prov_type))
  # not given (NG) ; everything else
all_data7$prov_type <- ifelse(all_data7$prov_type!= "W" &
  all_data7$prov_type != "Z" & all_data7$prov_type != "H" &
  all_data7$prov_type != "N" & all_data7$prov_type != "U",
  "NG",all_data7$prov_type)
all_data7$prov_type[which(is.na(all_data7$prov_type))] <- "NG"

# check results
table(all_data7$prov_type)
#    H   N    NG     U     W    Z
# 10628  7 13541  9765 12055  627

##
## B) Number of Individuals
##

all_data7$num_indiv <- str_to_lower(all_data7$num_indiv)
sort(unique(all_data7$num_indiv))
  ## IF NEEDED: replace unwanted characters
  all_data7$num_indiv <- mgsub(all_data7$num_indiv,
    c(" \\(all dead\\)"," \\(removed 2017\\)","\\?",","," plants"," or 5",
      " pieces","ca\\. ","ca\\."," alive","\\*","1 \\(","\\)","\\+",
      " in terra","alive;","a a a a a a a a a a a",";[0-9][0-9]",";1",
      ";2",";3",";4",";5",";7",";0",";"),
      c(""), fixed=F)
  sort(unique(all_data7$num_indiv))
  # change type to numeric and replace NA with 1
  # "NAs introduced by coercion" warning ok
all_data7$num_indiv <- as.numeric(all_data7$num_indiv)
all_data7$num_indiv[which(is.na(all_data7$num_indiv))] <- 1

# check results
sort(unique(all_data7$num_indiv))
nrow(all_data7) #46623

# remove records with no individuals; save as separate file
no_indiv <- all_data7[which(all_data7$num_indiv == 0),]
nrow(no_indiv) #746
write.csv(no_indiv, file.path(main_dir,data_out,
  paste0("ExSitu_Dead_", Sys.Date(), ".csv")),row.names = F)
  # save to in situ data folder as well
  #if(!dir.exists(file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ")))
  #  dir.create(file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ"),
  #  recursive=T)
  #write.csv(no_indiv, file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ",
  #  paste0("ExSitu_Dead_", Sys.Date(), ".csv")),row.names = F)
all_data7 <- all_data7[which(all_data7$num_indiv > 0),]
nrow(all_data7) #45877

##
## c) Combine individuals (same institution and accession number)
##      so that everything is (hopefully) at the accession-level
##

# preserve original acc_num before removing indiviual signifiers
all_data7$orig_acc_num <- all_data7$acc_num

# combine duplicates (same acc num)
all_data7 <- all_data7 %>%
  group_by(inst_short,acc_num,taxon_name_accepted) %>%
  mutate(num_indiv = sum(as.numeric(num_indiv)),
                germ_type = paste(unique(germ_type),collapse="; "),
                garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  ungroup() %>%
  distinct(inst_short,acc_num,taxon_name_accepted,.keep_all=T)
nrow(all_data7) #43284

# fill spaces in acc_num with dashes
all_data7$acc_num <- gsub(" ","-",all_data7$acc_num)

# now we will try to find patterns that identify accession numbers
# with individual identifiers - we want everything to be at the
# acccession level, not at the individual level for this analysis

# can look at what will be removed in the acc_num before we
# actually combine things; these patterns seem to work for all
all_data7[which(grepl("\\*",all_data7$acc_num)),]$acc_num
all_data7[which(grepl("/[1-9]$",all_data7$acc_num)),]$acc_num
  # these don't seem to work well with this dataset but may work for others:
#all_data7[which(grepl("_",all_data7$acc_num)),]$acc_num

# remove individual-specific identifiers (to combine dup accessions)
# "Expected 1 pieces" warning here is fine.
all_data7 <- all_data7 %>%
  separate("acc_num","acc_num",
    sep="\\*|/[1-9]$",remove=F) %>%
  group_by(inst_short,acc_num,taxon_name_accepted) %>%
  mutate(num_indiv = sum(as.numeric(num_indiv)),
                germ_type = paste(unique(germ_type),collapse="; "),
                garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  ungroup() %>%
  distinct(inst_short,acc_num,taxon_name_accepted,.keep_all=T)
nrow(all_data7) #39927

# create subset of records with acc_num longer than 9 characters
#   (these are usually the ones with plant identifiers; some are missed
#    but this gets most of them)
check_accnum <- all_data7[which(nchar(all_data7$acc_num)>9),]
nrow(check_accnum)
no_check_accnum <- setdiff(all_data7,check_accnum)
nrow(no_check_accnum)

# can look at what will be removed in the acc_num
sort(check_accnum[which(grepl("[A-F]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("-[1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("\\.[0][1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("/[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("\\.[0][0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
  # these don't seem to work well with this dataset but may work for others:
#sort(check_accnum[which(grepl("-[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)

# can check which records contain specific elements from above patterns:
#as.data.frame(all_data7[which(grepl("-DCH-18A;-TF-6A",all_data7$acc_num)),])
#as.data.frame(all_data7[which(all_data7$acc_num == "1955-0001*ES"),])

# remove individual-specific identifiers (to combine dup accessions)
#   edit the patterns below based on what you found works well in
#   above tests; each pattern is separated by the "or" (pipe) symbol "|"
check_accnum <- check_accnum %>%
  separate("acc_num","acc_num",
    sep="[A-F]$|-[1-9]$|\\.[0][1-9]$|/[0-9][1-9]$|\\.[0][0-9][1-9]$",
    remove=F) %>%
  group_by(inst_short,acc_num,taxon_name_accepted) %>%
  mutate(num_indiv = sum(as.numeric(num_indiv)),
                germ_type = paste(unique(germ_type),collapse="; "),
                garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  ungroup() %>%
  distinct(inst_short,acc_num,taxon_name_accepted,.keep_all=T)
nrow(check_accnum) #5964

all_data8 <- full_join(check_accnum,no_check_accnum)
nrow(all_data8) #39386

##
## ** ADD Unique ID Column
##

# add unique ID column
  # create UID with institution name, acc num, provenance type, and taxon name
  # also remove duplicates based on new UID and sum individuals
  # now create UID and remove dups
need_id <- all_data8[which(is.na(all_data8$UID)),]
dont_need_id <- all_data8[which(!is.na(all_data8$UID)),]
nms <- names(need_id)
nrow(need_id)
need_id <- need_id %>%
  arrange(orig_lat,locality) %>%
  mutate(UID = paste(inst_short,acc_num,prov_type,taxon_name_accepted,sep="~")) %>%
  group_by(UID) %>%
  mutate(num_indiv = sum(as.numeric(num_indiv))) %>%
  distinct(UID,.keep_all=T) %>%
  ungroup() %>%
  select(c("UID",all_of(nms)))
all_data8 <- rbind(need_id,dont_need_id)
nrow(all_data8) #39367

##
## C) Latitude and Longitude
##

# preserve original lat and long columns
all_data8$lat_dd <- all_data8$orig_lat
all_data8$long_dd <- all_data8$orig_long

# fix some common lat/long character issues
all_data8$lat_dd <- mgsub(all_data8$lat_dd,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
all_data8$long_dd <- mgsub(all_data8$long_dd,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
# replace comma with decimal (european notation)
all_data8$lat_dd <- mgsub(all_data8$lat_dd, c(","), ".")
all_data8$long_dd <- mgsub(all_data8$long_dd, c(","), ".")
# separate values if lat and long both ended up in the lat column
all_data8[which(grepl("\\. ",all_data8$lat_dd)),] <-
  separate(all_data8[which(grepl("\\. ",all_data8$lat_dd)),], col = lat_dd,
           into = c("lat_dd","long_dd"), sep = "\\. ", remove = FALSE)

# replace unwanted characters
  ## latitude
  # replace random unnecessary characters
all_data8$lat_dd <- mgsub(all_data8$lat_dd,
  c("N","\\","/","M","A",": ","E","AZ","R","d","a"," \\."," W")," ")
    # remove leading zero
all_data8$lat_dd[which(grepl("^ *[0][1-9]+",all_data8$lat_dd))] <- gsub(
  "^ *[0]","",all_data8$lat_dd[which(grepl("^ *[0][1-9]+",all_data8$lat_dd))])
all_data8$lat_dd[which(grepl("^S *[0][1-9]+",all_data8$lat_dd))] <- gsub(
  "^S *[0]","-",all_data8$lat_dd[which(grepl("^S *[0][1-9]+",all_data8$lat_dd))])
    # add negative sign if south and remove "S"
all_data8$lat_dd[grep("S",all_data8$lat_dd,ignore.case=T)] <-
  paste("-",all_data8$lat_dd[grep("S",all_data8$lat_dd,ignore.case=T)],sep="")
all_data8$lat_dd <- gsub("S","",all_data8$lat_dd)
all_data8$lat_dd <- gsub("--","-",all_data8$lat_dd)
    # remove double spaces or leading/trailing whitespace
all_data8$lat_dd <- str_squish(all_data8$lat_dd)
sort(unique(all_data8$lat_dd))
  # can check source of specific values that aren't formatted correctly
#all_data8[which(all_data8$lat_dd == "422538"),]
  ## longitude
all_data8$long_dd <- replace_non_ascii(all_data8$long_dd,
  replacement=" ", remove.nonconverted=T)
all_data8$long_dd <- mgsub(all_data8$long_dd,
  c("E","\\","/","NR","d","A","a"," .","o","O")," ")
all_data8$long_dd[which(grepl("^ *[0][1-9]+",all_data8$long_dd))] <- gsub(
  "^ *[0]","",all_data8$long_dd[which(grepl("^ *[0][1-9]+",all_data8$long_dd))])
all_data8$long_dd[which(grepl("^W *[0][1-9]+",all_data8$long_dd))] <- gsub(
  "^W *[0]","-",all_data8$long_dd[which(grepl("^W *[0][1-9]+",
    all_data8$long_dd))])
all_data8$long_dd[grep("W",all_data8$long_dd,ignore.case=T)] <-
  paste("-",all_data8$long_dd[grep("W",all_data8$long_dd,ignore.case=T)],sep="")
all_data8$long_dd <- gsub("W","",all_data8$long_dd)
all_data8$long_dd <- mgsub(all_data8$long_dd,c("--","- "),"-")
all_data8$long_dd <- str_squish(all_data8$long_dd)
sort(unique(all_data8$long_dd))

# convert decimal-minutes-seconds (dms) to decimal degrees (dd)
#   [d, m, and s must be in the same cell, with 1 space between each value]
#   format = ## ## ## (DMS) OR ## ##.### (DM)
  # mark rows that need to be converted
convert <- all_data8[which(grepl(" ",all_data8$lat_dd) |
  grepl(" ",all_data8$long_dd)),]
  nrow(convert) #224
unique(convert$lat_dd)
good <- anti_join(all_data8, convert)
  # separate by dec_min_sec and deg_dec_min then convert to decimal degrees
    # latitude
dms <- convert[which(str_count(convert$lat_dd," ") == 2),]; nrow(dms)
ddm <- convert[which(str_count(convert$lat_dd," ") == 1),]; nrow(ddm)
other <- convert[which((str_count(convert$lat_dd," ") != 1 &
  str_count(convert$lat_dd," ") != 2) | is.na(str_count(convert$lat_dd," "))),]
  nrow(other)
dms$lat_dd = measurements::conv_unit(dms$lat_dd, from = 'deg_min_sec',
  to = 'dec_deg')
ddm$lat_dd = measurements::conv_unit(ddm$lat_dd, from = 'deg_dec_min',
  to = 'dec_deg')
convert <- rbind(dms,ddm,other); nrow(convert)
    # longitude
dms <- convert[which(str_count(convert$long_dd," ") == 2),]; nrow(dms)
ddm <- convert[which(str_count(convert$long_dd," ") == 1),]; nrow(ddm)
other <- convert[which((str_count(convert$long_dd," ") != 1 &
  str_count(convert$long_dd," ") != 2) | is.na(str_count(convert$long_dd," "))),]
  nrow(other)
  dms$long_dd = measurements::conv_unit(dms$long_dd, from = 'deg_min_sec',
    to = 'dec_deg')
  ddm$long_dd = measurements::conv_unit(ddm$long_dd, from = 'deg_dec_min',
    to = 'dec_deg')
  convert <- rbind(dms,ddm,other); nrow(convert)
  # join everything back together
all_data8 <- rbind(good,convert); nrow(all_data8)

# check validity of lat and long
all_data8$lat_dd <- as.numeric(all_data8$lat_dd)
  #sort(unique(all_data8$lat_dd))
all_data8$long_dd <- as.numeric(all_data8$long_dd)
  #sort(unique(all_data8$long_dd))
  # if coords are both 0, set to NA
zero <- which(all_data8$lat_dd == 0 & all_data8$long_dd == 0)
all_data8$lat_dd[zero] <- NA; all_data8$long_dd[zero] <- NA
  # flag non-numeric and not available coordinates and lat > 90, lat < -90,
  # lon > 180, and lon < -180
coord_test <- cc_val(all_data8, lon = "long_dd",lat = "lat_dd",
  value = "flagged", verbose = TRUE) #Flagged 31530 records.
  # try switching lat and long for invalid points and check validity again
all_data8[!coord_test,c("lat_dd","long_dd")] <-
  all_data8[!coord_test,c("long_dd","lat_dd")]
coord_test <- cc_val(all_data8,lon = "long_dd",lat = "lat_dd",
  value = "flagged",verbose = TRUE) #Flagged 29630 records.

# check longitude values that are positive, since
# sometimes the negative gets left off
#pos_long <- all_data8[which(all_data8$long_dd > 0),]
unique(all_data8$taxon_region)
# if the taxon is native to N/S America, make long value negative:
setDT(all_data8)[long_dd > 0 & taxon_region != "Europe" &
                 taxon_region != "Asia" & taxon_region != "Asia & Europe",
                 long_dd := as.numeric(paste0("-",as.character(long_dd)))]

# make coords NA if they are still flagged
coord_test <- cc_val(all_data8,lon = "long_dd",lat = "lat_dd",
  value = "flagged",verbose = TRUE) #Flagged 29630 records.
all_data8[!coord_test,"lat_dd"] <- NA
all_data8[!coord_test,"long_dd"] <- NA
nrow(all_data8) #39367

### now we'll work just with the geolocated points for a bit
all_data8$flag <- ""
have_coord <- all_data8 %>% filter(!is.na(lat_dd) & !is.na(long_dd))
nrow(have_coord) #9737
no_coord <- anti_join(all_data8,have_coord)
no_coord$latlong_country <- ""
# add country-level information to fix potentially-switched lat/longs
  # turn occurrence point data into a spatial object
crdref <- "+proj=longlat +datum=WGS84"
geo_pts_spatial <- vect(cbind(have_coord$long_dd,have_coord$lat_dd),
  atts=have_coord, crs=crdref)
world_polygons <- vect(file.path(main_dir,
  "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))
# add country polygon data to each point based on lat-long location
geo_pts <- terra::intersect(geo_pts_spatial,world_polygons)
# try switching lat and long for points in Antarctica
on_land <- as.data.frame(geo_pts); nrow(on_land)
on_land[which(on_land$COUNTRY == "Antarctica"),c("lat_dd","long_dd")] <-
  on_land[which(on_land$COUNTRY == "Antarctica"),c("long_dd","lat_dd")]
head(on_land)
# remove columns from first join and join to countries again
on_land <- on_land %>% select(UID:flag)
geo_pts_spatial <- vect(cbind(on_land$long_dd,on_land$lat_dd),
  atts=on_land, crs=crdref)
geo_pts <- terra::intersect(geo_pts_spatial,world_polygons)
on_land <- geo_pts %>%
  select(UID:flag,COUNTRY) %>%
  rename(latlong_country = COUNTRY)
on_land <- as.data.frame(on_land)
nrow(on_land) #7285
land_uid <- unique(on_land$UID)
# check if geolocated points are in water and mark
  # get points that have coords but didn't fall in a country;
  # these are in the water
in_water <- have_coord %>% filter(!(UID %in% land_uid))
in_water <- as.data.frame(in_water)
in_water$flag <- "Given lat-long is in water"
in_water$latlong_country <- ""
nrow(in_water) #2452
# bind all the points back together
all_data9 <- rbind(on_land,in_water,no_coord)
nrow(all_data9)

# mark lat-long for records with same inst lat-long and wild lat-long
all_data9$lat_round <- round(all_data9$lat_dd,digits=1)
all_data9$long_round <- round(all_data9$long_dd,digits=1)
all_data9$inst_lat_round <- round(all_data9$inst_lat,digits=1)
all_data9$inst_long_round <- round(all_data9$inst_long,digits=1)
garden_latlong <- all_data9 %>% filter(lat_round == inst_lat_round &
  long_round == inst_long_round & prov_type != "N")
unique(garden_latlong$inst_short)
nrow(garden_latlong)
all_data9[which(all_data9$UID %in% garden_latlong$UID),]$flag <-
  "Given lat-long is at institution, use only if native to grounds"
#all_data9[all_data9$UID %in% garden_latlong$UID,]$lat_dd <- NA
#all_data9[all_data9$UID %in% garden_latlong$UID,]$long_dd <- NA
table(all_data9$flag) #grounds = 3064; water = 2452

# add gps_det (gps determination) column
all_data9$gps_det[which(all_data9$prov_type == "H")] <- "N/A (horticultural)"
all_data9$gps_det[which(!is.na(all_data9$lat_dd) &
  !is.na(all_data9$long_dd))] <- "Given"
all_data9$gps_det[which(all_data9$gps_det == "")] <- NA
table(all_data9$gps_det)
# Given     N/A (horticultural)
# 9737      5989

# where prov_type is "N/A (horticultural)" but lat-long is given, change to "H?"
  # create new prov type column
all_data9$prov_type[which(all_data9$gps_det == "Given" &
  all_data9$prov_type == "H")] <- "H?"
table(all_data9$prov_type)
#    H   H?    N    NG     U     W     Z
# 5989  3157   7 12189  8200  9293   532

##
## D) Collection year
##

# not cleaning up this column right now
sort(unique(all_data9$coll_year))
## IF NEEDED: replace non-year words/characters
#all_data9$coll_year <- mgsub(all_data9$coll_year,
#  c("([0-9]+);","about ","ca.","Unknown","original","Estate","estate"),"")
#all_data9$coll_year[which(all_data9$coll_year == "")] <- NA

# remove extra elements so its just year
#all_data9$coll_year <- gsub(";[1-2][0-9][0-9][0-9]","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[0-9][0-9]/[0-9][0-9]/","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[0-9]/[0-9][0-9]/","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[0-9][0-9]/[0-9]/","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[0-9]/[0-9]/","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[1-9] [A-Z][a-z][a-z] ","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[A-Z][a-z][a-z] ","",all_data9$coll_year)
#all_data9$coll_year <- gsub("[1-9]-[A-Z][a-z][a-z]-","",all_data9$coll_year)

# make column numeric
#all_data9$coll_year <- as.numeric(all_data9$coll_year)
## IF NEEDED: add first two numbers in year
#  # assume 2000s if values is less than 21
#all_data9$coll_year[which(all_data9$coll_year < 10)] <-
#  paste0("200",as.character(all_data9$coll_year[which(all_data9$coll_year < 10)]))
#all_data9$coll_year <- as.numeric(all_data9$coll_year)
#all_data9$coll_year[which(all_data9$coll_year < 21)] <-
#  paste0("20",as.character(all_data9$coll_year[which(all_data9$coll_year < 21)]))
#all_data9$coll_year <- as.numeric(all_data9$coll_year)
#  # assume 1900s if values is greater than or equal to 21
#all_data9$coll_year[which(all_data9$coll_year < 100)] <-
#  paste0("19",as.character(all_data9$coll_year[which(all_data9$coll_year < 100)]))
#all_data9$coll_year <- as.numeric(all_data9$coll_year)
#sort(unique(all_data9$coll_year))

##
## E) Lineage number
##

# remove lin_num when same as acc_num
all_data9[which(all_data9$acc_num == all_data9$lin_num),]$lin_num <- NA

##
## F) Locality
##

# create all_locality column
all_data9$latitude <- round(all_data9$lat_dd,digits=4)
all_data9$longitude <- round(all_data9$long_dd,digits=4)
all_data9 <- unite(all_data9, "all_locality",
  c(locality,municipality,county,state,country,orig_source,
    lin_num,coll_num,coll_name,coll_year,
    latitude,longitude,notes),sep = " | ",remove = F)
# remove NA in concatenated locality column
all_data9$all_locality <- gsub("NA","",all_data9$all_locality)
# if no locality info at all, make it NA
all_data9$all_locality[which(all_data9$all_locality ==
  " |  |  |  |  |  |  |  |  |  |  |  | ")] <- NA

##
## G) Institution type
##

# add inst_type for gene bank data (didn't do this for Quercus)
#all_data9$inst_type[which(is.na(all_data9$inst_type))] <- "Gene/Seed Bank"
#table(all_data9$inst_type)
# Botanic Garden    Botanic Garden; Gene/Seed Bank     Gene/Seed Bank
# 6195              159                                3463
table(all_data9$data_source)
# ex_situ_BG_survey  FAO-WIEWS   Genesys
# 38947              0           420

################################################################################
# 7. Summary statistics
################################################################################

## FINAL VERSION: KEEP ONLY NECESSARY COLUMNS
keep_col <- c(
  # key data
  "UID","inst_short","taxon_name_accepted","acc_num","prov_type","num_indiv",
  # locality
  "lat_dd","long_dd","flag","gps_det","uncertainty","locality",
  "municipality","county","state","country","latlong_country","assoc_sp",
  # source
  "orig_source","lin_num","coll_name","coll_num","coll_year",
  # material info
  "germ_type","garden_loc","rec_as",
  # other metadata
  "notes","filename","submission_year","dataset_year","data_source",
  # taxon name details
  "taxon_full_name_orig","taxon_full_name_concat",
  "genus","hybrid","species","infra_rank","infra_name","cultivar",
  "taxon_name_status","author","taxon_verif",
  # institution metadata
  "inst_country","inst_lat","inst_long",
  # original versions of columns, for reference
  "orig_prov_type","orig_acc_num","orig_lat","orig_long",
  # OPTIONAL additional taxon metadata
  "taxon_region","iucnredlist_category"
)
data_sel <- all_data9[,keep_col]

# write file
write.csv(data_sel, file.path(main_dir,data_out,
  paste0("All_ExSitu_Compiled_", Sys.Date(), ".csv")),row.names = F)

# [OPTIONAL] regional subsets for easier use:
meso <- data_sel %>% filter(grepl("Mesoamerica",taxon_region))
nrow(meso) #6514
write.csv(meso, file.path(main_dir,data_out,
  paste0("Mesoamerican-spp_ExSitu_Compiled_", Sys.Date(), ".csv")),row.names = F)
us <- data_sel %>% filter(grepl("US",taxon_region))
nrow(us) #24029
write.csv(us, file.path(main_dir,data_out,
  paste0("UnitedStates-spp_ExSitu_Compiled_", Sys.Date(), ".csv")),row.names = F)

################################################################################
# One more thing for Quercus...
################################################################################

### Compile 2020 geolocated CSVs for U.S. species into one dataframe for
### eaiser reference

# read in geolocated CSVs
file_list <- list.files(path=file.path(main_dir,data_out,
  "Regional subsets (for faster use)",
  "U.S. Species of Concern - Geolocated CSV files 2020"),
  pattern=".csv",full.names=TRUE)
file_dfs <- lapply(file_list,read.csv,header=TRUE,fileEncoding="LATIN1",
  strip.white=TRUE,stringsAsFactors=F,na.strings=c("","NA"))
# bind geolocated CSVs together
geolocated <- Reduce(rbind.fill, file_dfs)
head(geolocated)

# separate to accession level again
geolocated2 <- separate_rows(geolocated, c(acc_num_CONCAT,UID_CONCAT), sep="; ")
head(as.data.frame(geolocated2))

# keep what's needed and rename columns
geolocated2 <- geolocated2 %>%
  rename(locality = locality.string,
         lat_dd = latitude,
         long_dd = longitude,
         acc_num = acc_num_CONCAT,
         UID = UID_CONCAT,
         taxon_name_accepted = species_name_acc) %>%
  select(UID,inst_short,taxon_name_accepted,acc_num,prov_type,
         lat_dd,long_dd,flag,gps_det,uncertainty,locality,
         coll_year,coll_name,submission_year)

# make gps_det "X" if NA -- not doing this now to distinquish records
geolocated2[which(is.na(geolocated2$gps_det)),]$gps_det <- "X"
head(as.data.frame(geolocated2))

# write file
write.csv(geolocated2, file.path(main_dir,data_out,
  "Regional subsets (for faster use)",
  "UnitedStates-spp_Geolocated_Records_2020.csv"),
  row.names = F)







### STOPPED UPDATING HERE !!! NEXT SECTIONS ARE OLD !!! ###


################################################################################
# 8. (Optionally) Explore geoerferencing needs
################################################################################

### explore georeferencing needs
# table with...
#   species name
#   num wild acc
#   num non-H acc w/ coords
#   num non-H acc with no coords & yes locality info
geo_needs <- data_sel %>%
  group_by(taxon_name_acc) %>%
  summarize(
    num_acc = sum(!is.na(taxon_name_acc)),
    num_wild = sum(prov_type == "W"),
    NotH_YesCoords = sum(!is.na(lat_dd) & prov_type != "H"),
    NotH_NoCoords_YesLocality = sum(is.na(lat_dd) & !is.na(all_locality) & prov_type != "H"),
    Percent_NonH_NeedGeo = (sum(is.na(lat_dd) & !is.na(all_locality) & prov_type != "H") / sum(prov_type != "H")*100)
  )
head(geo_needs,n=20)
# write file
write.csv(geo_needs, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_Geolocation_Needs_Summary_", Sys.Date(), ".csv")),row.names = F)

# records that may need geolocation
#   (no lat-long, yes locality, prov type not H)
need_geo <- data_sel %>%
  filter(is.na(lat_dd) & !is.na(all_locality) & prov_type != "H")
nrow(need_geo) #4282
# add a couple more columns for keeping notes while geolocating
need_geo$uncertainty <- NA
need_geo$geolocated_by <- NA
need_geo$gps_notes <- NA
# condense all_locality duplicates
need_geo <- need_geo %>%
  group_by(prov_type,lat_dd,long_dd,gps_det,all_locality) %>%
  mutate(UID = paste0(UID,collapse=" | "),
         inst_short = paste0(unique(inst_short),collapse=" | "),
         taxon_name_acc = paste0(unique(taxon_name_acc),collapse=" | ")) %>%
  ungroup() %>%
  distinct(UID,inst_short,taxon_name_acc,prov_type,lat_dd,long_dd,uncertainty,
    gps_det,geolocated_by,gps_notes,all_locality,county,state,country) %>%
  select(UID,inst_short,taxon_name_acc,prov_type,lat_dd,long_dd,uncertainty,
    gps_det,geolocated_by,gps_notes,all_locality,county,state,country)
nrow(need_geo) #2511
head(need_geo)
# flag records for species that are high priority...
#   threatened and/or have less than 15 wild accessions
#   (you can choose whatever threshold you want)
  # thresholds
rl_threat <- c("CR (Critically Endangered)","EN (Endangered)","VU (Vulnerable)")
ns_threat <- c("G1 (Critically Imperiled)","G2 (Imperiled)","G3 (Vulerable)")
few_wild <- geo_needs[geo_needs$num_wild<15,]$taxon_name_acc
  # flag records
flag_taxa <- taxon_list %>%
  filter(iucnredlist_category %in% rl_threat |
         natureserve_rank %in% ns_threat |
         taxon_name %in% few_wild) %>%
  select(taxon_name_acc)
flag_taxa$flag <- "Priority"
need_geo <- left_join(need_geo,flag_taxa)
table(need_geo$flag) #368

# write file
write.csv(need_geo, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_Need_Geolocation_", Sys.Date(), ".csv")),row.names = F)

### LINK INSTRUCTIONS HERE FOR GEOLOCATING:
### !!XXXXXXXXXXXXXXXXXXX!!

################################################################################
# 9. Add geolocated data, after manual geolocation
################################################################################

# read in all compiled ex situ data (exported above)
exsitu <- read.csv(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  "ExSitu_Compiled_2022-10-17.csv"), header = T, colClasses="character")

# read in geolocated dataset
geo_raw <- read.csv(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  "ExSitu_Need_Geolocation_2022-10-17_Geolocated.csv"),
  header = T, colClasses="character")
head(geo_raw)

# add geolocated coordinates to ex situ data
  # separate UID row
geolocated <- separate_rows(geo_raw, UID, sep=" \\| ")
#OLD: geolocated <- separate_rows(geo_raw, UID, sep="\\|;\\|")
#OLD: geolocated$UID <- gsub("UCalifornia BGerkeley","UCaliforniaBGBerkeley",geolocated$UID)
  # keep only edited rows (lat, long, gps_det) and
  #   records that have gps_det filled in
geolocated <- geolocated %>%
  select(UID,lat_dd,long_dd,gps_det,uncertainty,geolocated_by,gps_notes) %>%
  filter(gps_det != "")
head(geolocated)
table(geolocated$gps_det)
      #   L    C    G    X
      #  149   36   23   270
  # select geolocated rows in full dataset and remove cols we want to add
exsitu_geo <- exsitu %>%
  filter(UID %in% geolocated$UID) %>%
  select(-lat_dd,-long_dd,-gps_det,-uncertainty)
    # these two values should be the same:
nrow(exsitu_geo)
nrow(geolocated)
  # add geolocation data
exsitu_geo <- full_join(exsitu_geo,geolocated)
  # join geolocated rows with rest of ex situ rows
exsitu_no_geo <- exsitu %>%
  filter(!(UID %in% exsitu_geo$UID))
nrow(exsitu_no_geo)
exsitu_all <- rbind.fill(exsitu_no_geo,exsitu_geo)
nrow(exsitu_all)
table(exsitu_all$gps_det)

# write new file
write.csv(exsitu_all, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_Compiled_Post-Geolocation_", Sys.Date(), ".csv")), row.names = F)
# write to in situ folder also
write.csv(exsitu_all, file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ",
  paste0("ExSitu_Compiled_Post-Geolocation_", Sys.Date(), ".csv")), row.names = F)




# Used the above geolocation process, not the GEOLocate website......

################################################################################
# 9. (Optionally) Create file for GEOLocate
################################################################################

# read in data in again, if you didn't just run the whole script
#need_geo <- read.csv(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
#  "ExSitu_Need_Geolocation_2022-10-14.csv"), header = T, colClasses="character")

# add GEOLocate standard columns
need_geo$correction.status <- NA
need_geo$precision <- NA
need_geo$error.polygon <- NA
need_geo$multiple.results <- NA
need_geo$latitude <- NA
need_geo$longitude <- NA

# update column name and order for GEOLocate requirements
geolocate <- need_geo %>%
  # rename to GEOLocate standard columns
  rename(locality.string = all_locality) %>%
  # order records alphabetically
  arrange(taxon_name_acc,locality.string) %>%
  # reorder columns
  select(
    ## GeoLocate
    locality.string,country,state,county,latitude,longitude,
    correction.status,precision,error.polygon,multiple.results,uncertainty,
    ## metadata to include also
    taxon_name_acc,inst_short,prov_type,gps_det,geolocated_by,gps_notes,UID) %>%
  # replace NA with "" to make simpler to view in GEOLocate
  replace(., is.na(.), "")
head(geolocate)

# remove dot in column names (replace with space) for GEOLocate
names(geolocate) <- gsub(x = names(geolocate),pattern = "\\.",replacement = " ")
str(geolocate)
head(as.data.frame(geolocate))

# split by species, so each species has separate file;
#   this step makes it easier to use the GEOLocate tool because it doesn't
#   handle lots of records very well

# create one CSV for each target species
# !!! this doesn't work well if you combined species above to remove locality duplicates !!!
sp_split <- split(geolocate, as.factor(geolocate$taxon_name_acc))
names(sp_split) <- gsub(" ","_",names(sp_split))

# write files
if(!dir.exists(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R","files_for_GEOLocate")))
  dir.create(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R","files_for_GEOLocate"),recursive=T)
lapply(seq_along(sp_split), function(i) write.csv(sp_split[[i]],
  file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R","files_for_GEOLocate",
  paste0(names(sp_split)[[i]],".csv")),row.names = F))


### !!!
### ! NOW GO GEOLOCATE !
### !!!



# Other random old code chunk:

### TESTING GEONAMES PACKAGE ###

# read in compiled ex situ data
df <- read.csv(file.path(main_dir,"outputs",
  "exsitu_compiled_standardized_2021-02-11_firstpassGpsDet.csv"),
  header = T, colClasses="character")
str(df)
unique(df$gps_det)
tail(df[which(df$gps_det == "L"),])

# create login: https://www.geonames.org/login
# https://www.geonames.org/manageaccount
#   the account is not yet enabled to use the free web services. Click here to enable.
usethis::edit_r_environ()
devtools::install_github("ropensci/geonames")
library(geonames)

# https://rstudio-pubs-static.s3.amazonaws.com/489236_0259d8532a354ad6945d818bc4c052f1.html

login <- read_lines(log_loc)
username  <- login[4]
options(geonamesUsername = username)
head(GNsearch(name = "Chicago"))
#  adminCode1       lng geonameId               toponymName countryId fcl population countryCode                      name           fclName adminCodes1.ISO3166_2   countryName                                      fcodeName adminName1      lat fcode
#1         IL -87.65005   4887398                   Chicago   6252001   P    2720546          US                   Chicago city, village,...                    IL United States seat of a second-order administrative division   Illinois 41.85003 PPLA2
#2         IL -87.89062   6955104 Chicago-Naperville-Joliet   6252001   L    9569624          US Chicago-Naperville-Joliet   parks,area, ...                    IL United States                                economic region   Illinois 41.70778  RGNE
#3         IL  -87.6979  12070033 Chicago metropolitan area   6252001   L    9472676          US Chicago metropolitan area   parks,area, ...                    IL United States                                         region   Illinois  41.8902   RGN
#4         OH -82.72629   5176830                   Willard   6252001   P       6063          US                   Willard city, village,...                    OH United States                                populated place       Ohio 41.05311   PPL
#5         IL -87.69644   4887463              Chicago Lawn   6252001   P      55551          US              Chicago Lawn city, village,...                    IL United States                     section of populated place   Illinois 41.77503  PPLX
#6         IL -87.55425   4911863             South Chicago   6252001   P      28095          US             South Chicago city, village,...                    IL United States                                populated place   Illinois 41.73977   PPL

lanc_coords <- lanc_df[1, c("lng", "lat")]
