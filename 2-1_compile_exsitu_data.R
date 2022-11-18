################################################################################

## 2-1_compile_exsitu_data_CWRcopy.R

### Author: Emily Beckman Bruns
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Smaller edits were added with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 30 March 2022
### Last updated: 26 May 2022

### R version 4.1.3

### DESCRIPTION:
  #
  # This script takes a folder of CSV files representing accessions data from
  #   different institutions, combines them into one dataset, and standardizes
  #   some important fields.
  #

### DATA IN:
  #
  # 1. Folder ("exsitu_standard_column_names") of CSV files whose column names
  #     have already be standardized by hand using the
  #     "standardizing_accessions_data_fields.xlsx" template [**ADD TO REPO**]
  #
  # 2. Table with metadata for institutions who provided accessions data
  #     (respondent_institution_data_table.csv); required columns include:
  #     | inst_short                               | inst_id         | inst_type      | inst_lat  | inst_long | inst_country |date_data_received  |
  #     |------------------------------------------|-----------------|----------------|-----------|-----------|--------------|--------------------|
  #     |nickname for institution (see template    |institution      |institution type|institution|institution|institution   |date the accessions |
  #     |standardizing_accessions_data_fields.xlsx)|identifier from  |(Botanic Garden,|latitude   |longitude  |country       |data were received  |
  #     |                                          |BGCI GardenSearch| Gene/Seed Bank)|           |           |              |from the institution|
  #
  # 3. Target taxa list (target_taxa_with_synonyms.csv); required columns include:
  #     | taxon_name_acc | taxon_name              | genus | genus_species      | taxon_name_status  |
  #     |----------------|-------------------------|-------|--------------------|--------------------|
  #     |accepted taxon  |taxon name (either same  |genus  |species name in the |status of taxon_name|
  #     |name            |as accepted or a synonym)|name   |form 'Genus species'|(Accepted, Synonym) |
  #
  # 4. Accession-level data downloads from international crop genebank databases
  # (a) Genesys [Global Crop Diversity Trust]: Download data
  #       <https://www.genesys-pgr.org/a/overview> for target genera and name
  #       the folder "genesys-accessions"
  # (b) WIEWS [FAO's World information and early warning system on plant genetic
  #     resources for food and agriculture]: Download data
  #       <https://www.fao.org/wiews/data/ex-situ-sdg-251/search/en/?no_cache=1>
  #       for target genera, and name the file "Wiews_Exsitu.csv"
  #
  # 5. Polygons...

### DATA OUT:
  # ExSitu_AllDataRaw_GenusFilterOnly_[[Sys.Date]].csv
  #

################################################################################
# Load libraries
################################################################################

#rm(list=ls())
my.packages <- c('plyr','tidyverse', 'data.table', 'textclean',
                 'measurements', 'CoordinateCleaner','rnaturalearth',
                 'rnaturalearthdata','maps','raster','spatialEco','leaflet')
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

select <- dplyr::select

################################################################################
# Set working directory
################################################################################

# Google Drive folder with accession data
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

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

# create buffers around points, using specified projection
create.buffers <- function(df,radius,pt_proj,buff_proj){
	# select coordinate columns
	latlong <- df %>% select(long_dd,lat_dd)
	# turn occurrence point data into a SpatialPointsDataFrame
	sp_df <- SpatialPointsDataFrame(latlong, df, proj4string = pt_proj)
	# reproject SpatialPointsDataFrame to specified projection
	proj_df <- spTransform(sp_df,buff_proj)
	# place buffer around each point
	buffers <- buffer(proj_df,width=radius,dissolve=T)
	# return buffer polygons
	return(buffers)
}

# map individual species with raster data and ex situ points + buffers
map_sp <- function(wild_dist,state_bound,sp_pts){
    # create palette
  pal <- colorNumeric(c("#6ca874"), values(wild_dist), na.color = "transparent") #57945f
    # get CRS
  wgs.proj <- wild_dist@crs
    # create map
  map <- leaflet() %>%
  	## background
  	addProviderTiles("CartoDB.PositronNoLabels",#"Esri.WorldGrayCanvas","Esri.WorldShadedRelief","Esri.WorldTerrain",
  		options = providerTileOptions(maxZoom = 10)) %>%
  	## EPA Level IV ecoregions
  	#addPolygons(
  	#	data = ecol4_sel,
  	#	fillColor = ~ecol4_pal(ecol4_sel@data$US_L4NAME), fillOpacity = 0.6,
  	#	color = ~ecol4_pal(ecol4_sel@data$US_L4NAME), weight = 0.5, opacity = 1) %>%
  	## in situ distribution (raster)
    addRasterImage(wild_dist, colors = pal, opacity = 0.7) %>%
  	## state boundaries
  	addPolygons(data = state_bound,
  		fillOpacity = 0, color = "#969696", weight = 1.5, opacity = 1) %>%
    ## ex situ buffers
  	addPolygons(data = create.buffers(sp_pts,50000,wgs.proj,wgs.proj),
  		smoothFactor = 0.5,	weight = 1.5, opacity = 1, color = "#3e146e",
  		fillOpacity = 0.3) %>%
  	## in situ points
  	#addCircleMarkers(data = insitu,
  	#	lng = ~decimalLongitude, lat = ~decimalLatitude,
  	#	#popup = ~paste("Source(s):", all_source_databases, UID),
  	#	radius = 4, fillOpacity = 1, stroke = F, color = "#1c1c1b") %>%
  	## ex situ points
  	addCircleMarkers(data = sp_pts,
  		lng = ~long_dd, lat = ~lat_dd,
  		#popup = ~paste("Garden:", datasetName, UID),
  		radius = 4, fillOpacity = 1, stroke = F, color = "#3e146e") %>%
    ## title
  	#addControl("In situ distribution and representation in ex situ collections"),
  	#	position = "topright") %>%
  	## legend
  	addLegend(labels =
  		c("Predicted wild distribution (PNAS 2020)",
  			"Populations sampled for ex situ (wild collection locations and 50 km buffers)"),
  		colors = c("#6ca874","#3e146e"), title = "Legend",
  		position = "bottomright", opacity = 0.8) %>%
  	## scale bar
  	#addScaleBar(position = "bottomleft",
  	#	options = scaleBarOptions(maxWidth = 150)) %>%
  	## pick coords for center of frame and zoom level when map first appears
  	setView(-97, 40, zoom = 5)
  return(map)
}

################################################################################
################################################################################
# 1. Read in and stack all ex situ accessions data from survey of gardens
################################################################################

### FIRST: After receiving accession data from institutions, you need to process
#   it manually. See here for instructions:
#   https://docs.google.com/spreadsheets/d/1p5HAS7vIE-3CbQcUmwrnuBGv5324dXy-42Iu6LlbX0E/edit?usp=sharing

## read in data from multiple surveys and stack, or just read in from one folder
##    this function also adds columns for 1) the file name [often equivalent to the
##    "inst_short" institution nickname] 2) a sumbission year, 3) an accession
##    number if one isn't given
#raw_2020 <- read.exsitu.csv(file.path(main_dir,"inputs",
#  "exsitu_standard_column_names","data_2020"), "2020")
raw_2022 <- read.exsitu.csv(file.path(main_dir,"ex-situ_data",
  "exsitu_standard_column_names"), "2022")
# stack all data if you had multiple years:
#to_stack <- list(raw_2020,raw_2022)
#all_data <- Reduce(rbind.fill,to_stack)

# create new version before big changes, so can easily go back to original
all_data <- raw_2022

# replace all non-ascii characters
all_data <- as.data.frame(lapply(all_data,replace_non_ascii),stringsAsFactors=F)

# check out column names
sort(colnames(all_data))
## IF NEEDED: separate column into multiple
#all_data <- all_data %>% separate("specific",
#  c("infra_rank_add","infra_name_add"),sep=" ",remove=T,fill="right")
## IF NEEDED: see which datasets have extraneous columns so you can fix manually
##  in the raw data as desired; update line below
#unique(all_data$filename[all_data$locality.1 !=""])
## IF NEEDED: merge similar columns (you may not need to do this if no schema
##  mistakes were made when manually editing column names)
all_data <- tidyr::unite(all_data,"cultivar", c("cultivar","cultivsr"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"genus", c("genus","gensu"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"hybrid", c("hybrid","hyrbid"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"notes", c("notes","Notes","notes2"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"taxon_full_name", c("taxon_full_name","taxon_name_full","taxon.full_name"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"taxon_verif", c("taxon_verif","Taxonomic.Verification","taxon_det"),sep="; ",remove=T,na.rm=T)
all_data <- tidyr::unite(all_data,"rec_year", c("rec_year","planted_year"),sep="; ",remove=T,na.rm=T)

## IF NEEDED: remove unused columns or rename columns
  unique(all_data$condition) # make sure nothing is "dead" or "removed"
all_data <- all_data[ , -which(names(all_data) %in%
  c("condition","infra_author","viable_indiv","voucher"))]

### CHECK THINGS OUT ###
sort(colnames(all_data)); ncol(all_data)
# There should be max of 36 columns, including no more than:
  # acc_num,assoc_sp,author,coll_date,coll_name,coll_num,coll_year,country,
  # county,cultivar,filename,garden_loc,genus,germ_type,hybrid,infra_name,
  # infra_rank,inst_short,inst_short2,lin_num,locality,municipality,notes,
  # num_indiv,orig_lat,orig_long,orig_source,prov_type,rec_as,rec_date,
  # rec_year,species,state,submission_year,taxon_full_name,taxon_verif

# fill in inst_short column with filename if none provided
all_data$inst_short[is.na(all_data$inst_short)] <-
  all_data[is.na(all_data$inst_short),]$filename
nrow(all_data) #26463
# remove rows with no inst_short
all_data <- all_data %>% filter(inst_short!="" & inst_short!="ManualEntry")
nrow(all_data) #25651

### CHECK ALL INSTITUTIONS ARE HERE ###
sort(unique(all_data$inst_short)) #197

# remove leading, trailing, and middle (e.g., double space) whitespace,
#   to prevent future errors
all_data <- as.data.frame(lapply(all_data, function(x) str_squish(x)),
  stringsAsFactors=F)
# replace "" cells with NA in whole dataset
all_data[all_data == ""] <- NA

# add data source colun
all_data$data_source <- "ex_situ_BG_survey"

nrow(all_data) #25651


##### REMOVE DUPLICATE DATA
#remove data for rows in the respondent_institution_data_table that say “No Quercus” (based on inst_short and survey_year)
#if there is newer data from 2019 or 2020, I’ll remove the old data from 2017 and/or 2019
#for 2021 and 2022, I’ll only remove older data for target species, based on this list (Mesoamerican sp only) for 2021 data and this list (Mesoamerican sp plus US sp of concern) for 2022 data. Therefore it’s important that the 2021 data and 2022 data are in the correct folder/ labelled correctly in the survey_year column of the respondent_institution_data_table
#if there is PCN data from the same year that the garden provided their own data, I’ll remove the PCN data and keep what the garden provided directly



################################################################################
# 2. Compile Genesys data
################################################################################

# field metadata: https://www.genesys-pgr.org/documentation/basics

# load in data
  # collection site description
genPath <- paste0(main_dir,"/ex-situ_data/genesys-accessions/coll.csv")
gen_col <- data.table::fread(file = genPath, header = TRUE)
  str(gen_col)
  nrow(gen_col) ; length(unique(gen_col$genesysId)) #49998
  # metadata about record
genPath <- paste0(main_dir,"/ex-situ_data/genesys-accessions/core.csv")
gen_core <- data.table::fread(genPath, header = TRUE)
  str(gen_core)
  nrow(gen_core) ; length(unique(gen_core$genesysId)) #78432
  # spatial info
  #   for some reason the headers dont read in correctly with fread;
  #   will use read.csv just to get the headers and add them to the fread df
genPath <- paste0(main_dir,"/ex-situ_data/genesys-accessions/geo.csv")
gen_geo <- data.table::fread(genPath, header = TRUE)
  str(gen_geo)
  colnames(gen_geo)
gen_geo_head <- read.csv(genPath)
  colnames(gen_geo_head)
gen_geo <- gen_geo %>% select(-V8) # remove the last column (nothing in it)
colnames(gen_geo) <- colnames(gen_geo_head)
  str(gen_geo)
  nrow(gen_geo) ; length(unique(gen_geo$genesysId)) #78431

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

# rename columns to match ex situ data and select only those we need
genesys_sel <- genesys %>%
  dplyr::rename(taxon_full_name = fullTaxa,
         inst_short2 = duplSite,
         coll_num = collNumb,
         coll_date = collDate,
         locality = collSite,
         notes = uuid,
         inst_short = instCode,
         acc_num = acceNumb,
         infra_rank = subtaxa,
         country = origCty,
         rec_date = acqDate,
         orig_lat = latitude,
         orig_long = longitude,
         gps_det = method,
         germ_type = storage,
         orig_source = collSrc,
         prov_type = sampStat) %>%
  select(UID,taxon_full_name,inst_short2,coll_num,coll_date,locality,notes,
    inst_short,acc_num,infra_rank,country,rec_date,orig_lat,orig_long,gps_det,
    germ_type,orig_source,prov_type,genus,species,uncertainty)

# filter by target taxa & add taxon data
## We do this later!
#genesys_sel <- left_join(genesys_sel,taxon_list)
#genesys_target <- genesys_sel %>% filter(!is.na(taxon_name_acc))
#nrow(genesys_target) #5233

# US institutions covered by GRIN (I think); we will remove these
#   pulled this list from CWR Gap Analysis 2020 (Khoury et al.)...
USDAcodes <- c("USA003" ,"USA004", "USA005" ,"USA016" ,"USA020",
"USA022", "USA026", "USA028", "USA029", "USA042" ,"USA047", "USA049",
"USA074", "USA108", "USA133", "USA148", "USA151", "USA167", "USA176",
 "USA390", "USA955", "USA956", "USA970", "USA971", "USA995")
genesys_target <- genesys_sel %>% filter(!(inst_short %in% USDAcodes))
nrow(genesys_target) #63935

# add data source colun
genesys_target$data_source <- "Genesys"

# join to exsitu data
all_data <- rbind.fill(all_data,genesys_target)
nrow(all_data) #89586

################################################################################
# 3. Compile WIEWS data
################################################################################

# read in
wiewsPath <- paste0(main_dir,"/ex-situ_data/Wiews_Exsitu.csv")
#wiews <- data.table::fread(wiewsPath, header = TRUE)
wiews <- read.csv(wiewsPath)
str(wiews) #95551

# filter out Genesys data
wiews <- wiews %>%
  filter(Source.of.information != "Genesys (https://www.genesys-pgr.org)")
nrow(wiews) #81051

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
  dplyr::rename(inst_short = Holding.institute.code,
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
nrow(wiews_sel) #81051

# add data source colun
wiews_sel$data_source <- "FAO-WIEWS"

# join to exsitu data
all_data <- rbind.fill(all_data,wiews_sel)
nrow(all_data) #170637

################################################################################
# 4. Save raw output of compiled data for target genera
#     (can be used to look for hybrids/cultivars, which are removed in next step)
################################################################################

all_data2 <- all_data
nrow(all_data2) #170637

# preserve original taxon name
all_data2$taxon_full_name_orig <- all_data2$taxon_full_name

# fill genus column if not already filled
  # this warning is ok: "Expected 1 pieces. Additional pieces discarded..."
all_data2 <- all_data2 %>% separate("taxon_full_name","genus_temp",sep=" ",remove=F)
all_data2[which(is.na(all_data2$genus)),]$genus <- all_data2[which(is.na(all_data2$genus)),]$genus_temp
# standardize capitalization
all_data2$genus <- str_to_title(all_data2$genus)

### MAKE SURE NO GENUS MISSPELLINGS OR ABBREVIATIONS ###
sort(unique(all_data2$genus))
# make corrections if necessary:
#all_data2$genus <- mgsub(all_data2$genus,
#  c("^Q$","Querucs","Cyclobalanopsis"),"Quercus",fixed=F)

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "target_taxa_with_synonyms.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
str(taxon_list) #284

# remove rows not in target genus/genera
target_genera <- unique(taxon_list$genus)
all_data3 <- all_data2 %>% filter(genus %in% target_genera)
nrow(all_data2); nrow(all_data3) #170637 ; 157687

### CHECK OUT THE HYBRID COLUMN ###
sort(unique(all_data3$hybrid))
# standardize a bit
all_data3$hybrid <- mgsub(all_data3$hybrid,
  c(" A ","^A ","^A$"," X ","^X ","^X$"," _ ","^_ ","^_$"),
  " x ", fixed=F)
all_data3$hybrid <- str_squish(all_data3$hybrid)
# make sure everything has an " x " in it somewhere or that the cell starts
#   with "x" (important later)
sort(unique(all_data3$hybrid))
all_data3$hybrid <- mgsub(all_data3$hybrid,
  c("M. ioensisxpumila","M. baccataxspectabilis",
      "M. atrosanguineaxM. pumilavar.niedzwetzkyana","M. 'Dorothea'"),
  c("M. ioensis x pumila","M. baccata x spectabilis",
      "M. atrosanguinea x M. pumila var. niedzwetzkyana"," x M. 'Dorothea'" ))

# create concatenated taxon_full_name column
all_data3 <- tidyr::unite(all_data3, "taxon_full_name_concat",
  c(genus,hybrid,species,infra_rank,infra_name,cultivar), sep=" ", remove=F,
  na.rm=T)

# when blank, fill taxon_full_name column with concatenated full name
all_data3[is.na(all_data3$taxon_full_name),]$taxon_full_name <-
  all_data3[is.na(all_data3$taxon_full_name),]$taxon_full_name_concat
unique(all_data3$taxon_full_name)

# standardize common hybrid signifiers in taxon_full_name
all_data3$taxon_full_name <- mgsub(all_data3$taxon_full_name,
  c(" A ","A ","_"," X ","X ","\\*","×"," х ")," x ",fixed=T)
# make sure the author is separated in taxon_full_name
all_data3$taxon_full_name <- gsub("\\("," (",all_data3$taxon_full_name)
all_data3$taxon_full_name <- str_squish(all_data3$taxon_full_name)
unique(all_data3$taxon_full_name)

# create folder for output data
if(!dir.exists(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R")))
  dir.create(file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R"), recursive=T)

# write copy of all data
  # select columns
all_data_export1 <- all_data3 %>%
  select(inst_short,inst_short2,taxon_full_name,genus,species,infra_rank,
    infra_name,hybrid,cultivar,taxon_full_name_orig,prov_type,orig_lat,orig_long,
    country,state,municipality,locality,assoc_sp,acc_num,lin_num,orig_source,
    rec_as,rec_year,rec_date,num_indiv,germ_type,garden_loc,coll_num,coll_name,
    coll_year,coll_date,taxon_verif,notes,filename,submission_year,data_source)
write.csv(all_data_export1, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_AllDataRaw_GenusFilterOnly_", Sys.Date(), ".csv")),row.names = F)

# summary of genera for each institution
gen_summary <- all_data3 %>%
  arrange(genus) %>%
  rename(genera = genus) %>%
  group_by(inst_short) %>%
  mutate(
    genera = paste(unique(genera), collapse = '; ')) %>%
  ungroup() %>%
  select(inst_short,genera) %>%
  distinct(inst_short,.keep_all=T)
head(as.data.frame(gen_summary))
write.csv(gen_summary, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("Genera_Reported_by_Institutions_", Sys.Date(), ".csv")),row.names = F)

################################################################################
# 3. Further standardize taxon name, then keep data for target taxa only
#     (removes hybrids and cultivars without specific epithet)
################################################################################

## FINISH STANDARDIZING TAXON NAMES

# add space after periods in taxon_full_name
all_data3$taxon_full_name <- gsub(".",". ",all_data3$taxon_full_name,fixed=T)
all_data3$taxon_full_name <- str_squish(all_data3$taxon_full_name)

## (OPTIONALLY) REMOVE HYBRIDS

# if you DO want hybrids:
all_data4 <- all_data3
# if you DON'T want hybrids:
# remove hybrids based on " x " in taxon_full_name_concat and/or taxon_full_name
#all_data4 <- all_data3 %>% filter(!grepl(" x ",taxon_full_name_concat))
#nrow(all_data3); nrow(all_data4) #45057 ; 43363
#all_data4 <- all_data4 %>% filter(!grepl(" x ",taxon_full_name))
#nrow(all_data4) #41936
# see hybrids removed:
#sort(unique(anti_join(all_data3,all_data4)$taxon_full_name))

## CREATE NEW SEPARATED TAXON NAME COLUMNS

# fixing some taxon name issues I noticed:
all_data4$taxon_full_name <- mgsub(all_data4$taxon_full_name,
  c("Corylus Acolurnoides","Prunus Ayedoensis","Pyrus (Malus) fusca"),
  c("Corylus x colurnoides","Prunus x yedoensis","Malus fusca"))
## first change hybrid notation temporarily: remove space so it stays together
all_data4$taxon_full_name <- gsub(" x "," +",all_data4$taxon_full_name)
# separate out taxon full name and trim whitespace again
    # this warning is ok: "Expected 1 pieces. Additional pieces discarded..."
all_data4 <- all_data4 %>% separate("taxon_full_name",
  c("genus_new","species_new","extra1","extra2",
    "extra3","extra4","extra5","extra6","extra7"),sep=" ",extra="warn",
    remove=F,fill="right")
all_data4 <- as.data.frame(lapply(all_data4,str_squish),stringsAsFactors=F)
# replace genus_new with genus, since we fixed that up in the previous section
all_data4$genus_new <- all_data4$genus

## REMOVE RECORDS WITHOUT SPECIES NAME

# remove records with no/non-standard specific epithet (mostly cultivars or 'sp.')
#   (by looking in species name column)
all_data5 <- all_data4 %>%
  dplyr::filter(!grepl("\"",species_new) &
                !grepl("\'",species_new) &
                !grepl("\\[",species_new) &
                !grepl("\\(",species_new) &
                !grepl("\\.",species_new) &
                !grepl("[A-Z]",species_new) &
                !grepl("[0-9]",species_new) &
                !grepl("\\?",species_new) &
                !is.na(species_new))
nrow(all_data5) #150636
# see records removed; can add anything you want to fix to the
#   "fixing some taxon name issues I noticed" section above:
sort(unique(anti_join(all_data4,all_data5)$taxon_full_name))

## FIND INFRATAXA

## look for infrataxa key words
# make data in all "extra" columns lower case
sp_col <- grep("^species_new$", colnames(all_data5))
all_data5[,sp_col:(sp_col+5)] <- as.data.frame(sapply(
  all_data5[,sp_col:(sp_col+5)], tolower), stringsAsFactors=F)
# create matrix of all "extra" species name columns, to search for
#   infraspecific key words
search.col <- matrix(cbind(all_data5$extra1,all_data5$extra2,all_data5$extra3,
  all_data5$extra4,all_data5$extra5,all_data5$extra6,all_data5$extra7),
  nrow=nrow(all_data5))
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
all_data5$infra_rank_new <- NA
all_data5$infra_rank_new[matches_i] <- all_data5[matches_i]
#unique(all_data5$infra_rank_new) # check results

# create new infra_name column and fill with next column over from "extra"
#   contents that matched infraspecific key word
all_data5$infra_name_new <- NA
matches_i[,2] <- matches_i[,2]+1
all_data5$infra_name_new[matches_i] <- all_data5[matches_i]
#sort(unique(all_data5$infra_name_new))

# standardize infraspecific rank names
all_data5$infra_rank_new <- replace(all_data5$infra_rank_new,
  grep("^v$|^v.$|^var$|^variety$|^va$",all_data5$infra_rank_new), "var.")
all_data5$infra_rank_new <- replace(all_data5$infra_rank_new,
  grep("^subspecies$|^subsp$|^ssp$|^ssp.$|^subs.$|^spp.$|^sub.$",
  all_data5$infra_rank_new), "subsp.")
all_data5$infra_rank_new <- replace(all_data5$infra_rank_new,
 grep("^forma$|^form$|^fma$|^fo$|^fo.$|^f$",all_data5$infra_rank_new), "f.")
unique(all_data5$infra_rank_new)

## CREATE FINAL TAXON FULL NAME FOR FILTERING

# create new taxon full name column
all_data5$taxon_full_name <- NA
  # select rows with infraspecific name and concatenate
yes_infra <- which(!is.na(all_data5$infra_rank_new) &
  !is.na(all_data5$infra_name_new))
all_data5$taxon_full_name[yes_infra] <- paste(all_data5$genus_new[yes_infra],
  all_data5$species_new[yes_infra], all_data5$infra_rank_new[yes_infra],
  all_data5$infra_name_new[yes_infra],sep=" ")
  # select rows without infraspecific name and concatenate
all_data5$taxon_full_name[-yes_infra] <- paste(all_data5$genus_new[-yes_infra],
  all_data5$species_new[-yes_infra],sep=" ")
# check out results
sort(unique(all_data5$taxon_full_name))
# switch hybrid symbol back to " x "
all_data5$taxon_full_name <- gsub(" \\+"," x ",all_data5$taxon_full_name)
sort(unique(all_data5$taxon_full_name))

## FILTER OUT NON-TARGET TAXA

# rename some taxon name columns to preserve originals
all_data5 <- all_data5 %>%
  dplyr::rename(taxon_name = taxon_full_name,
                genus_orig = genus,
                species_orig = species,
                infra_rank_orig = infra_rank,
                infra_name_orig = infra_name)
all_data5 <- all_data5 %>%
  dplyr::rename(genus = genus_new,
                species = species_new,
                infra_rank = infra_rank_new,
                infra_name = infra_name_new)

# join dataset to taxa list
  # join by taxon name
taxon_list <- taxon_list %>% select(-genus)
all_data6 <- left_join(all_data5,taxon_list)
  # if no taxon match, join again just by species name
need_match <- all_data6[which(is.na(all_data6$taxon_name_status)),]
  nrow(need_match) #139666
    # remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data6)-ncol(taxon_list)+1)]
    # remove taxon_name col from taxon data so it doesn't match
taxon_list_sp <- taxon_list %>% select(-taxon_name)
  # create genus_species column
need_match$genus_species <- paste(need_match$genus,need_match$species)
    # new join by genus_species
need_match <- left_join(need_match,taxon_list_sp)
  # bind together new matches and previously matched
matched <- all_data6[which(!is.na(all_data6$taxon_name_status)),]
all_data6 <- rbind(matched,need_match)
  table(all_data6$taxon_name_status) # Accepted: 12341 | Synonym: 1487
  # see how many rows have taxon name match
nrow(all_data6[which(!is.na(all_data6$taxon_name_status)),]) #13828

### CHECK UNMATCHED SPECIES, TO ADD TO SYNONYM LIST AS NECESSARY ###
check <- all_data6 %>% filter(is.na(taxon_name_status))
check <- data.frame(taxon_name = sort(unique(check$taxon_name)))
head(check)
# write file for checking, as desired
  # IF YOU FIND MISSPELLINGS AND/OR ADDITIONAL SYNONYMS, YOU CAN ADD THEM TO
  #   YOUR TARGET TAXA LIST AND GO BACK TO THE LINE WHERE WE "read in target
  #   taxa list" (SECTION 2) AND RUN AGAIN FROM THERE
write.csv(check, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  "ExSitu_UnmatchedSpecies.csv"),row.names = F)

# keep only matched names
all_data7 <- all_data6 %>% dplyr::filter(!is.na(taxon_name_status))
nrow(all_data7) #13828

################################################################################
# 4. Standardize important columns
################################################################################

# keep only necessary columns
all_data8 <- all_data7 %>% select(
  # key data
  UID,inst_short,inst_short2,taxon_name_acc,species_name_acc,prov_type,
  num_indiv,acc_num,
  # locality
  orig_lat,orig_long,locality,municipality,county,state,country,assoc_sp,
  gps_det,uncertainty,
  # source
  orig_source,lin_num,coll_num,coll_name,coll_year,coll_date,
  # material info
  germ_type,garden_loc,rec_as,rec_year,rec_date,
  # other metadata
  notes,filename,submission_year,data_source,taxon_name_status,fruit_nut,
  iucnredlist_category,natureserve_rank,
  # taxon name details
  taxon_name,genus,species,infra_rank,infra_name,cultivar,hybrid,
  taxon_full_name_orig,taxon_full_name_concat,taxon_verif)

# add institution metadata
inst_data <- read.csv(file.path(main_dir,"ex-situ_data",
  "respondent_institution_data_table.csv"), stringsAsFactors = F)
str(inst_data)
all_data9 <- left_join(all_data8,inst_data)
str(all_data9)

##
## Provenance type
##

# look at column contents and change below phrases as needed
all_data9$orig_prov_type <- all_data9$prov_type
all_data9$prov_type <- str_to_lower(all_data9$prov_type)
sort(unique(all_data9$prov_type))
## IF NEEDED: transfer contents of one column to another column, if data
#   needs to be preserved but is in wrong place
all_data9[which(all_data9$prov_type=="130"|grepl("^130) ",all_data9$prov_type)),]$notes <- "Semi-natural/sown"
all_data9[which(all_data9$prov_type=="410"|grepl("^410) ",all_data9$prov_type)),]$notes <- "Breeding/research material: Breeder's line"
all_data9[which(all_data9$prov_type=="300"|grepl("^300) ",all_data9$prov_type)),]$notes <- "Traditional cultivar/landrace"
all_data9[which(all_data9$prov_type=="400"|grepl("^400) ",all_data9$prov_type)),]$notes <- "Breeding/research material"
all_data9[which(all_data9$prov_type=="500"|grepl("^500) ",all_data9$prov_type)),]$notes <- "Advanced or improved cultivar (conventional breeding methods)"

# standardize column by searching for keywords and replacing with standard value
  # remove confusing words/phrases
all_data9$prov_type <- mgsub(all_data9$prov_type,
  c("not of known wild origin","could be cultivated"), "")
  # ex wild (Z)
all_data9$prov_type <- ifelse(grepl(paste(
  c("indirect","ex wild","^z$","cultivated from wild"),
  collapse = "|"), all_data9$prov_type),"Z",all_data9$prov_type)
  # wild (W)
all_data9$prov_type <- ifelse(grepl(paste(
  c("wild","wld","collect","^w$","^w\\*","^\\(w\\)$","wd","w\\?","genetic",
    "100","110","130"),
  collapse = "|"), all_data9$prov_type),"W",all_data9$prov_type)
  # native to site (N)
all_data9$prov_type <- ifelse(grepl(paste(
  c("original to site","spontaneous","^n$"),
  collapse = "|"), all_data9$prov_type),"N",all_data9$prov_type)
  # unknown (U)
all_data9$prov_type <- ifelse(grepl(paste(
  c("^\\(u\\)$","^u$","^u ","unsure","insufficient data","unknown","\\?","un",
    "need to confirm","breeding","400","410"),
  collapse = "|"), all_data9$prov_type),"U",all_data9$prov_type)
  # cultivated (H)
all_data9$prov_type <- ifelse(grepl(paste(
  c("cultiva","garden","^c$","^g$","^g ","^h$","horticult","landrace","clone",
    "300","500"),
  collapse = "|"), all_data9$prov_type),"H",all_data9$prov_type)
# check one last time
sort(unique(all_data9$prov_type))
  # not given (NG) ; everything else
all_data9$prov_type <- ifelse(all_data9$prov_type!= "W" &
  all_data9$prov_type != "Z" & all_data9$prov_type != "H" &
  all_data9$prov_type != "N" & all_data9$prov_type != "U",
  "NG",all_data9$prov_type)
all_data9$prov_type[which(is.na(all_data9$prov_type))] <- "NG"

# check results
table(all_data9$prov_type)
#    H    N   NG    U    W    Z
# 4731   25 2994 1736 4141  201

##
## B) Number of Individuals
##

sort(unique(all_data9$num_indiv))
  ## IF NEEDED: replace unwanted characters
  all_data9$num_indiv <- mgsub(all_data9$num_indiv,
    c(" \\(all dead\\)"," \\(removed 2017\\)","\\?",","," Plants"," or 5",
      " pieces","ca\\. ","ca\\."), c(""), fixed=F)
  sort(unique(all_data9$num_indiv))
  # change type to numeric and replace NA with 1
all_data9$num_indiv <- as.numeric(all_data9$num_indiv)
all_data9$num_indiv[which(is.na(all_data9$num_indiv))] <- 1

# check results
sort(unique(all_data9$num_indiv))
nrow(all_data9) #13828

# remove records with no individuals; save as separate file
no_indiv <- all_data9[which(all_data9$num_indiv == 0),]
nrow(no_indiv) #936
write.csv(no_indiv, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_Dead_", Sys.Date(), ".csv")),row.names = F)
  # save to in situ data folder as well
  if(!dir.exists(file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ")))
    dir.create(file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ"),
    recursive=T)
  write.csv(no_indiv, file.path(main_dir,"occurrence_points","raw_occurrence_data","Ex-situ",
    paste0("ExSitu_Dead_", Sys.Date(), ".csv")),row.names = F)
all_data9 <- all_data9[which(all_data9$num_indiv > 0),]
nrow(all_data9) #12892

##
## c) Combine individuals (same institution and accession number)
##      so that everything is (hopefully) at the accession-level
##

# preserve original acc_num before removing indiviual signifiers
all_data9$orig_acc_num <- all_data9$acc_num

# combine duplicates (same acc num)
all_data9 <- all_data9 %>%
  dplyr::group_by(inst_short,acc_num,taxon_name_acc) %>%
  dplyr::mutate(num_indiv = sum(as.numeric(num_indiv)),
         germ_type = paste(unique(germ_type),collapse="; "),
         garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(inst_short,acc_num,taxon_name_acc,.keep_all=T)
nrow(all_data9) #10684

# can look at what will be removed in the acc_num;
#   these patterns seem to work for all
all_data9[which(grepl("\\*",all_data9$acc_num)),]$acc_num
all_data9[which(grepl("/[1-9]$",all_data9$acc_num)),]$acc_num
  # this doesn't seem to work well with this dataset but may work for others:
#all_data9[which(grepl("_",all_data9$acc_num)),]$acc_num

# remove individual-specific identifiers (to combine dup accessions)
all_data9 <- all_data9 %>%
  separate("acc_num","acc_num",
    sep="\\*|/[1-9]$",remove=F) %>%
  dplyr::group_by(inst_short,acc_num,taxon_name_acc) %>%
  dplyr::mutate(num_indiv = sum(as.numeric(num_indiv)),
         germ_type = paste(unique(germ_type),collapse="; "),
         garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(inst_short,acc_num,taxon_name_acc,.keep_all=T)
nrow(all_data9) #9888

# create subset of records with acc_num longer than 9 characters
#   (these are usually the ones with plant identifiers; some are missed
#    but this gets most of them)
check_accnum <- all_data9[which(nchar(all_data9$acc_num)>9),]
nrow(check_accnum)
no_check_accnum <- setdiff(all_data9,check_accnum)
nrow(no_check_accnum)

# can look at what will be removed in the acc_num
sort(check_accnum[which(grepl("[A-F]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("-[1-9]$",check_accnum$acc_num)),]$acc_num)
  # these don't seem to work well with this dataset but may work for others:
#sort(check_accnum[which(grepl("\\.[0-9][0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
#sort(check_accnum[which(grepl("\\.[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
#sort(check_accnum[which(grepl("/[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
#sort(check_accnum[which(grepl("-[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)

# can check which records contain specific elements from above patterns:
#as.data.frame(all_data9[which(grepl("2005.0005.",all_data9$acc_num)),])

# remove individual-specific identifiers (to combine dup accessions)
check_accnum <- check_accnum %>%
  separate("acc_num","acc_num",
    sep="[A-F]$|-[1-9]$",
    remove=F) %>%
  dplyr::group_by(inst_short,acc_num,taxon_name_acc) %>%
  dplyr::mutate(num_indiv = sum(as.numeric(num_indiv)),
         germ_type = paste(unique(germ_type),collapse="; "),
         garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(inst_short,acc_num,taxon_name_acc,.keep_all=T)
nrow(check_accnum) #1344

all_data10 <- full_join(check_accnum,no_check_accnum)
nrow(all_data10) #9817

##
## ** ADD Unique ID Column
##

# add unique ID column
  # create UID with institution name, acc num, provenance type, and taxon name
  # also remove duplicates based on new UID and sum individuals
  # now create UID and remove dups
need_id <- all_data10[which(is.na(all_data10$UID)),]
dont_need_id <- all_data10[which(!is.na(all_data10$UID)),]
nms <- names(need_id)
nrow(need_id)
need_id <- need_id %>%
  dplyr::arrange(orig_lat,locality) %>%
  dplyr::mutate(UID = paste(inst_short,acc_num,prov_type,taxon_name_acc,sep="~")) %>%
  dplyr::group_by(UID) %>%
  dplyr::mutate(num_indiv = sum(as.numeric(num_indiv))) %>%
  dplyr::distinct(UID,.keep_all=T) %>%
  dplyr::ungroup() %>%
  select(c("UID",all_of(nms)))
all_data10 <- rbind(need_id,dont_need_id)
nrow(all_data10) #9817

##
## C) Latitude and Longitude
##

# preserve original lat and long columns
all_data10$lat_dd <- all_data10$orig_lat
all_data10$long_dd <- all_data10$orig_long

# fix some common lat/long character issues
all_data10$lat_dd <- mgsub(all_data10$lat_dd,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
all_data10$long_dd <- mgsub(all_data10$long_dd,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
# replace comma with decimal (european notation)
all_data10$lat_dd <- mgsub(all_data10$lat_dd, c(","), ".")
all_data10$long_dd <- mgsub(all_data10$long_dd, c(","), ".")
# separate values if lat and long both ended up in the lat column
all_data10[which(grepl("\\. ",all_data10$lat_dd)),] <-
  separate(all_data10[which(grepl("\\. ",all_data10$lat_dd)),], col = lat_dd,
           into = c("lat_dd","long_dd"), sep = "\\. ", remove = FALSE)

# replace unwanted characters
  ## latitude
  # replace random unnecessary characters
all_data10$lat_dd <- mgsub(all_data10$lat_dd,
  c("N","\\","/","M","A",": ","E","AZ","R","d","a"," \\."," W")," ")
    # remove leading zero
all_data10$lat_dd[which(grepl("^ *[0][1-9]+",all_data10$lat_dd))] <- gsub(
  "^ *[0]","",all_data10$lat_dd[which(grepl("^ *[0][1-9]+",all_data10$lat_dd))])
all_data10$lat_dd[which(grepl("^S *[0][1-9]+",all_data10$lat_dd))] <- gsub(
  "^S *[0]","-",all_data10$lat_dd[which(grepl("^S *[0][1-9]+",all_data10$lat_dd))])
    # add negative sign if south and remove "S"
all_data10$lat_dd[grep("S",all_data10$lat_dd,ignore.case=T)] <-
  paste("-",all_data10$lat_dd[grep("S",all_data10$lat_dd,ignore.case=T)],sep="")
all_data10$lat_dd <- gsub("S","",all_data10$lat_dd)
all_data10$lat_dd <- gsub("--","-",all_data10$lat_dd)
    # remove double spaces or leading/trailing whitespace
all_data10$lat_dd <- str_squish(all_data10$lat_dd)
sort(unique(all_data10$lat_dd))
  # can check source of specific values that aren't formatted correctly
#all_data10[which(all_data10$lat_dd == "422538"),]
  ## longitude
all_data10$long_dd <- replace_non_ascii(all_data10$long_dd,
  replacement=" ", remove.nonconverted=T)
all_data10$long_dd <- mgsub(all_data10$long_dd,
  c("E","\\","/","NR","d","A","a"," .","o","O")," ")
all_data10$long_dd[which(grepl("^ *[0][1-9]+",all_data10$long_dd))] <- gsub(
  "^ *[0]","",all_data10$long_dd[which(grepl("^ *[0][1-9]+",all_data10$long_dd))])
all_data10$long_dd[which(grepl("^W *[0][1-9]+",all_data10$long_dd))] <- gsub(
  "^W *[0]","-",all_data10$long_dd[which(grepl("^W *[0][1-9]+",
    all_data10$long_dd))])
all_data10$long_dd[grep("W",all_data10$long_dd,ignore.case=T)] <-
  paste("-",all_data10$long_dd[grep("W",all_data10$long_dd,ignore.case=T)],sep="")
all_data10$long_dd <- gsub("W","",all_data10$long_dd)
all_data10$long_dd <- mgsub(all_data10$long_dd,c("--","- "),"-")
all_data10$long_dd <- str_squish(all_data10$long_dd)
sort(unique(all_data10$long_dd))

# convert decimal-minutes-seconds (dms) to decimal degrees (dd)
#   [d, m, and s must be in the same cell, with 1 space between each value]
#   format = ## ## ## (DMS) OR ## ##.### (DM)
  # mark rows that need to be converted
convert <- all_data10[which(grepl(" ",all_data10$lat_dd) |
  grepl(" ",all_data10$long_dd)),]
  nrow(convert) #134
unique(convert$lat_dd)
good <- anti_join(all_data10, convert)
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
all_data10 <- rbind(good,convert); nrow(all_data10)

# check validity of lat and long
all_data10$lat_dd <- as.numeric(all_data10$lat_dd)
  #sort(unique(all_data10$lat_dd))
all_data10$long_dd <- as.numeric(all_data10$long_dd)
  #sort(unique(all_data10$long_dd))
  # if coords are both 0, set to NA
zero <- which(all_data10$lat_dd == 0 & all_data10$long_dd == 0)
all_data10$lat_dd[zero] <- NA; all_data10$long_dd[zero] <- NA
  # flag non-numeric and not available coordinates and lat > 90, lat < -90,
  # lon > 180, and lon < -180
coord_test <- cc_val(all_data10, lon = "long_dd",lat = "lat_dd",
  value = "flagged", verbose = TRUE) #Flagged 8522 records.
  # try switching lat and long for invalid points and check validity again
all_data10[!coord_test,c("lat_dd","long_dd")] <-
  all_data10[!coord_test,c("long_dd","lat_dd")]
coord_test <- cc_val(all_data10,lon = "long_dd",lat = "lat_dd",
  value = "flagged",verbose = TRUE) #Flagged 8522 records.

# check longitude values that are positive
pos_long <- all_data10[which(all_data10$long_dd > 0),]
pos_long <- pos_long %>% select(inst_short,country,lat_dd,long_dd)
unique(pos_long)
# make long value negative if in North America
setDT(all_data10)[long_dd > 0 &
                    (country == "US" | country == "United States of America" |
                     country == "United States" | country == "U.S.A." |
                     country == "Canada" | country == "USA"),
                  long_dd := as.numeric(paste0("-",as.character(long_dd)))]

# make coords NA if they are still flagged
coord_test <- cc_val(all_data10,lon = "long_dd",lat = "lat_dd",
  value = "flagged",verbose = TRUE) #Flagged 7951 records.
all_data10[!coord_test,"lat_dd"] <- NA
all_data10[!coord_test,"long_dd"] <- NA

# check if geolocated points are in water and mark
world_polygons <- ne_countries(type = 'countries', scale = 'large')
geo_pts <- all_data10 %>% filter(!is.na(lat_dd) & !is.na(long_dd))
in_water <- geo_pts[is.na(map.where(world_polygons,
  geo_pts$long_dd,geo_pts$lat_dd)),]
nrow(in_water)
all_data10$flag <- ""
all_data10[which(all_data10$UID %in% in_water$UID),]$flag <-
  "Given lat-long is in water"
table(all_data10$flag)
#all_data10[which(all_data10$UID %in% in_water$UID),]$lat_dd <- NA
#all_data10[which(all_data10$UID %in% in_water$UID),]$long_dd <- NA

# mark lat-long for records with same inst lat-long and wild lat-long
all_data10$lat_round <- round(all_data10$lat_dd,digits=1)
all_data10$long_round <- round(all_data10$long_dd,digits=1)
all_data10$inst_lat_round <- round(all_data10$inst_lat,digits=1)
all_data10$inst_long_round <- round(all_data10$inst_long,digits=1)
garden_latlong <- all_data10 %>% filter(lat_round == inst_lat_round &
  long_round == inst_long_round & prov_type != "N")
unique(garden_latlong$inst_short)
nrow(garden_latlong)
all_data10[which(all_data10$UID %in% garden_latlong$UID),]$flag <-
  "Given lat-long is at institution, use only if native to grounds"
#all_data10[all_data10$UID %in% garden_latlong$UID,]$lat_dd <- NA
#all_data10[all_data10$UID %in% garden_latlong$UID,]$long_dd <- NA
table(all_data10$flag) #grounds = 107; water = 102

# add country-level information to check if lat-long in right spot
# create SpatialPointsDataFrame
proj4string4poly <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
geo_pts_spatial <- SpatialPointsDataFrame(geo_pts[,c("long_dd",
  "lat_dd")], geo_pts, proj4string = CRS(proj4string4poly))
# add country polygon data to each point based on lat-long location
geo_pts <- point.in.poly(geo_pts_spatial, world_polygons, sp=TRUE)@data
# try switching lat and long for points in Antarctica
geo_pts[which(geo_pts$country.name == "Antarctica"),c("lat_dd","long_dd")]<-
  geo_pts[which(geo_pts$country.name == "Antarctica"),c("long_dd","lat_dd")]
# round 2: add country-level information to check if lat-long in right spot
geo_pts <- geo_pts %>% select(UID:long_dd)
geo_pts_spatial <- SpatialPointsDataFrame(geo_pts[,c("long_dd",
  "lat_dd")], geo_pts, proj4string = CRS(proj4string4poly))
geo_pts <- point.in.poly(geo_pts_spatial, world_polygons, sp=TRUE)@data
geo_pts <- geo_pts %>% select(UID,admin) %>%
  dplyr::rename(latlong_country = admin)
all_data10 <- full_join(all_data10,geo_pts)

# add gps_det (gps determination) column
#all_data10$gps_det <- NA
all_data10$gps_det[which(all_data10$prov_type == "H")] <- "N/A (horticultural)"
all_data10$gps_det[which(!is.na(all_data10$lat_dd) &
  !is.na(all_data10$long_dd))] <- "Given"
all_data10$gps_det[which(all_data10$gps_det == "")] <- NA
table(all_data10$gps_det)
# Given     N/A (horticultural)
# 1295      3250

# where prov_type is "N/A (horticultural)" but lat-long is given, change to "H?"
  # create new prov type column
all_data10$prov_type[which(all_data10$gps_det == "Given" &
  all_data10$prov_type == "H")] <- "H?"
table(all_data10$prov_type)
#    H   H?    N   NG    U    W    Z
# 3250  123   25 1963 1456 2845  155

##
## D) Collection year
##

# not cleaning up this column right now
sort(unique(all_data10$coll_year))
## IF NEEDED: replace non-year words/characters
#all_data10$coll_year <- mgsub(all_data10$coll_year,
#  c("([0-9]+);","about ","ca.","Unknown","original","Estate","estate"),"")
#all_data10$coll_year[which(all_data10$coll_year == "")] <- NA

# remove extra elements so its just year
#all_data10$coll_year <- gsub(";[1-2][0-9][0-9][0-9]","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9][0-9]/[0-9][0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9]/[0-9][0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9][0-9]/[0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9]/[0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[1-9] [A-Z][a-z][a-z] ","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[A-Z][a-z][a-z] ","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[1-9]-[A-Z][a-z][a-z]-","",all_data10$coll_year)

# make column numeric
#all_data10$coll_year <- as.numeric(all_data10$coll_year)
## IF NEEDED: add first two numbers in year
#  # assume 2000s if values is less than 21
#all_data10$coll_year[which(all_data10$coll_year < 10)] <-
#  paste0("200",as.character(all_data10$coll_year[which(all_data10$coll_year < 10)]))
#all_data10$coll_year <- as.numeric(all_data10$coll_year)
#all_data10$coll_year[which(all_data10$coll_year < 21)] <-
#  paste0("20",as.character(all_data10$coll_year[which(all_data10$coll_year < 21)]))
#all_data10$coll_year <- as.numeric(all_data10$coll_year)
#  # assume 1900s if values is greater than or equal to 21
#all_data10$coll_year[which(all_data10$coll_year < 100)] <-
#  paste0("19",as.character(all_data10$coll_year[which(all_data10$coll_year < 100)]))
#all_data10$coll_year <- as.numeric(all_data10$coll_year)
#sort(unique(all_data10$coll_year))

##
## E) Lineage number
##

# remove lin_num when same as acc_num
all_data10[which(all_data10$acc_num == all_data10$lin_num),]$lin_num <- NA

##
## F) Locality
##

# create all_locality column
all_data10$latitude <- round(all_data10$lat_dd,digits=3)
all_data10$longitude <- round(all_data10$long_dd,digits=3)
all_data10 <- unite(all_data10, "all_locality",
  c(locality,municipality,county,state,country,orig_source,
    lin_num,coll_num,coll_name,coll_year,
    latitude,longitude,notes),sep = " | ",remove = F)
# remove NA in concatenated locality column
all_data10$all_locality <- gsub("NA","",all_data10$all_locality)
# if no locality info at all, make it NA
all_data10$all_locality[which(all_data10$all_locality ==
  " |  |  |  |  |  |  |  |  |  |  |  | ")] <- NA

##
## G) Institution type
##

# add inst_type for gene bank data
all_data10$inst_type[which(is.na(all_data10$inst_type))] <- "Gene/Seed Bank"
table(all_data10$inst_type)
# Botanic Garden    Botanic Garden; Gene/Seed Bank     Gene/Seed Bank
# 6195              159                                3463
table(all_data10$data_source)
# ex_situ_BG_survey  FAO-WIEWS   Genesys
# 7551               2060        206

################################################################################
# 7. Summary statistics
################################################################################

# keep only necessary columns
data_sel <- all_data10 %>%
  select(
    # key data
    UID,inst_short,inst_short2,taxon_name_acc,species_name_acc,
    prov_type,num_indiv,acc_num,lat_dd,long_dd,flag,gps_det,
    # locality
    all_locality,locality,municipality,county,state,country,assoc_sp,
    latlong_country,orig_lat,orig_long,uncertainty,orig_prov_type,
    # source
    orig_source,lin_num,coll_num,coll_name,coll_year,coll_date,
    # material info
    germ_type,garden_loc,rec_as,rec_year,rec_date,
    # other metadata
    notes,filename,submission_year,data_source,taxon_name_status,orig_acc_num,
    fruit_nut,iucnredlist_category,natureserve_rank,
    # taxon name details
    taxon_name,genus,species,infra_rank,infra_name,cultivar,hybrid,
    taxon_full_name_orig,taxon_full_name_concat,taxon_verif,
    # institution metadata
    inst_id,inst_type,inst_lat,inst_long,inst_country,date_data_received)

# write file
write.csv(data_sel, file.path(main_dir,"ex-situ_data","OUTPUTS_FROM_R",
  paste0("ExSitu_Compiled_", Sys.Date(), ".csv")),row.names = F)

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
  dplyr::group_by(taxon_name_acc) %>%
  dplyr::summarize(
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
  dplyr::filter(is.na(lat_dd) & !is.na(all_locality) & prov_type != "H")
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


### !!!! ENDED HERE !!!!









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
