################################################################################

## 2-1_compile_exsitu_data_CWRcopy.R

### Author: Emily Beckman Bruns
### Funding: Cooperative agreement between United States Botanic Garden and
#    San Diego Botanic Garden; subcontract with The Morton Arboretum.
#    With support from Botanic Gardens Conservation International U.S.

### Creation date: 30 March 2022
### Last updated: 23 May 2022

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
  #     | inst_short                               | inst_id  | inst_type      | inst_lat  | inst_long | inst_country |
  #     |------------------------------------------|----------|----------------|-----------|-----------|--------------|
  #     |nickname for institution (see template    |instition |institution type|institution|institution|institution   |
  #     |standardizing_accessions_data_fields.xlsx)|identifier|(Botanic Garden,|latitude   |longitude  |country       |
  #     |                                          |          | Gene/Seed Bank)|           |           |              |
  #
  # 3. Target taxa list (target_taxa_with_syn.csv); required columns include:
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
my.packages <- c('plyr','tidyverse', 'data.table', #'anchors',
                 'textclean',
                 'measurements', #'naniar',
                 'CoordinateCleaner','rnaturalearth',
                 'rnaturalearthdata','maps','raster','spatialEco','leaflet'
                )
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

# Google Drive folder with accession data
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Ex situ - G - records"
# Google Drive folder with polygon data
poly_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"

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
##    this function also adds columns for 1) the file name [equivalent to the
##    "inst_short" institution nickname] 2) a sumbission year, 3) an accession
##    number if one isn't given
#raw_2020 <- read.exsitu.csv(file.path(main_dir,"inputs",
#  "exsitu_standard_column_names","data_2020"), "2020")
raw_2022 <- read.exsitu.csv(file.path(main_dir,"inputs",
  "exsitu_standard_column_names"), "2022")
# stack all data
#to_stack <- list(raw_2020,raw_2022)
#all_data_raw <- Reduce(rbind.fill,to_stack)

# create new version before big changes, so can easily go back to original
all_data <- raw_2022

# replace non-ascii characters
  # first fix some common lat/long character issues
all_data$orig_lat <- mgsub(all_data$orig_lat,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
all_data$orig_long <- mgsub(all_data$orig_long,
  c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´","*","À","?","`")," ")
  # replace all non-ascii
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

################################################################################
# 2. Compile Genesys data
################################################################################

# field metadata: https://www.genesys-pgr.org/documentation/basics

# load in data
  # collection site description
genPath <- paste0(main_dir,"/inputs/genesys-accessions/coll.csv")
gen_col <- data.table::fread(file = genPath, header = TRUE)
  str(gen_col)
  nrow(gen_col) ; length(unique(gen_col$genesysId)) #49998
  # metadata about record
genPath <- paste0(main_dir,"/inputs/genesys-accessions/core.csv")
gen_core <- data.table::fread(genPath, header = TRUE)
  str(gen_core)
  nrow(gen_core) ; length(unique(gen_core$genesysId)) #78432
  # spatial info
  #   for some reason the headers dont read in correctly with fread;
  #   will use read.csv just to get the headers and add them to the fread df
genPath <- paste0(main_dir,"/inputs/genesys-accessions/geo.csv")
gen_geo <- data.table::fread(genPath, header = TRUE)
  str(gen_geo)
  colnames(gen_geo)
gen_geo_head <- read.csv(genPath)
  colnames(gen_geo_head)
gen_geo <- gen_geo[,1:7] #remove last column (nothing in it)
colnames(gen_geo) <- colnames(gen_geo_head)
  str(gen_geo)
  nrow(gen_geo) ; length(unique(gen_geo$genesysId)) #78431
  # this has mult. entires for each ID
  #   will concatenate them so one entry per ID
  # Actually, I dont think this has useful info for us; not using
#genPath <- paste0(main_dir,"/inputs/genesys-accessions/names.csv")
#gen_names <- data.table::fread(genPath, header = TRUE)
#  str(gen_names)
#  nrow(gen_names) ; length(unique(gen_names$genesysId))
#gen_names <- gen_names %>%
#  group_by(genesysId) %>%
#  summarize(
#    instCode = paste(instCode, collapse = '; '),
#    name = paste(name, collapse = '; '),
#    aliasType = paste(aliasType, collapse = '; '),
#    version = paste(version, collapse = '; ')) %>%
#  ungroup()
#  str(gen_names)
#  nrow(gen_names) ; length(unique(gen_names$genesysId))

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
  # add Genesys to ID
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
  dplyr::select(UID,taxon_full_name,inst_short2,coll_num,coll_date,locality,notes,
    inst_short,acc_num,infra_rank,country,rec_date,orig_lat,orig_long,gps_det,
    germ_type,orig_source,prov_type,genus,species,uncertainty)

# filter by target taxa & add taxon data
## We do this later!
#genesys_sel <- left_join(genesys_sel,taxon_list)
#genesys_target <- genesys_sel %>% filter(!is.na(taxon_name_acc))
#nrow(genesys_target) #5233

# US institutions covered by GRIN? We will remove these
#   list from CWR Gap Analysis 2020 (PNAS)...
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
wiewsPath <- paste0(main_dir,"/inputs/Wiews_Exsitu.csv")
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
  dplyr::select(inst_short,acc_num,taxon_full_name,genus,species,rec_date,country,
         prov_type,inst_short2,orig_lat,orig_long,orig_source,germ_type)

# !!! JUST FOR NOW; need to check exact to filter out !!!
#   US institutions (covered by GRIN?)
# looks like they are all removed by removing Genesys anyways?
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
all_data2 <- all_data2 %>% separate("taxon_full_name","genus_temp",sep=" ",remove=F)
all_data2[which(is.na(all_data2$genus)),]$genus <- all_data2[which(is.na(all_data2$genus)),]$genus_temp
# standardize capitalization
all_data2$genus <- str_to_title(all_data2$genus)

### MAKE SURE NO GENUS MISSPELLINGS OR ABBREVIATIONS ###
sort(unique(all_data2$genus))
#all_data2$genus <- mgsub(all_data2$genus,
#  c("^Q$","Querucs","Cyclobalanopsis"),"Quercus",fixed=F)

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "inputs", "target_taxa_with_syn.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
str(taxon_list) #266

# remove rows not in target genus/genera
target_genera <- unique(taxon_list$genus)
all_data3 <- all_data2 %>% filter(genus %in% target_genera)
nrow(all_data2); nrow(all_data3) #170637 ; 157687

### CHECK OUT THE HYBRID COLUMN ###
# standardize a bit
sort(unique(all_data3$hybrid))
#all_data3 <- replace.value(all_data3, "hybrid", c("species","unknown"), NA)
all_data3$hybrid <- mgsub(all_data3$hybrid,
  c(" A ","^A ","^A$"," X ","^X ","^X$"," _ ","^_ ","^_$"),
  " x ", fixed=F)
all_data3$hybrid <- str_squish(all_data3$hybrid)
# make sure everything has an " x " in it somewhere (important later)
all_data3$hybrid <- mgsub(all_data3$hybrid,
  c("M. ioensisxpumila","M. baccataxspectabilis"),
  c("M. ioensis x pumila","M. baccata x spectabilis"))

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
# make sure author is separated in taxon_full_name
all_data3$taxon_full_name <- gsub("\\("," (",all_data3$taxon_full_name)
all_data3$taxon_full_name <- str_squish(all_data3$taxon_full_name)

# write copy of all data
  # select columns
all_data_export1 <- all_data3 %>%
  dplyr::select(inst_short,inst_short2,taxon_full_name,genus,species,infra_rank,
    infra_name,hybrid,cultivar,taxon_full_name_orig,prov_type,orig_lat,orig_long,
    country,state,municipality,locality,assoc_sp,acc_num,lin_num,orig_source,
    rec_as,rec_year,rec_date,num_indiv,germ_type,garden_loc,coll_num,coll_name,
    coll_year,coll_date,taxon_verif,notes,filename,submission_year,data_source)
write.csv(all_data_export1, file.path(main_dir,"outputs",
  paste0("ExSitu_AllDataRaw_GenusFilterOnly_", Sys.Date(), ".csv")),row.names = F)

# summary of genera for each institution
gen_summary <- all_data3 %>%
  dplyr::arrange(genus) %>%
  dplyr::rename(genera = genus) %>%
  dplyr::group_by(inst_short) %>%
  dplyr::mutate(
    genera = paste(unique(genera), collapse = '; ')) %>%
  dplyr::ungroup() %>%
  dplyr::select(inst_short,genera) %>%
  dplyr::distinct(inst_short,.keep_all=T)
head(as.data.frame(gen_summary))
write.csv(gen_summary, file.path(main_dir,"outputs",
  paste0("Genera_Reported_by_Institutions_", Sys.Date(), ".csv")),row.names = F)

################################################################################
# 3. Further standardize taxon name, then keep data for target taxa only
#     (removes hybrids and cultivars without specific epithet)
################################################################################

## FINISH STANDARDIZING TAXON NAMES

# add space after periods in taxon_full_name
all_data3$taxon_full_name <- gsub(".",". ",all_data3$taxon_full_name,fixed=T)
# replace unwanted characters in taxon_full_name
#all_data3$taxon_full_name <- mgsub(all_data3$taxon_full_name,
#  c("(",")",";","[","]",",","^","#")," ")
all_data3$taxon_full_name <- str_squish(all_data3$taxon_full_name)

## REMOVE HYBRIDS

# remove hybrids based on " x " in taxon_full_name_concat and/or taxon_full_name
#all_data4 <- all_data3 %>% filter(!grepl(" x ",taxon_full_name_concat))
#nrow(all_data3); nrow(all_data4) #45057 ; 43363
#all_data4 <- all_data4 %>% filter(!grepl(" x ",taxon_full_name))
#nrow(all_data4) #41936
# see hybrids removed:
#sort(unique(anti_join(all_data3,all_data4)$taxon_full_name))

## CREATE NEW SEPARATED TAXON NAME COLUMNS

# fix anything before removing:
all_data3$taxon_full_name <- mgsub(all_data3$taxon_full_name,
  c("Corylus Acolurnoides","Prunus Ayedoensis","Pyrus (Malus) fusca"),
  c("Corylus x colurnoides","Prunus x yedoensis","Malus fusca"))
## first change species hybrids temporarily: remove space so stays together
all_data3$taxon_full_name <- gsub(" x "," +",all_data3$taxon_full_name)
# separate out taxon full name and trim whitespace again
all_data3 <- all_data3 %>% separate("taxon_full_name",
  c("genus_new","species_new","extra1","extra2",
    "extra3","extra4","extra5","extra6","extra7"),sep=" ",extra="warn",
    remove=F,fill="right")
all_data3 <- as.data.frame(lapply(all_data3,str_squish),stringsAsFactors=F)
# replace genus_new with genus, since we fixed that up in the previous section
all_data3$genus_new <- all_data3$genus

## REMOVE RECORDS WITHOUT SPECIES NAME

# remove records with no/non-standard specific epithet
#   (by looking in species name column)
all_data4 <- all_data3 %>%
  dplyr::filter(!grepl("\"",species_new) &
                !grepl("\'",species_new) &
                !grepl("\\[",species_new) &
                !grepl("\\(",species_new) &
                !grepl("\\.",species_new) &
                !grepl("[A-Z]",species_new) &
                !grepl("[0-9]",species_new) &
                !grepl("\\?",species_new) &
                !is.na(species_new))
nrow(all_data4) #
# see records removed; can add anything you want to fix to the
#   "fix anything before removing" section above:
sort(unique(anti_join(all_data3,all_data4)$taxon_full_name))

## FIND INFRATAXA

## look for infrataxa key words
# make data in all "extra" columns lower case
sp_col <- grep("^species_new$", colnames(all_data4))
all_data4[,sp_col:(sp_col+5)] <- as.data.frame(sapply(
  all_data4[,sp_col:(sp_col+5)], tolower), stringsAsFactors=F)
# create matrix of all "extra" species name columns, to search for
#   infraspecific key words
search.col <- matrix(cbind(all_data4$extra1,all_data4$extra2,all_data4$extra3,
  all_data4$extra4,all_data4$extra5,all_data4$extra6,all_data4$extra7),
  nrow=nrow(all_data4))
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
all_data4$infra_rank_new <- NA
all_data4$infra_rank_new[matches_i] <- all_data4[matches_i]
#unique(all_data4$infra_rank_new) # check results

# create new infra_name column and fill with next column over from "extra"
#   contents that matched infraspecific key word
all_data4$infra_name_new <- NA
matches_i[,2] <- matches_i[,2]+1
all_data4$infra_name_new[matches_i] <- all_data4[matches_i]
#sort(unique(all_data4$infra_name_new))

# standardize infraspecific rank names
all_data4$infra_rank_new <- replace(all_data4$infra_rank_new,
  grep("^v$|^v.$|^var$|^variety$|^va$",all_data4$infra_rank_new), "var.")
all_data4$infra_rank_new <- replace(all_data4$infra_rank_new,
  grep("^subspecies$|^subsp$|^ssp$|^ssp.$|^subs.$|^spp.$|^sub.$",
  all_data4$infra_rank_new), "subsp.")
all_data4$infra_rank_new <- replace(all_data4$infra_rank_new,
 grep("^forma$|^form$|^fma$|^fo$|^fo.$|^f$",all_data4$infra_rank_new), "f.")
unique(all_data4$infra_rank_new)

## CREATE FINAL TAXON FULL NAME FOR FILTERING

# create new taxon full name column
all_data4$taxon_full_name <- NA
  # select rows with infraspecific name and concatenate
yes_infra <- which(!is.na(all_data4$infra_rank_new) &
  !is.na(all_data4$infra_name_new))
all_data4$taxon_full_name[yes_infra] <- paste(all_data4$genus_new[yes_infra],
  all_data4$species_new[yes_infra], all_data4$infra_rank_new[yes_infra],
  all_data4$infra_name_new[yes_infra],sep=" ")
  # select rows without infraspecific name and concatenate
all_data4$taxon_full_name[-yes_infra] <- paste(all_data4$genus_new[-yes_infra],
  all_data4$species_new[-yes_infra],sep=" ")
# check out results
sort(unique(all_data4$taxon_full_name))
# switch hybrid symbol back to " x "
#all_data4$taxon_full_name <- gsub(" \\+"," x ",all_data4$taxon_full_name)
#sort(unique(all_data4$taxon_full_name))

## FILTER OUT NON-TARGET TAXA

# rename some taxon name columns to preserve originals
all_data4 <- all_data4 %>%
  dplyr::rename(taxon_name = taxon_full_name,
                genus_orig = genus,
                species_orig = species,
                infra_rank_orig = infra_rank,
                infra_name_orig = infra_name)
all_data4 <- all_data4 %>%
  dplyr::rename(genus = genus_new,
                species = species_new,
                infra_rank = infra_rank_new,
                infra_name = infra_name_new)

# join dataset to taxa list
  # join by taxon name
all_data5 <- left_join(all_data4,taxon_list)
  # if no taxon match, join again just by species name
need_match <- all_data5[which(is.na(all_data5$taxon_name_status)),]
  nrow(need_match) #138992 ; 139837
    # remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data5)-ncol(taxon_list)+2)]
    # remove taxon_name col from taxon data so it doesn't match
taxon_list_sp <- taxon_list %>% select(-taxon_name)
  # create genus_species column
need_match$genus_species <- paste(need_match$genus,need_match$species)
    # new join
need_match <- left_join(need_match,taxon_list_sp)
  # bind together new matches and previously matched
matched <- all_data5[which(!is.na(all_data5$taxon_name_status)),]
all_data6 <- rbind(matched,need_match)
  table(all_data6$taxon_name_status) # Accepted: 11836 | Synonym: 1821
  # see how many rows have taxon name match
nrow(all_data6[which(!is.na(all_data6$taxon_name_status)),]) #13532 ; 13657; 13828

### CHECK UNMATCHED SPECIES, TO ADD TO SYNONYM LIST AS NECESSARY ###
check <- all_data6 %>%
  dplyr::filter(is.na(taxon_name_status)) %>%
  dplyr::select(taxon_name,taxon_full_name_orig,taxon_full_name_concat)
head(check)
check <- data.frame(
  taxon_name = sort(unique(check$taxon_name)),
  taxon_name_acc = rep(NA),
  fruit_nut = rep(NA))
# write file for checking, as desired
write.csv(check, file.path(main_dir,"outputs",
  "ExSitu_UnmatchedSpecies.csv"),row.names = F)
  ## fill taxon_name_acc column with correct accepted name, where applicable
  ## save as ExSitu_UnmatchedSpecies-filled.csv

# [OPTIONAL] read in additional corrected names (saved in previous step)
add_names <- read.csv(file.path(main_dir,"outputs","ExSitu_UnmatchedSpecies-filled.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
add_names <- add_names %>%
  dplyr::filter(!is.na(taxon_name_acc)) %>%
  separate("taxon_name",c("genus","species"),sep=" ",remove=F) %>%
  separate("taxon_name_acc",c("genus_acc","species_acc"),sep=" ",remove=F)
add_names$genus_species <- paste(add_names$genus,add_names$species)
add_names$species_name_acc <- paste(add_names$genus_acc,add_names$species_acc)
add_names$taxon_name_status <- "Synonym"
add_names <- add_names %>%
  dplyr::select(taxon_name_acc,species_name_acc,taxon_name,genus,genus_species,taxon_name_status,fruit_nut)
add_names$rl_category <- "Not Evaluated"
add_names$rl_year <- NA
add_names
taxon_list <- full_join(taxon_list,add_names)
## !!!!
## !!!!!! now run the "join dataset to taxa list" section again !!!!!
## !!!!

# keep only matched names
all_data7 <- all_data6 %>% dplyr::filter(!is.na(taxon_name_status))
nrow(all_data7) #10533 ; 13828

################################################################################
# 4. Standardize important columns
################################################################################

# keep only necessary columns
all_data8 <- all_data7 %>% dplyr::select(
  # key data
  UID,inst_short,inst_short2,taxon_name_acc,species_name_acc,rl_category,
  rl_year,prov_type,num_indiv,acc_num,
  # locality
  orig_lat,orig_long,locality,municipality,county,state,country,assoc_sp,
  gps_det,uncertainty,
  # source
  orig_source,lin_num,coll_num,coll_name,coll_year,coll_date,
  # material info
  germ_type,garden_loc,rec_as,rec_year,rec_date,
  # other metadata
  notes,filename,submission_year,data_source,taxon_name_status,fruit_nut,
  # taxon name details
  taxon_name,genus,species,infra_rank,infra_name,cultivar,hybrid,
  taxon_full_name_orig,taxon_full_name_concat,taxon_verif)

# add institution metadata
inst_data <- read.csv(file.path(main_dir,"inputs",
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
# 4635   24 2957 1729 4111  201

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
write.csv(no_indiv, file.path(main_dir,"outputs",
  paste0("ExSitu_Dead_", Sys.Date(), ".csv")),row.names = F)
all_data9 <- all_data9[which(all_data9$num_indiv > 0),]
nrow(all_data9) #12892

##
## c) Combine duplicates (same institution and accession number)
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
nrow(all_data9) #10088

# can look at what will be removed in the acc_num;
#   these patterns seem to work for all
all_data9[which(grepl("\\*",all_data9$acc_num)),]$acc_num
#all_data9[which(grepl("_",all_data9$acc_num)),]$acc_num #this doesn't work well
all_data9[which(grepl("/[1-9]$",all_data9$acc_num)),]$acc_num

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
nrow(all_data9) #9394

# create subset of records with acc_num longer than 9 characters
#   (these are usually the ones with plant identifiers; some are missed
#    but this gets most of them)
check_accnum <- all_data9[which(nchar(all_data9$acc_num)>9),]
nrow(check_accnum) #1411
no_check_accnum <- setdiff(all_data9,check_accnum)
nrow(no_check_accnum) #7983

# can look at what will be removed in the acc_num
sort(check_accnum[which(grepl("/[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("\\.[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("\\.[0-9][0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("[A-F]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("-[1-9]$",check_accnum$acc_num)),]$acc_num)
sort(check_accnum[which(grepl("-[0-9][1-9]$",check_accnum$acc_num)),]$acc_num)
  #as.data.frame(check_accnum[which(grepl("A 1971-432",check_accnum$acc_num)),])

# remove individual-specific identifiers (to combine dup accessions)
check_accnum <- check_accnum %>%
  separate("acc_num","acc_num",
    sep="/[0-9][1-9]$|\\.[0-9][1-9]$|\\.[0-9][0-9][1-9]$|[A-F]$|-[1-9]$|-[0-9][1-9]$",
    remove=F) %>%
  dplyr::group_by(inst_short,acc_num,taxon_name_acc) %>%
  dplyr::mutate(num_indiv = sum(as.numeric(num_indiv)),
         germ_type = paste(unique(germ_type),collapse="; "),
         garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(inst_short,acc_num,taxon_name_acc,.keep_all=T)
nrow(check_accnum) #1315

all_data10 <- full_join(check_accnum,no_check_accnum)
nrow(all_data10) #9298

# look at acc_num with potential qualifiers that were not removed;
#   can fix manually if desired
#all_data10[which(grepl("\\*",all_data10$acc_num)),]$acc_num
#all_data10[which(grepl("_",all_data10$acc_num)),]$acc_num
#all_data10[which(grepl("/[1-9]$",all_data10$acc_num)),]$acc_num
all_data10[which(grepl("/[0-9][1-9]$",all_data10$acc_num)),]$acc_num
all_data10[which(grepl("\\.[0-9][1-9]$",all_data10$acc_num)),]$acc_num
all_data10[which(grepl("\\.[0-9][0-9][1-9]$",all_data10$acc_num)),]$acc_num
all_data10[which(grepl("[A-F]$",all_data10$acc_num)),]$acc_num
#  all_data10 <- all_data10 %>% separate("acc_num","acc_num",sep="[A-F]$",remove=F)
  #as.data.frame(all_data10[which(grepl("159[A-F]",all_data10$acc_num)),])
#sort(all_data10[which(grepl("-[1-9]$",all_data10$acc_num) & nchar(all_data10$acc_num)>6),]$acc_num)
#  all_data10[which(grepl("10796-",all_data10$acc_num)),]$acc_num <- "10796"
#  all_data10[which(grepl("88I54-",all_data10$acc_num)),]$acc_num <- "88I54"
all_data10[which(grepl("-[0-9][1-9]$",all_data10$acc_num)),]$acc_num

# combine duplicates one final time
#all_data10 <- all_data10 %>%
#  group_by(inst_short,acc_num,taxon_name_acc) %>%
#  mutate(num_indiv = sum(as.numeric(num_indiv)),
#         germ_type = paste(unique(germ_type),collapse="; "),
#         garden_loc = paste(unique(garden_loc),collapse="; ")) %>%
#  ungroup() %>%
#  distinct(inst_short,acc_num,taxon_name_acc,.keep_all=T)
#nrow(all_data10)

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
  dplyr::select(c("UID",all_of(nms)))
all_data10 <- rbind(need_id,dont_need_id)
nrow(all_data10) #9297

##
## C) Latitude and Longitude
##

# preserve original lat and long columns
all_data10$lat_dd <- all_data10$orig_lat
all_data10$long_dd <- all_data10$orig_long

# replace comma with decimal (european notation)
all_data10$lat_dd <- mgsub(all_data10$lat_dd, c(","), ".")
all_data10$long_dd <- mgsub(all_data10$long_dd, c(","), ".")

# replace unwanted characters
  ## latitude
  # replace random unnecessary characters
all_data10$lat_dd <- mgsub(all_data10$lat_dd,
  c("N","\\","/","M","A",": ","E","AZ","R","d","a"," .")," ")
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
#sort(unique(all_data10$lat_dd))
  # check source of specific values that aren't formatted correctly
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
#sort(unique(all_data10$long_dd))

# convert decimal-minutes-seconds (dms) to decimal degrees (dd)
#   [d, m, and s must be in the same cell, with 1 space between each value]
#   format = ## ## ## (DMS) OR ## ##.### (DM)
  # mark rows that need to be converted
convert <- all_data10[which(grepl(" ",all_data10$lat_dd) |
  grepl(" ",all_data10$long_dd)),]
  nrow(convert) #140
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
  convert <- rbind(dms,ddm,other); nrow(convert) #140
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
  value = "flagged", verbose = TRUE) #Flagged 8034 records.
  # try switching lat and long for invalid points and check validity again
all_data10[!coord_test,c("lat_dd","long_dd")] <-
  all_data10[!coord_test,c("long_dd","lat_dd")]
coord_test <- cc_val(all_data10,lon = "long_dd",lat = "lat_dd",
  value = "flagged",verbose = TRUE) #Flagged 8034 records.

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
  value = "flagged",verbose = TRUE) #Flagged 7932 records.
all_data10[!coord_test,"lat_dd"] <- NA
all_data10[!coord_test,"long_dd"] <- NA

# check if geolocated points are in water and mark
world_polygons <- ne_countries(type = 'countries', scale = 'medium')
geo_pts <- all_data10 %>% filter(!is.na(lat_dd) & !is.na(long_dd))
in_water <- geo_pts[is.na(map.where(world_polygons,
  geo_pts$long_dd,geo_pts$lat_dd)),]
nrow(in_water)
all_data10$flag <- ""
all_data10[which(all_data10$UID %in% in_water$UID),]$flag <-
  "Given lat-long is in water"
table(all_data10$flag) #110
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
table(all_data10$flag) #grounds = 89; water = 110

# add country-level information to check if lat-long in right spot
# create SpatialPointsDataFrame
proj4string4poly <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
geo_pts_spatial <- SpatialPointsDataFrame(geo_pts[,c("long_dd",
  "lat_dd")], geo_pts, proj4string = CRS(proj4string4poly))
# add country polygon data to each point based on lat-long location
load(file.path(poly_dir, "inputs", "gis_data", "admin_shapefiles.RData"))
geo_pts <- point.in.poly(geo_pts_spatial, adm0.poly, sp=TRUE)@data
# try switching lat and long for points in Antarctica
geo_pts[which(geo_pts$country.name == "Antarctica"),c("lat_dd","long_dd")]<-
  geo_pts[which(geo_pts$country.name == "Antarctica"),c("long_dd","lat_dd")]
# round 2: add country-level information to check if lat-long in right spot
geo_pts <- geo_pts %>%
  dplyr::select(-country.name,-country.iso_a2,-country.iso_a3,-country.continent)
geo_pts_spatial <- SpatialPointsDataFrame(geo_pts[,c("long_dd",
  "lat_dd")], geo_pts, proj4string = CRS(proj4string4poly))
geo_pts <- point.in.poly(geo_pts_spatial, adm0.poly, sp=TRUE)@data
geo_pts <- geo_pts %>% dplyr::select(UID,country.name) %>%
  dplyr::rename(latlong_country = country.name)
all_data10 <- full_join(all_data10,geo_pts)

# add gps_det (gps determination) column
#all_data10$gps_det <- NA
all_data10$gps_det[which(all_data10$prov_type == "H")] <- "N/A (horticultural)"
all_data10$gps_det[which(!is.na(all_data10$lat_dd) &
  !is.na(all_data10$long_dd))] <- "Given"
all_data10$gps_det[which(all_data10$gps_det == "")] <- "Unknown"
table(all_data10$gps_det)
# Given     N/A (horticultural)     Unknown
# 1263      2788                    109

# where prov_type is "N/A (horticultural)" but lat-long is given, change to "H?"
  # create new prov type column
all_data10$prov_type[which(all_data10$gps_det == "Given" &
  all_data10$prov_type == "H")] <- "H?"
table(all_data10$prov_type)
#    H   H?    N   NG    U    W    Z
# 2788  123   25 1945 1453 2808  155



######################################################
######## !!! SKIPPING THIS FOR NOW !!!################

##
## D) Collection year
##

sort(unique(all_data10$coll_year))
## IF NEEDED: replace non-year words/characters
#all_data10$coll_year <- mgsub(all_data10$coll_year,
#  c("([0-9]+);","about ","ca.","Unknown","original","Estate","estate"),"")
## IF NEEDED: replace non-year words/characters

#all_data10$coll_year[which(all_data10$coll_year == "")] <- NA

# remove extra elements so its just year
all_data10$coll_year <- gsub(";[1-2][0-9][0-9][0-9]","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9][0-9]/[0-9][0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9]/[0-9][0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9][0-9]/[0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[0-9]/[0-9]/","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[1-9] [A-Z][a-z][a-z] ","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[A-Z][a-z][a-z] ","",all_data10$coll_year)
#all_data10$coll_year <- gsub("[1-9]-[A-Z][a-z][a-z]-","",all_data10$coll_year)

# make column numeric
all_data10$coll_year <- as.numeric(all_data10$coll_year)
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
sort(unique(all_data10$coll_year))

######################################################
######################################################


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
  c(locality,municipality,county,state,country,orig_source,#notes,
    lin_num,coll_num,coll_name,coll_year,
    latitude,longitude),sep = " | ",remove = F)
# remove NA in concatenated locality column
all_data10$all_locality <- gsub("NA","",all_data10$all_locality)
# if no locality info at all, make it NA
all_data10$all_locality[which(all_data10$all_locality ==
  " |  |  |  |  |  |  |  |  |  |  | ")] <- NA

##
## G) Institution type
##

# add inst_type for gene bank data
all_data10$inst_type[which(is.na(all_data10$inst_type))] <- "Gene/Seed Bank"
table(all_data10$inst_type)
# Botanic Garden    Botanic Garden; Gene/Seed Bank     Gene/Seed Bank
# 6137              157                                3003
table(all_data10$data_source)
# ex_situ_BG_survey  FAO-WIEWS   Genesys
# 7485               1612        200

################################################################################
# 7. Summary statistics
################################################################################

# keep only necessary columns
data_sel <- all_data10 %>%
  dplyr::select(
    # key data
    UID,inst_short,inst_short2,taxon_name_acc,species_name_acc,rl_category,rl_year,
    prov_type,num_indiv,acc_num,lat_dd,long_dd,flag,gps_det,
    # locality
    all_locality,locality,municipality,county,state,country,assoc_sp,
    latlong_country,orig_lat,orig_long,uncertainty,orig_prov_type,
    # source
    orig_source,lin_num,coll_num,coll_name,coll_year,coll_date,
    # material info
    germ_type,garden_loc,rec_as,rec_year,rec_date,
    # other metadata
    notes,filename,submission_year,data_source,taxon_name_status,orig_acc_num,fruit_nut,
    # taxon name details
    taxon_name,genus,species,infra_rank,infra_name,cultivar,hybrid,
    taxon_full_name_orig,taxon_full_name_concat,taxon_verif,
    # institution metadata
    inst_id,inst_type,inst_lat,inst_long,inst_country)

# write file
write.csv(data_sel, file.path(main_dir,"outputs",
  paste0("ExSitu_Compiled_", Sys.Date(), ".csv")),row.names = F)


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
    NotH_NoCoords_YesLocal = sum(is.na(lat_dd) & !is.na(all_locality) & prov_type != "H")
  )
head(geo_needs)
# write file
write.csv(geo_needs, file.path(main_dir,"outputs",
  paste0("ExSitu_Geolocation_Needs_Summary_", Sys.Date(), ".csv")),row.names = F)






















##
##### just fruit accessions
##

fruit <- data_sel %>% filter(fruit_nut == "Fruit"); nrow(fruit) #5424
sort(unique(fruit$taxon_name_acc))

# number of institutions who responded to survey
length(sort(unique(fruit$inst_short[which(!grepl("[0-9]",fruit$inst_short))]))) #175
  # number in the US/Canada
NA_inst <- fruit %>% distinct(inst_short,inst_country) %>%
  filter(!grepl("[0-9]",inst_short)) %>% filter(grepl("US|CA",inst_country))
  nrow(NA_inst) #96
# total number (including genebanks)
length(sort(unique(fruit$inst_short))) #205

# map institutions
inst_fruit <- fruit %>%
  select(inst_short,inst_lat,inst_long) %>% distinct()

  map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = inst_fruit,
      lng = ~inst_long, lat = ~inst_lat,
      radius = 9,color = "purple",fillOpacity = 0.6,stroke = F) %>%
    addLegend(
      labels = "Institutions that provided fruit data",
      colors = "purple",position = "bottomright",opacity = 1)
  map

## accession-level stats
  # by taxon and provenance type
fruit_acc <- fruit %>%
  mutate(prov_type = recode(prov_type,
    "N" = "Wild",
    "H?" = "Hort./Unknown",
    "W" = "Wild",
    "Z" = "Wild",
    "NG" = "Hort./Unknown",
    "H" = "Hort./Unknown",
    "U" = "Hort./Unknown")) %>%
  group_by(taxon_name_acc,prov_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
  # add taxon list, to see which have no accessions
fruit_names <- data.frame(taxon_name_acc = sort(unique(taxon_list[taxon_list$fruit_nut=="Fruit",]$taxon_name_acc)))
fruit_names
fruit_acc <- Reduce(full_join,list(fruit_acc,fruit_names))
fruit_acc$num_records[is.na(fruit_acc$num_records)] <- 0
fruit_acc$prov_type[is.na(fruit_acc$prov_type)] <- "Wild"
fruit_acc$taxon_name_acc <- gsub("\\+","x ",fruit_acc$taxon_name_acc)
  # plot
fruit_acc %>%
  mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=prov_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR fruit taxa",
         fill = "Provenance Type")
  #ggsave(file.path(main_dir,"outputs","fruit_prov_barchart.png")#,
  #width = NA, height = NA)
  # look at top species
top_sp <- fruit_acc %>%
  group_by(taxon_name_acc) %>%
  summarize(num_records = sum(num_records)) %>%
  arrange(desc(num_records))
top_sp <- as.data.frame(top_sp)[1:5,1]
for(i in 1:length(top_sp)){
  print(top_sp[i])
  top <- fruit %>%
    filter(taxon_name_acc == top_sp[i])
  print(nrow(top))
  top <- top %>%
    group_by(inst_short) %>% count() %>% ungroup() %>%
    rename(num_records = n) %>% arrange(desc(num_records))
  print(head(as.data.frame(top)))
}
## create table to compare to PNAS results
pnas <- fruit %>%
  group_by(taxon_name_acc) %>%
  count(coord_na = is.na(lat_dd)) %>%
  ungroup()
pnas_wide <- pnas %>%
  spread("coord_na","n")
as.data.frame(pnas_wide)
## by taxon and institution type
fruit_inst <- fruit %>%
  mutate(inst_type = recode(inst_type,
    "Gene/Seed Bank - Genesys/WEIWS" = "Gene/Seed Bank",
    "Botanic Garden; Gene/Seed Bank" = "Unknown")) %>%
  group_by(taxon_name_acc,inst_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
fruit_names
fruit_inst <- Reduce(full_join,list(fruit_inst,fruit_names))
fruit_inst$num_records[is.na(fruit_inst$num_records)] <- 0
fruit_inst$inst_type[is.na(fruit_inst$inst_type)] <- "Unknown"
fruit_inst$taxon_name_acc <- gsub("\\+","x ",fruit_inst$taxon_name_acc)
fruit_inst <- fruit_inst %>% separate("taxon_name_acc","genus",sep=" ",remove=F)
  # plot
fruit_inst %>%
  #mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=inst_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    #facet_wrap(vars(genus), scales="free_y") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR fruit taxa",
         fill = "Institution Type")
## map all points
fruit_pts <- fruit %>%
  distinct(taxon_name_acc,lat_dd,long_dd,.keep_all=T) %>%
  filter(UID != "MEX208~1031700021~NG~Prunus serotina")
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = fruit_pts,
      lng = ~long_dd, lat = ~lat_dd,
      radius = 6,color = "green",fillOpacity = 0.6,stroke = F,
      popup = fruit_pts$UID) %>%
    addControl("North American fruit tree CWR: All wild collection locations",
      position = "topright")
  map
## map specific species as examples
  # get state boundary data
state_bound <- ne_states(country="united states of america")
### Asimina triloba
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Asimina triloba" & !is.na(lat_dd))
  nrow(sp_pts) #29
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Asimina triloba__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Diospyros texana
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Diospyros texana" & !is.na(lat_dd))
  nrow(sp_pts) #12
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Diospyros texana__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Persea borbonia
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Persea borbonia" & !is.na(lat_dd))
  nrow(sp_pts) #4
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Persea borbonia__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Prunus virginiana
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Prunus virginiana" & !is.na(lat_dd))
  nrow(sp_pts) #34
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Prunus virginiana__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)


##
##### just nut accessions
##

nut <- data_sel %>% filter(fruit_nut == "Nut"); nrow(nut) #3829

# number of institutions who responded to survey
length(sort(unique(nut$inst_short[which(!grepl("[0-9]",nut$inst_short))]))) #171
  # number in the US/Canada
NA_inst <- nut %>% distinct(inst_short,inst_country) %>%
  filter(!grepl("[0-9]",inst_short)) %>% filter(grepl("US|CA",inst_country))
  nrow(NA_inst) #93
# total number (including genebanks)
length(sort(unique(nut$inst_short))) #187

# map institutions
inst_nut <- nut %>%
  select(inst_short,inst_lat,inst_long) %>% distinct()

  map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = inst_nut,
      lng = ~inst_long, lat = ~inst_lat,
      radius = 9,color = "purple",fillOpacity = 0.6,stroke = F) %>%
    addLegend(
      labels = "Institutions that provided nut data",
      colors = "purple",position = "bottomright",opacity = 1)
map

  # by taxon and provenance type
nut_acc <- nut %>%
  mutate(prov_type = recode(prov_type,
    "N" = "Wild",
    "H?" = "Hort./Unknown",
    "W" = "Wild",
    "Z" = "Wild",
    "NG" = "Hort./Unknown",
    "H" = "Hort./Unknown",
    "U" = "Hort./Unknown")) %>%
  group_by(taxon_name_acc,prov_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
  # add taxon list, to see which have no accessions
nut_names <- data.frame(taxon_name_acc = sort(unique(taxon_list[taxon_list$fruit_nut=="Nut",]$taxon_name_acc)))
nut_names
nut_acc <- Reduce(full_join,list(nut_acc,nut_names))
nut_acc$num_records[is.na(nut_acc$num_records)] <- 0
nut_acc$prov_type[is.na(nut_acc$prov_type)] <- "Wild"
nut_acc$taxon_name_acc <- gsub("\\+","x ",nut_acc$taxon_name_acc)
nut_inst <- nut_inst %>% separate("taxon_name_acc","genus",sep=" ",remove=F)
  # plot
nut_acc %>%
  mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=prov_type, y=num_records, x=taxon_name_acc)) +
    #geom_bar(position="stack", stat="identity") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR nut taxa",
         fill = "Provenance Type")
  # look at top species
top_sp <- nut_acc %>%
  group_by(taxon_name_acc) %>%
  summarize(num_records = sum(num_records)) %>%
  arrange(desc(num_records))
top_sp <- as.data.frame(top_sp)[1:5,1]
for(i in 1:length(top_sp)){
  print(top_sp[i])
  top <- nut %>%
    filter(taxon_name_acc == top_sp[i])
  print(nrow(top))
  top <- top %>%
    group_by(inst_short) %>% count() %>% ungroup() %>%
    rename(num_records = n) %>% arrange(desc(num_records))
  print(head(as.data.frame(top)))
}
  # create table to compare to PNAS results
pnas <- nut %>%
  group_by(taxon_name_acc) %>%
  count(coord_na = is.na(lat_dd)) %>%
  ungroup()
pnas_wide <- pnas %>%
  spread("coord_na","n")
as.data.frame(pnas_wide)
  # by taxon and instition type
nut_inst <- nut %>%
  mutate(inst_type = recode(inst_type,
    "Gene/Seed Bank - Genesys/WEIWS" = "Gene/Seed Bank",
    "Botanic Garden; Gene/Seed Bank" = "Unknown")) %>%
  group_by(taxon_name_acc,inst_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
nut_names
nut_inst <- Reduce(full_join,list(nut_inst,nut_names))
nut_inst$num_records[is.na(nut_inst$num_records)] <- 0
nut_inst$inst_type[is.na(nut_inst$inst_type)] <- "Unknown"
nut_inst$taxon_name_acc <- gsub("\\+","x ",nut_inst$taxon_name_acc)
nut_inst <- nut_inst %>% separate()
  # plot
nut_inst %>%
  #mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=inst_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    #facet_grid(vars(genus), scales="free_y") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR nut taxa",
         fill = "Institution Type")
  # map all points
nut_pts <- nut %>%
  distinct(taxon_name_acc,lat_dd,long_dd,.keep_all=T) #%>%
  #filter(UID != "MEX208~1031700021~NG~Prunus serotina")
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = nut_pts,
      lng = ~long_dd, lat = ~lat_dd,
      radius = 6,color = "green",fillOpacity = 0.6,stroke = F,
      popup = nut_pts$UID)
  map
## map specific species as examples
  # get state boundary data
#state_bound <- ne_states(country="united states of america")
### Juglans microcarpa
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Juglans microcarpa" & !is.na(lat_dd))
  nrow(sp_pts) #10
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Juglans microcarpa__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Corylus cornuta
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Corylus cornuta" & !is.na(lat_dd))
  nrow(sp_pts) #18
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Corylus cornuta__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Castanea pumila
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Castanea pumila" & !is.na(lat_dd))
  nrow(sp_pts) #9
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Castanea pumila__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Carya glabra
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Carya glabra" & !is.na(lat_dd))
  nrow(sp_pts) #26
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Carya glabra__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)





















### !!!! ENDED HERE FOR NOW !!!!



##
## SELECT AND ORDER FINAL COLUMNS
##

all_data10 <- as.data.frame(lapply(all_data10, function(x) str_squish(x)),
  stringsAsFactors=F)
all_data10 <- as.data.frame(lapply(all_data10, function(x) gsub(",",";",x)),
  stringsAsFactors=F)

all_data13 <- all_data10 %>%
  ### combine duplicates at all_locality level
  group_by(inst_short,species_name_acc,prov_type,all_locality) %>%
  mutate(
    UID = paste(UID, collapse="|"),
    notes = paste(unique(notes), collapse="; "),
    assoc_sp = paste(unique(assoc_sp), collapse="; "),
    acc_num = paste(acc_num, collapse="|"),
    sum_num_indiv = sum(as.numeric(num_indiv)),
    germ_type = paste(unique(germ_type), collapse="; "),
    garden_loc = paste(unique(garden_loc), collapse="; "),
    rec_as = paste(unique(rec_as), collapse="; "),
    taxon_det = paste(unique(taxon_det), collapse="; "),
    taxon_name_acc = paste(unique(taxon_name_acc), collapse="; "),
    taxon_full_name = paste(unique(taxon_full_name), collapse="; "),
    taxon_full_name_orig = paste(unique(taxon_full_name_orig), collapse="; "),
    taxon_full_name_concat = paste(unique(taxon_full_name_concat), collapse="; "),
    cultivar = paste(unique(cultivar), collapse="; "),
    sum_num_acc = n()) %>%
  ungroup() %>%
  distinct(inst_short,species_name_acc,prov_type,all_locality,.keep_all=T) %>%
  dplyr::select(
    # grouping data
    inst_short,species_name_acc,prov_type,all_locality,
    # key data
    UID,gps_det,flag,lat_dd,long_dd,
    # locality
    locality,municipality,county,state,country,latlong_country,
    orig_source,notes,orig_lat,orig_long,assoc_sp,
    # source
    acc_num,lin_num,coll_num,coll_name,coll_year,
    # material info
    sum_num_indiv,sum_num_acc,germ_type,garden_loc,rec_as,taxon_det,
    # taxon name
    list,taxon_name_acc,taxon_full_name,genus,
    taxon_full_name_orig,taxon_full_name_concat,cultivar,
    # species metadata
    rl_year,rl_category,
    # institution metadata
    inst_name,inst_country,inst_lat,inst_long,filename)
nrow(all_data13) #13867
head(as.data.frame(all_data13))

# write file
write.csv(all_data13, file.path(main_dir,"outputs",
  paste0("ExSitu_Compiled_Standardized_", Sys.Date(), ".csv")),row.names = F)










### NOT USING YET/CURRENTLY ###


##
## RENAME FOR GEOLOCATE AND (optionally.. SPLIT BY SPECIES)
##

# FIRST CHECK TO BE SURE THIS IS ZERO !!
all_data10[which(grepl("\\|",all_data10$acc_num)),]

  # add GEOLocate standard columns
all_data10$correction.status <- NA
all_data10$precision <- NA
all_data10$error.polygon <- NA
all_data10$multiple.results <- NA
all_data10$uncertainty <- NA

all_data14 <- all_data10 %>%
  # filter to remove cultivated records and those without locality info
  #filter(rl_category == "CR" | rl_category == "EN" |
  #       rl_category == "VU" | rl_category == "NT") %>%
  #filter(prov_type != "H") %>%
  #filter(!is.na(all_locality)) %>%
  # rename to GEOLocate standard columns
  rename(locality.string = all_locality) %>%
  # order with NA lat-long records on top
  arrange(locality.string) %>%
  arrange(!is.na(latitude),latitude) %>%
  # replace NA with "" to make simpler to view in GEOLocate
  replace(., is.na(.), "") %>%
  # group by all non-ID fields
  group_by(
    locality.string,country,state,county,latitude,longitude,
    flag,gps_det,prov_type,lin_num,coll_num,coll_name,coll_year,
    inst_short,filename,inst_lat,inst_long,
    list,species_name_acc,taxon_full_name) %>%
  # concatenate values in ID fields
  mutate(
    UID = paste(UID, collapse="|"),
    acc_num = paste(acc_num, collapse="|"),
    sum_num_indiv = sum(as.numeric(num_indiv))) %>%
  ungroup() %>%
  # remove duplicates
  distinct(
    locality.string,country,state,county,latitude,longitude,
    flag,gps_det,prov_type,lin_num,coll_num,coll_name,coll_year,
    inst_short,filename,inst_lat,inst_long,
    list,species_name_acc,taxon_full_name,
    .keep_all=T) %>%
  # reorder columns
  dplyr::select(
    ## GeoLocate
    locality.string,country,state,county,latitude,longitude,
    correction.status,precision,error.polygon,multiple.results,uncertainty,
    ## record metadata
    flag,gps_det,prov_type,acc_num,lin_num,coll_num,coll_name,coll_year,sum_num_indiv,
    ## institituion metadata
    inst_short,filename,inst_lat,inst_long,
    ## taxon name & record ID
    list,species_name_acc,taxon_full_name,UID) %>%
  # rename concatenated fields to make that clear
  rename(
    acc_num_CAT = acc_num,
    UID_CAT = UID)

# remove dot in column names (replace with space) for GEOLocate
names(all_data14) <- gsub(x = names(all_data14),pattern = "\\.",
  replacement = " ")
str(all_data14)
head(as.data.frame(all_data14),n=30)

# write file
write.csv(all_data14, file.path(main_dir,"outputs","to_geolocate",
  paste0("To_Geolocate_CR-EN-VU-NT_", Sys.Date(), ".csv")),row.names = F)


# create one CSV for each target species
sp_split <- split(all_data14, as.factor(all_data14$species_name_acc))
names(sp_split) <- gsub(" ","_",names(sp_split))

# write files
if(!dir.exists(file.path(main_dir,"outputs","to_geolocate")))
  dir.create(file.path(main_dir,"outputs","to_geolocate"),
  recursive=T)
lapply(seq_along(sp_split), function(i) write.csv(sp_split[[i]],
  file.path(main_dir,"outputs","to_geolocate",
  paste0(names(sp_split)[[i]], ".csv")),row.names = F))




### Read geolocated CSVs back in and add geolocate info to rest of data

# read in geolocated CSVs
file_list <- list.files(
  path=file.path(main_dir,"Compiled ex situ data","Geolocated_CSV_by_target_species",target_genus),
  pattern=".csv",full.names=TRUE)
file_dfs <- lapply(file_list,read.csv,header=TRUE,fileEncoding="LATIN1",
  strip.white=TRUE,stringsAsFactors=F,na.strings=c("","NA"))
# check that geolocated CSVs each have same number of rows as inital export;
# if they don't, you'll need to manually see where the mistake is
for(i in 1:length(sp_split)){
  for(j in 1:nrow(sp_split[[i]])){
    if(sp_split[[i]]$inst_short[j] !=file_dfs[[i]]$inst_short[j]){
      print(sp_split[[i]]$species_name_acc[1])
    } else {
      print(i)
    }
  }
}

# bind geolocated CSVs together
post_geo <- Reduce(rbind.fill, file_dfs)
# fix a few inconsistencies
  ## provenance type column
unique(post_geo$prov_type)
    # check "H?" rows to see if should be "W" and
    # if all are "X" gps_det, change prov_type to "H"
post_geo[which(post_geo$prov_type == "H?"),]
post_geo[which(post_geo$prov_type == "H?"),]$prov_type <- "H"
    # check prov_type for rows with coordinates
unique(post_geo[which(!is.na(post_geo$latitude)),]$prov_type)
post_geo[which(!is.na(post_geo$latitude) & post_geo$prov_type == "U"),]
  ## gps determination column
unique(post_geo$gps_det)
post_geo[which(is.na(post_geo$gps_det)),]$latitude
post_geo[which(is.na(post_geo$gps_det)),]$gps_det <- "X"
  ## uncertainty column
post_geo$uncertainty <- gsub(" m","",post_geo$uncertainty)
post_geo[which(post_geo$uncertainty == "0"),]$uncertainty <- NA
post_geo$uncertainty <- as.numeric(post_geo$uncertainty)
  ## lat and long
sort(unique(post_geo$latitude))
sort(unique(post_geo$longitude))

# keep only edited columns
geolocated <- post_geo %>%
  rename(UID = UID_CAT) %>%
  dplyr::select(country,state,county,latitude,longitude,precision,
    uncertainty,gps_det,prov_type,UID) #multiple.results
head(as.data.frame(geolocated),n = 40)

## for Quercus, because geolocated before changed UID system
# bind itial export together, to join UID to geolocated data
#pre_geo <- Reduce(rbind.fill, sp_split)
#pre_geo <- pre_geo %>% dplyr::select(UID_CAT)
# add UID to geolocated rows and separate to accession level again
#geolocated <- cbind(post_geo,pre_geo)

# separate to accession level again
geolocated2 <- separate_rows(geolocated, UID, sep="\\|")
geolocated2 <- separate_rows(geolocated2, UID, sep="; ")
geolocated2 <- geolocated2 %>%
  rename(lat_dd = latitude, long_dd = longitude,
    coord_precision = uncertainty)

# see if all UIDs are matching up
yes_geo <- all_data13 %>% filter(UID %in% geolocated2$UID)
  # should be character(0)
setdiff(geolocated2$UID,yes_geo$UID)
# add gelocated rows to data
yes_geo <- yes_geo %>% dplyr::select(-prov_type,-gps_det,-lat_dd,
  -long_dd,-county,-state,-country,-coord_precision)
yes_geo <- full_join(yes_geo,geolocated2)
# bind all data together (geolocated and garden origin)
no_geo <- all_data13 %>% filter(!(UID %in% geolocated2$UID))
no_geo$precision <- NA
all_data15 <- rbind(yes_geo,no_geo)

# make gps_det "X" if NA -- not doing this now to distinquish records
#   that have been checked to see if can geolocate versus those that haven't
#all_data15[which(is.na(all_data15$gps_det)),]$gps_det <- "X"
unique(all_data15$gps_det)

# arrange columns
all_data15 <- all_data15 %>%
  arrange(species_name_acc,UID) %>%
  rename(gps_notes = flag) %>%
  dplyr::select(
    # key data
    UID,inst_short,submission_year,species_name_acc,target_species,
    prov_type,gps_det,lat_dd,long_dd,coord_precision,precision,gps_notes,
    all_locality,
    # locality
    locality,municipality,county,state,country,latlong_country,
    orig_lat,orig_long,
    orig_source,notes,assoc_sp,habitat,num_indiv,acc_num,
    # source
    lin_num,coll_num,coll_name,coll_year,
    # material info
    germ_type,garden_loc,rec_as,condition,name_determ,
    # other metadata
    dataset_year,private,filename,list,
    # taxon name
    taxon_name_acc,taxon_full_name,genus,species,infra_rank,infra_name,
    taxon_full_name_orig,taxon_full_name_concat,cultivar,trade_name,
    # institution metadata
    inst_name,inst_country,inst_lat,inst_long)
all_data15[which(all_data15$gps_notes == ""),]$gps_notes <- NA
str(all_data15)

# write file
write.csv(all_data15, file.path(main_dir,"Compiled ex situ data",
  "ALL DATA - POST GEOLOCATE",
  paste0(target_genus,"_POSTGEO_exsitu_compiled_standardized.csv")),
  row.names = F)

# create one CSV for each target species
all_data_target <- all_data15 %>% filter(target_species == "Y")
sp_split2 <- split(all_data_target,
  as.factor(all_data_target$species_name_acc))
names(sp_split2) <- gsub(" ","_",names(sp_split2))

# write files
if(!dir.exists(file.path(main_dir,"Compiled ex situ data",
  "ALL DATA - POST GEOLOCATE","target_species")))
  dir.create(file.path(main_dir,"Compiled ex situ data",
  "ALL DATA - POST GEOLOCATE","target_species"),recursive=T)
lapply(seq_along(sp_split2), function(i) write.csv(sp_split2[[i]],
  file.path(main_dir,"Compiled ex situ data","ALL DATA - POST GEOLOCATE",
  "target_species",
  paste0(names(sp_split2)[[i]], "_ALL_POSTGEO.csv")),row.names = F))












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
