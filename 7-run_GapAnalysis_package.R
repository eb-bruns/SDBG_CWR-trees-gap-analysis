################################################################################

## 7-run_GapAnalysis_package.R

### Author: Emily Beckman Bruns
### Funding: Cooperative agreement between the United States Botanic Garden and 
# San Diego Botanic Garden (subcontracted to The Morton Arboretum), with 
# support from Botanic Gardens Conservation International U.S.

### Last updated: 29 January 2023
### R version 4.2.2

### DESCRIPTION:
# 

### INPUTS:
# 

### OUTPUTS:
# 

################################################################################
# Load libraries
################################################################################

my.packages <- c(# additional packages for mapping and other visualizations
                  'leaflet','RColorBrewer','ggplot2',
                  'forcats','grid','downloadthis',#'png','gridExtra',
                 # suggested packages for GapAnalysis workflow
                 'GapAnalysis',
                 'dplyr', 'sp', 'tmap', 'data.table', 'sf', 'methods',
                 'geosphere', 'fasterize', 'rmarkdown', 'knitr', 
                 'kableExtra', 'DT',#'rgdal', 'rgeos', 
                 # additional packages for mapping and other visualizations
                  'terra','raster'
                 )
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)

# source my updated versions of GapAnalysis functions
source("/Users/emily/Documents/GitHub/GapAnalysis/R/SRSin.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/GRSin.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/ERSin.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/SRSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/GRSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/ERSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSin.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSc_mean.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/SummaryHTML.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/OccurrenceCounts.R")

source("/Users/emily/Documents/GitHub/GapAnalysis/R/LeafletMapGRSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/SummaryScoresChart.R")

# source function for filtering occurrence points based on specific flags
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/filter_points.R")


################################################################################
# Set working directory
################################################################################

# assign main working directory
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")

# set up additional file paths
path.pts.raw <- file.path(main_dir,occ_dir,"raw_occurrence_data")
path.pts <- file.path(main_dir,occ_dir,"standardized_occurrence_data",
                      "taxon_edited_points")
path.sdm <- file.path(main_dir,"sdm_rasters")
path.rastbuff <- file.path(main_dir,"rasterized_buffers")

################################################################################
# Read in data for all target taxa
################################################################################

## first we will define which columns to keep in raw output occurrence datasets
## (the ones that will be downloadable as a whole for all taxa)
keep_col <- c( #data source and unique ID
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
  #location flags
  "latlong_countryCode",".cen",".urb",".inst",".outl",".nativectry",
  ".yr1950",".yr1980",".yrna")

### Protected areas raster
# this seems to download the 10 arc minutes version instead of 2.5 arc min...
#data(ProtectedAreas)
#str(ProtectedAreas)
#res(ProtectedAreas) # resolution = 0.1666667 0.1666667
# ...so I downloaded from the source to read that in instead...
# read in PA file downloaded from Dataverse (https://dataverse.harvard.edu/dataverse/GapAnalysis)
ProtectedAreas <- raster(file.path(main_dir,gis_dir,"wdpa_reclass.tif"))
res(ProtectedAreas) # resolution = 0.04166667 0.04166667

### Global ecoregions shapefile
# this seems to download a version that is cropped to the western US ?...
#data(ecoregions)
#head(ecoregions)
# ...so I downloaded from the source to read that in instead...
# read in as simple features
ecoregions_sf <- sf::read_sf(file.path(
  main_dir,gis_dir,"terr-ecoregions-TNC/tnc_terr_ecoregions.shp"))
# convert to SpatialPolygonsDataFrame
ecoregions <- as_Spatial(ecoregions_sf)
# crop to target regions, to make a little smaller
unique(ecoregions$WWF_REALM2)
ecoregions <- ecoregions[ecoregions$WWF_REALM2 == "Nearctic" | 
                         ecoregions$WWF_REALM2 == "Neotropic",]
# can check visually if desired: plot(ecoregions)

## Land boundary shapefile for North America 
# (created from ecoregions with Great Lakes clipped out)
# we use this to clip buffers when mapping so they're not in the water
boundary.poly <- vect(file.path(main_dir,gis_dir,"NorthAm_land_boundary",
                                "NorthAm_land_boundary.shp"))

### Occurrence data
# from package example: data(CucurbitaData)
# my occurrence data compiled through full workflow:
## read in raw occurrence data (before keeping only geolocated & removing dups)
#  and edit to match format needed for GapAnalysis
all_occ_raw <-  read.csv(file.path(path.pts.raw,
                                  "all_occurrence_data_raw_2023-01-04.csv"))
  ## save version without ex situ for people to download
  #all_occ_raw_noG <- all_occ_raw %>% dplyr::filter(database!="Ex_situ")
  #rm_col <- c("latlong_countryCode",".cen",".urb",".inst",".outl",".nativectry",
  #             ".yr1950",".yr1980",".yrna","all_source_databases",
  #             "countryCode_standard")  
  #raw_keep <- keep_col[!keep_col %in% rm_col]
  #all_occ_raw_noG <- all_occ_raw_noG[,raw_keep]
  #write.csv(all_occ_raw_noG, file.path(
   # main_dir,occ_dir,"raw_occurrence_data",
  #  "All_RefRecords_Jan2023_NAFruitNutCWR.csv"), row.names = F)
all_occ_raw <- all_occ_raw %>%
  dplyr::mutate(type = recode(database,
                              "Ex_situ" = "G",
                              .default = "H")) %>%
  dplyr::select(-species) %>%
  dplyr::rename(species = taxon_name_accepted,
                latitude = decimalLatitude,
                longitude = decimalLongitude) %>%
  dplyr::select(species,latitude,longitude,type,
                UID,database,
                taxon_name,
                basisOfRecord,year,
                nativeDatabaseID,institutionCode,datasetName,datasetName,
                references,recordedBy,
                establishmentMeans,individualCount,
                coordinateUncertaintyInMeters,localityDescription,
                iucnredlist_category,natureserve_rank)
all_occ_raw$species <- gsub(" ","_",all_occ_raw$species)
## read in flagged occurrence data (ready for spatial analyses)
occ_files <- list.files(path.pts, pattern = ".csv", full.names = T)
occ_dfs <- lapply(occ_files, read.csv)
all_occ <- Reduce(rbind, occ_dfs)
rm(occ_files,occ_dfs)
# read in manual edits to occurrence data (flagging additional bad points, etc)
manual_pt_edits <- read.csv(file.path(main_dir,occ_dir,
                                      "manual_point_edits.csv"),
                            na.strings = c("NA",""),
                            colClasses = "character")
# filter points based on flagging columns and manual edits
all_occ <- filter.points(all_occ,manual_pt_edits)
  ## save version of filtered points for our records
  #write.csv(all_occ, file.path(
  #  main_dir,occ_dir,"standardized_occurrence_data",
  #  paste0("all_occurrence_data_filtered_",Sys.Date(),".csv")),row.names = F)
  ## save version without ex situ for people to download
  #all_occ_noG <- all_occ %>% dplyr::filter(database!="Ex_situ")
  #all_occ_noG <- all_occ_noG[,keep_col]
  #write.csv(all_occ_noG, file.path(
  #  main_dir,occ_dir,"standardized_occurrence_data",
  #  "All_RefCoordsMapped_Jan2023_NAFruitNutCWR.csv"), row.names = F)
# edit occurrence data to match format needed for GapAnalysis
all_occ <- all_occ %>%
  dplyr::mutate(type = recode(database,
                              "Ex_situ" = "G",
                              .default = "H")) %>%
  dplyr::rename(species = taxon_name_accepted,
                latitude = decimalLatitude,
                longitude = decimalLongitude) %>%
  dplyr::select(species,latitude,longitude,type,
                UID,database,all_source_databases,
                taxon_name,
                basisOfRecord,year,
                nativeDatabaseID,institutionCode,datasetName,datasetName,
                references,recordedBy,
                establishmentMeans,individualCount,
                coordinateUncertaintyInMeters,localityDescription,
                countryCode_standard)
all_occ$species <- gsub(" ","_",all_occ$species)
str(all_occ)
per_sp <- all_occ %>% count(species); per_sp

### Species distribution rasters
# from package example: data(CucurbitaRasters); CucurbitaRasters <- raster::unstack(CucurbitaRasters)
# read in own SDM rasters (created from forked CWR-of-the-USA-Gap-Analysis repo)
sdm_files <- list.files(path.sdm, pattern = ".tif", full.names = T)
sdm_list <- lapply(sdm_files, raster)

### Create species list
# from package example: speciesList <- unique(CucurbitaData$species)
# from own data
speciesList_sdm <- gsub("-spdist_thrsld_median.tif","",
                        list.files(path.sdm,pattern=".tif",full.names=F))

# read in own rasterized buffers (created in create_distribution_rasters.R)
#buff20_files <- list.files(file.path(path.rastbuff,"20km"),pattern=".tif",full.names=T)
#buff20_list <- lapply(buff20_files, raster)
buff50_files <- list.files(file.path(path.rastbuff,"50km"),pattern=".tif",full.names=T)
  # select just those that dont have SDM, for now:  
  buff50_files <- buff50_files[!(
    gsub("-rasterized_buffers_50km.tif","",list.files(
      file.path(path.rastbuff,"50km"),pattern=".tif",full.names=F)) 
    %in% speciesList_sdm)]
buff50_list <- lapply(buff50_files, raster)
#buff100_files <- list.files(file.path(path.rastbuff,"100km"),pattern=".tif",full.names=T)
#buff100_list <- lapply(buff100_files, raster)
speciesList_buff <- gsub("-rasterized_buffers_50km.tif","",
                        list.files(file.path(path.rastbuff,"50km"),
                                   pattern=".tif",full.names=F))
speciesList_buff <- speciesList_buff[!(speciesList_buff %in% speciesList_sdm)]

## stack those with SDM and those without (buffers instead)
speciesList <- c(speciesList_sdm,speciesList_buff)
rasterList <- c(sdm_list,buff50_list)

## make a list that corresponds with each raster layer, to tell if its an SDM or buffers
rasterType <- c(rep("SDM", 
                    times = length(speciesList_sdm)),
                rep("buffers", 
                    times = length(speciesList_buff)))

################################################################################
# Run SummaryHTML function
################################################################################

# create html outputs !
SummaryHTML(
  Species_list=speciesList, 
  Occurrence_data=all_occ, 
  Raster_list=rasterList,
  Buffer_distance=50000, 
  Ecoregions_shp=ecoregions,
  Pro_areas=ProtectedAreas, 
  Output_Folder=file.path(main_dir,"html_outputs"),
  writeRasters=T,
  Occurrence_data_raw=all_occ_raw,
  boundary=boundary.poly,
  Raster_type=rasterType
)

### USED FOR TESTING WITHIN Rmd AND INDIVIDUAL FUNCTIONS:
#Species_list <- speciesList
#Occurrence_data <- all_occ
#Occurrence_data_raw <- all_occ_raw
#Raster_list <- rasterList
#Buffer_distance <- 50000
#Output_Folder <- file.path(main_dir,"html_outputs")
#writeRasters <- T
#Raster_type <- rasterType

#Ecoregions_shp <- ecoregions
#Pro_areas <- ProtectedAreas
#pt.proj <<- "EPSG:4326"
#boundary <- boundary.poly
#Gap_Map <- T
#Buffer_distance <- 50000
#i <- 5
#Sl <- speciesList[i]
#Od <- all_occ[all_occ$species == speciesList[i], ]
#OdR <- all_occ_raw[all_occ_raw$species == speciesList[i], ]
#Rl <- rasterList[[i]]
#Rt <- rasterType[[i]]

################################################################################
# Get analysis results and visualize
################################################################################

## run using SDMs or buffers, depending on which there is for the taxon
  ## Ex situ analyses
FCSex_all <- FCSex(Species_list=speciesList,
                   Occurrence_data=all_occ,
                   Raster_list=rasterList,
                   Buffer_distance=50000,
                   Ecoregions_shp=ecoregions,
                   Gap_Map=FALSE,
                   Occurrence_data_raw=all_occ_raw)
#  FCSex_sdm$SpeciesDist <- "SDM"
#  FCSex_sdm$ExsituBuffer <- "50km"
  ## In situ analyses
FCSin_all <- FCSin(Species_list=speciesList,
                   Occurrence_data=all_occ,
                   Raster_list=rasterList,
                   Ecoregions_shp=ecoregions,
                   Pro_areas=ProtectedAreas,
                   Gap_Map=FALSE)
#  FCSin_sdm$SpeciesDist <- "SDM"
#  FCSin_sdm$ExsituBuffer <- "N/A"
  ## Calculate means then join everything together for summary figure
FCSc_mean <- FCSc_mean(FCSex_all, FCSin_all)  
FCSc_all <- Reduce(full_join,list(FCSc_mean_sdm,FCSex_sdm[,1:4],FCSin_sdm[,1:4]))
  # add star to those where buffer was used instead of SDM
  FCSc_all[which(FCSc_all$Taxon %in% speciesList_buff),"Taxon"] <- 
    paste0(FCSc_all[which(FCSc_all$Taxon %in% speciesList_buff),"Taxon"]," ★")
## Create summary figure
final_chart <- SummaryScoresChart(FCSc_all)
ggsave(file.path(main_dir,"All-GapAnalysis-Scores-Chart.png"),
       plot = final_chart, device = "png", units = "in", width = 10, height = 12, dpi = 320)


################################################################################
# Compare results using different 
################################################################################






  
  
  
  




# function to run gap analysis calculations for multiple buffer sizes and 
#   multiple species distribution rasters
rep_gap_analysis <- function(Species_list, Occurrence_data, Raster_list,
                             Ecoregions_shp, Pro_areas,
                             Occurrence_data_raw,
                             reps){
  all_results <- data.frame()
  for(i in 1:nrow(reps)){
    print(paste0("----- NOW STARTING ANALYSES USING: ",reps$SpeciesDist[i]," -----"))
    ## run all ex situ analyses
    ex_results <- data.frame()
    print("--- STARTING EX SITU ANALYSES ---")
    print(paste0("-- BUFFER SIZE FOR G POINTS: ",reps$GBufferSize[i]/1000,"KM --"))
    FCSex_df <- FCSex(Species_list=Species_list,
                      Occurrence_data=Occurrence_data,
                      Raster_list=reps$Raster_list[[i]],
                      Buffer_distance=reps$GBufferSize[i],
                      Ecoregions_shp=Ecoregions_shp,
                      Occurrence_data_raw=Occurrence_data_raw)
    FCSex_df$CalcType <- "Ex situ"
    FCSex_df$SpeciesDist <- reps$SpeciesDist[i] 
    FCSex_df$GBufferSizeKm <- (reps$GBufferSize[i]/1000)
    ## rename columns to remove ex (we have column for that instead)
    colnames(FCSex_df)[2:6] <- c("SRS","GRS","ERS","FCS","FCS_class")
    # rename columns to include buffer size used (this was for wide df)
    #colnames(FCSex_df)[2:6] <- c(
    #  paste0("SRSex",spdist_type[i]),
    #  paste0("GRSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
    #  paste0("ERSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
    #  paste0("FCSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
    #  paste0("FCSex_class",buffer_sizes[j]/1000,"_",spdist_type[i]))
    ## add to all results dataframe
    ex_results <- rbind(ex_results,FCSex_df)
    ## run all in situ analyses
    print("--- STARTING IN SITU ANALYSES ---")
    in_results <- FCSin(Species_list=Species_list,
                        Occurrence_data=Occurrence_data,
                        Raster_list=reps$Raster_list[[i]],
                        Ecoregions_shp=Ecoregions_shp,
                        Pro_areas=Pro_areas)
    in_results$CalcType <- "In situ"
    in_results$SpeciesDist <- reps$SpeciesDist[i] 
    in_results$GBufferSizeKm <- "N/A"
    ## rename columns to remove ex (we have column for that instead)
    colnames(in_results)[2:6] <- c("SRS","GRS","ERS","FCS","FCS_class")
    # rename columns to include buffer size used (this was for wide df)
    #colnames(in_results)[2:6] <- c(
    #  paste0("SRSin_",spdist_type[i]),
    #  paste0("GRSin_",spdist_type[i]),
    #  paste0("ERSin_",spdist_type[i]),
    #  paste0("FCSin_",spdist_type[i]),
    #  paste0("FCSin_class_",spdist_type)[i])
    # combine ex situ and in situ analyses in one df
    all_results_now <- rbind(ex_results,in_results)
    all_results <- rbind(all_results,all_results_now)
  }
  return(all_results)
}

## read in other buffer sizes (we just used 50km above)
  # 20km  
  buff20_files <- list.files(file.path(path.rastbuff,"20km"),pattern=".tif",full.names=T)
  buff20_files <- buff20_files[
    gsub("-rasterized_buffers_20km.tif","",list.files(
      file.path(path.rastbuff,"20km"),pattern=".tif",full.names=F)) 
    %in% speciesList_sdm]
  buff20_list <- lapply(buff20_files, raster)
  # 100km  
  buff100_files <- list.files(file.path(path.rastbuff,"100km"),pattern=".tif",full.names=T)
  buff100_files <- buff100_files[
    gsub("-rasterized_buffers_100km.tif","",list.files(
      file.path(path.rastbuff,"100km"),pattern=".tif",full.names=F)) 
    %in% speciesList_sdm]
  buff100_list <- lapply(buff100_files, raster)
  
## now run just for species with SDM
reps_matrix_yesSDM <- data.frame(
  SpeciesDist = c("20km buffers",
                  "50km buffers",
                  "100km buffers",
                  "SDM"),
  GBufferSize = c(20000,
                  50000,
                  100000,
                  50000))
  RasterList <- list(buff20_list,
                     buff50_list,
                     buff100_list,
                     sdm_list)

all_results_yesSDM_now <- rep_gap_analysis(
  Species_list = speciesList_sdm, 
  Occurrence_data = all_occ, 
  Rasterized_buffers_list = buff100_list,
  Ecoregions_shp = ecoregions, 
  Pro_areas = ProtectedAreas, 
  Occurrence_data_raw = all_occ_raw, 
  reps = reps_matrix_yesSDM)

### charts to visualize results

# all SRS results
ggplot(all_results_yesSDM, 
       aes(x = SRS, y = species)) + 
  geom_point(color="black") +
  xlim(0,100) +
  facet_grid(cols = vars(CalcType)) +
  labs(title = "Sampling Representativeness Score (SRS)")

# all FCS results
all_results_yesSDM %>% 
  mutate(SpeciesDist = fct_relevel(SpeciesDist, 
                                   "SDM","20km buffers",
                                   "50km buffers","100km buffers")) %>%
  ggplot(aes(x = FCS, y = species, col = SpeciesDist)) + 
    geom_point() +
    xlim(0,100) +
    facet_grid(cols = vars(CalcType)) +
    labs(title = "Final Conservation Score (FCS; average of SRS, GRS, and ERS)")
#ggplot(all_results_yesSDM, 
#       aes(x = FCS, y = SpeciesDist, col = GBufferSizeKm)) + 
#  geom_point() +
#  xlim(0,100) +
#  facet_grid(rows = vars(species), cols = vars(CalcType)) +
#  labs(title = "Final Conservation Score (FCS; average of SRS, GRS, and ERS)")

# all GRS results
all_results_yesSDM %>% 
  mutate(SpeciesDist = fct_relevel(SpeciesDist, 
                                   "SDM","20km buffers",
                                   "50km buffers","100km buffers")) %>%
  ggplot(aes(x = GRS, y = species, col = SpeciesDist)) + 
    geom_point() +
    xlim(0,100) +
    facet_grid(cols = vars(CalcType)) +
    labs(title = "Geographical Representativeness Score (GRS)")

# all ERS results
all_results_yesSDM %>% 
  mutate(SpeciesDist = fct_relevel(SpeciesDist, 
                                   "SDM","20km buffers",
                                   "50km buffers","100km buffers")) %>%
  ggplot(aes(x = ERS, y = species, col = SpeciesDist)) + 
    geom_point() +
    xlim(0,100) +
    facet_grid(cols = vars(CalcType)) +
    labs(title = "Ecological Representativeness Score (ERS)")

#all_results_ex <- all_results_yesSDM %>% filter(CalcType=="Ex situ")
#all_results_in <- all_results_yesSDM %>% filter(CalcType=="In situ")









################################################################################
# Functions
################################################################################

# !!! TESTING !!!
#  Species_list <- speciesList_yesSDM
#  Occurrence_data <- allOcc_yesSDM
#  Raster_list <- RasterList_yesSDM[[1]]
#  Rasterized_buffers_list <- buff50_list_yesSDM
#  Ecoregions_shp <- ecoregions
#  Pro_areas <- ProtectedAreas
#  Occurrence_data_raw <- allOccRaw_yesSDM
#  Select_database <- "GBIF"
#  reps <- reps_matrix_yesSDM
#  Gap_Map <- F
#  Buffer_distance <- 50000





















# for species without SDM, we will just the three rasterized buffer layers 
# created above:
#spdistType_noSDM <- c("A) 20km buffers","B) 50km buffers","C) 100km buffers")
#RasterList_noSDM <- list(buff20_list,buff50_list,buff100_list)






### TESTING WITH SDM

##Running all three ex situ gap analysis steps using FCSex function
## repeat with multiple buffer sizes, as desired
buffer_sizes <- c(20000,50000,100000)
ex_results_sdm <- data.frame(species=speciesList)
for(i in 1:length(buffer_sizes)){
  # run analyses
  FCSex_df <- FCSex(Species_list=speciesList,
                    Occurrence_data=JData,
                    Raster_list=JRasters,
                    Buffer_distance=buffer_sizes[i],
                    Ecoregions_shp=ecoregions,
                    Occurrence_data_raw=JDataRaw,
                    Select_database="IUCN_RedList")
  # rename columns to include buffer size used
  # we skip SRS because it doesn't use the buffer size (so same for all)
  colnames(FCSex_df)[3:6] <- c(
    paste0("GRSex_",buffer_sizes[i]/1000),
    paste0("ERSex_",buffer_sizes[i]/1000),
    paste0("FCSex_",buffer_sizes[i]/1000),
    paste0("FCSex_class_",buffer_sizes[i]/1000)
  )
  # add to all results dataframe
  ex_results_sdm <- full_join(ex_results_sdm,FCSex_df)
}

##Running all three in situ gap analysis steps using FCSin function
in_results_sdm <- FCSin(Species_list=speciesList,
                        Occurrence_data=JData,
                        Raster_list=JRasters,
                        Ecoregions_shp=ecoregions,
                        Pro_areas=ProtectedAreas,
                        Occurrence_data_raw=JDataRaw,
                        Select_database="IUCN_RedList")


### TESTING WITH RASTERIZED BUFFERS INSTEAD OF SDM

##Running all three ex situ gap analysis steps using FCSex function
## repeat with multiple buffer sizes, as desired
buffer_sizes <- c(20000,50000,100000)
#rasterized_buffers <- c(JDataBuff20,JDataBuff50,JDataBuff100)
ex_results_buff <- data.frame(species=speciesList)
for(i in 1:length(buffer_sizes)){
  # run analyses
  FCSex_df <- FCSex(Species_list=speciesList,
                    Occurrence_data=JData,
                    Raster_list=JDataBuff20, # !! CAN CHANGE THIS TOO
                    Buffer_distance=buffer_sizes[i],
                    Ecoregions_shp=ecoregions,
                    Occurrence_data_raw=JDataRaw,
                    Select_database="IUCN_RedList")
  # rename columns to include buffer size used
  # we skip SRS because it doesn't use the buffer size (so same for all)
  colnames(FCSex_df)[3:6] <- c(
    paste0("GRSex_",buffer_sizes[i]/1000),
    paste0("ERSex_",buffer_sizes[i]/1000),
    paste0("FCSex_",buffer_sizes[i]/1000),
    paste0("FCSex_class_",buffer_sizes[i]/1000)
  )
  # add to all results dataframe
  ex_results_buff <- full_join(ex_results_buff,FCSex_df)
}

##Running all three in situ gap analysis steps using FCSin function
in_results_buff <- FCSin(Species_list=speciesList,
                         Occurrence_data=JData,
                         Raster_list=JDataBuff50,
                         Ecoregions_shp=ecoregions,
                         Pro_areas=ProtectedAreas,
                         Occurrence_data_raw=JDataRaw,
                         Select_database="IUCN_RedList")


# view all results
ex_results_sdm
ex_results_buff

in_results_sdm
in_results_buff

# reformat results for plotting



ggplot(ex_results_sdm,
       aes(x=GRSex_20, y=ERSex_20,
       )) + #color=temperature_bin)) +
  geom_point(alpha = 0.5, color="red") +
  xlim(0, 100) + ylim(100, 0)


ggplot(ex_results_sdm, 
       aes(x=species, y=ERSex_20)) + #, fill=supp)) +
  geom_bar(stat="identity", position=position_dodge())







##Combine gap analysis metrics
#FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df, FCSin_df = FCSin_df)
#FCSc_mean_df

##Running Conservation indicator across taxa
#indicator_df <- indicator(FCSc_mean_df)
#indicator_df






##Generate summary HTML file with all result
#GetDatasets()

#####Adding test data for running summaryHTML.Rmd separately
#Sl <- speciesList
#Od <- JData
#Rl <- JRasters
#Buffer_distance <- 50000
#Ecoregions_shp <- ecoregions
#Pro_areas <- ProtectedAreas


source("/Users/emily/Documents/GitHub/GapAnalysis/R/SummaryHTML.R")
summaryHTML_file <- SummaryHTML(Species_list=speciesList,
                                Occurrence_data=JData,
                                Raster_list=JRasters,
                                Buffer_distance=50000,
                                Ecoregions_shp=ecoregions,
                                Pro_areas=ProtectedAreas,
                                Output_Folder=file.path(main_dir),
                                writeRasters=FALSE)





##Usage with different buffer distances for ex situ gap analysis
#Buffer distances for 5, 10, and 20 km respectively

buffer_distances <- c(5000,10000,20000)

SRSex_df <- SRSex(Species_list = speciesList,
                  Occurrence_data = CucurbitaData)

FCSex_df_list <- list()


##Running all three ex situ gap analysis steps using FCSex function

##Choose if gap maps are calculated for ex situ gap analysis using diferent buffer size
Gap_Map=FALSE

for(i in 1:length(speciesList)){
  
  FCSex_df_list[[i]] <- FCSex(Species_list=speciesList[i],
                              Occurrence_data=CucurbitaData,
                              Raster_list=CucurbitaRasters[i],
                              Buffer_distance=buffer_distances[i],
                              Ecoregions_shp=ecoregions,
                              Gap_Map=Gap_Map)
  
  
  
};rm(i)

##Returning FCSex object
if(Gap_Map==TRUE){
  FCSex_df <- list(FCSex=do.call(rbind,lapply(FCSex_df_list, `[[`, 1)),
                   GRSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 2)),
                   ERSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 3))
  )
} else {
  FCSex_df <- do.call(rbind,FCSex_df_list)
}


##Running all three in situ gap analysis steps using FCSin function

FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas,
                  Gap_Map = NULL)


##Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)


##Running Conservation indicator across taxa
indicator_df  <- indicator(FCSc_mean_df)
