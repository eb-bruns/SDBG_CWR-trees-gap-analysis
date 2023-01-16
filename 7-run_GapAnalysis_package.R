################################################################################

## 7-run_GapAnalysis_package.R

### Author: Emily Beckman Bruns
### Funding: Cooperative agreement between the United States Botanic Garden and 
# San Diego Botanic Garden (subcontracted to The Morton Arboretum), with 
# support from Botanic Gardens Conservation International U.S.

### Last updated: 10 January 2023
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

my.packages <- c('GapAnalysis',
                 # suggested packages for GapAnalysis workflow
                 'dplyr', 'sp', 'tmap', 'data.table', 'sf', 'methods',
                 'geosphere', 'fasterize', 'rmarkdown', 'knitr', 'rgdal',
                 'rgeos', 'kableExtra', 'DT',
                 # additional packages for mapping and other visualizations
                 'leaflet','RColorBrewer','raster','terra','ggplot2',
                 'forcats'
                 )
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)

# source my updated versions of GapAnalysis functions
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSin.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/SummaryHTML.R")
# source function for rasterizing buffers


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

################################################################################
# Functions
################################################################################

# !!! TESTING !!!
#  Species_list <- speciesList_yesSDM
#  Occurrence_data <- allOcc_yesSDM
#  Raster_list <- RasterList_yesSDM
#  Rasterized_buffers_list <- buff100_list_yesSDM
#  Ecoregions_shp <- ecoregions
#  Pro_areas <- ProtectedAreas
#  Occurrence_data_raw <- allOccRaw_yesSDM
#  Select_database <- "GBIF"
#  reps <- reps_matrix
#  Gap_Map <- F
#  Buffer_distance <- 50000

# function to run gap analysis for multiple buffer sizes and multiple species
# distribution rasters
rep_gap_analysis <- function(Species_list, Occurrence_data, Raster_list,
                             Rasterized_buffers_list,
                             Ecoregions_shp, Pro_areas,
                             Occurrence_data_raw, Select_database,
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
                      Raster_list=Raster_list[[i]],
                      Buffer_distance=reps$GBufferSize[i],
                      Ecoregions_shp=Ecoregions_shp,
                      Occurrence_data_raw=Occurrence_data_raw,
                      Select_database=Select_database)
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
                        Raster_list=Raster_list[[i]],
                        Rasterized_buffers_list=Rasterized_buffers_list,
                        Ecoregions_shp=Ecoregions_shp,
                        Pro_areas=Pro_areas,
                        Occurrence_data_raw=Occurrence_data_raw,
                        Select_database=Select_database)
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

################################################################################
# Read in data for all target taxa
################################################################################

### Protected areas raster
# this seems to download the 10 arc minutes version instead of 2.5 arc min...
#data(ProtectedAreas)
#str(ProtectedAreas)
#res(ProtectedAreas) # resolution = 0.1666667 0.1666667
# ...so I downloaded from the source to read that in instead...
# read in PA file downloaded from Dataverse (https://dataverse.harvard.edu/dataverse/GapAnalysis)
ProtectedAreas <- raster(file.path(main_dir,"gis_layers/wdpa_reclass.tif"))
res(ProtectedAreas) # resolution = 0.04166667 0.04166667

### Global ecoregions shapefile
# this seems to download a version that is cropped to the western US ?...
#data(ecoregions)
#head(ecoregions)
# ...so I downloaded from the source to read that in instead...
# read in as simple features
ecoregions <- sf::read_sf(file.path(
  main_dir,"gis_layers/terr-ecoregions-TNC/tnc_terr_ecoregions.shp"))
# convert to SpatialPolygonsDataFrame
ecoregions <- as_Spatial(ecoregions)

# crop to target regions, to make a little smaller
unique(ecoregions$WWF_REALM2)
ecoregions <- ecoregions[ecoregions$WWF_REALM2 == "Nearctic" | 
                           ecoregions$WWF_REALM2 == "Neotropic",]
# can check visually if desired: plot(ecoregions)

### Occurrence data
# from package example: data(CucurbitaData)
# my occurrence data compiled through full workflow:
# read in raw occurrence data (before keeping only geolocated & removing dups)
# and edit to match format needed for GapAnalysis
all_occ_raw <-  read.csv(file.path(path.pts.raw,
                                  "all_occurrence_data_raw_2023-01-04.csv"))
all_occ_raw <- all_occ_raw %>%
  dplyr::mutate(type = recode(database,
                              "Ex_situ" = "G",
                              .default = "H")) %>%
  dplyr::select(-species) %>%
  dplyr::rename(species = taxon_name_accepted,
                latitude = decimalLatitude,
                longitude = decimalLongitude) %>%
  dplyr::select(species,latitude,longitude,type,database)
all_occ_raw$species <- gsub(" ","_",all_occ_raw$species)
str(all_occ_raw)
# read in flagged occurrence data (ready for spatial analyses)
# and edit to match format needed for GapAnalysis
occ_files <- list.files(path.pts, pattern = ".csv", full.names = T)
occ_dfs <- lapply(occ_files, read.csv)
all_occ <- Reduce(rbind, occ_dfs); rm(occ_files,occ_dfs)
all_occ <- all_occ %>%
  dplyr::mutate(database = recode(database,
                                  "Ex_situ" = "G",
                                  .default = "H")) %>%
  dplyr::filter(
    database == "Ex_situ" | #keep all ex situ points, even if flagged
      # remove some flagged points 
      (.cen & .inst & .outl &
         #.con & .urb & .yr1950 & .yr1980 & .yrna &
         (.nativectry | is.na(.nativectry)) &
         basisOfRecord != "FOSSIL_SPECIMEN" & 
         basisOfRecord != "LIVING_SPECIMEN" &
         establishmentMeans != "INTRODUCED" & 
         establishmentMeans != "MANAGED" &
         establishmentMeans != "CULTIVATED"
      )) %>%
  dplyr::rename(species = taxon_name_accepted,
                latitude = decimalLatitude,
                longitude = decimalLongitude,
                type = database) %>%
  dplyr::select(species,latitude,longitude,type)
all_occ$species <- gsub(" ","_",all_occ$species)
str(all_occ)

### Create aggregated version of ecoregions for clipping buffers to land
pt.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ecoregions_sf <- as(ecoregions, "SpatVector")
ecoregions_proj <- project(ecoregions_sf,pt.proj)
boundary.poly <- aggregate(ecoregions_proj,dissolve = TRUE)

################################################################################
# Set up data for analysis & run, genus by genus
################################################################################

### loop through genus by genus
genera <- c("Asimina","Carya","Castanea","Corylus","Diospyros","Juglans",
            "Malus","Persea","Pistacia","Prunus")
  # df to aggregate data from all genera
all_results_yesSDM <- data.frame()

for(gen_now in 1:length(genera)){
  
  ### Species distribution rasters
  # from package example: data(CucurbitaRasters); CucurbitaRasters <- raster::unstack(CucurbitaRasters)
  # read in own SDM rasters (created from forked CWR-of-the-USA-Gap-Analysis repo)
  sdm_files <- list.files(file.path(path.sdm, genera[gen_now]), 
                          pattern = ".tif", full.names = T)
  sdm_list <- lapply(sdm_files, raster)
  
  ### Create species list from occurrence data
  # from package example: speciesList <- unique(CucurbitaData$species)
  # from own data
  #speciesList_all <- sort(unique(all_occ$species))
  speciesList_yesSDM <- list.files(file.path(path.sdm, genera[gen_now]), 
                                   pattern = ".tif", full.names = F)
  speciesList_yesSDM <- gsub("-spdist_thrsld_median.tif","",speciesList_yesSDM)
  print(paste0("Number of target taxa in ", genera[gen_now], ": ", 
               length(speciesList_yesSDM)))
  #speciesList_noSDM <- speciesList_all[!(speciesList_all %in% speciesList_yesSDM)]
  
  ### Filter occurrence data to just have the target species we want, 
  # depending on if SDM available or not
  allOcc_yesSDM <- all_occ %>% filter(species %in% speciesList_yesSDM)
  allOccRaw_yesSDM <- all_occ_raw %>% filter(species %in% speciesList_yesSDM)
  
  ### Create rastertized buffers around occurrence points, to use as species
  # distribution layer in analyses (next step).
  # do for multiple buffer sizes, here: 20km, 50km, 100km ; can edit as desired
  buff20_list_yesSDM <- list()
  buff50_list_yesSDM <- list()
  buff100_list_yesSDM <- list()
  sp_yesSDM <- split(allOcc_yesSDM, as.factor(allOcc_yesSDM$species))
  for(spp_now in 1:length(sp_yesSDM)){
    # 20km buffers
    spp.buff.20 <- rasterized.buffer(sp_yesSDM[[spp_now]],20000,pt.proj,
                                     pt.proj,boundary.poly,ProtectedAreas)
    buff20_list_yesSDM <- c(buff20_list_yesSDM,spp.buff.20)
    # 50km buffers
    spp.buff.50 <- rasterized.buffer(sp_yesSDM[[spp_now]],50000,pt.proj,
                                     pt.proj,boundary.poly,ProtectedAreas)
    buff50_list_yesSDM <- c(buff50_list_yesSDM,spp.buff.50)
    # 100km buffers
    spp.buff.100 <- rasterized.buffer(sp_yesSDM[[spp_now]],100000,pt.proj,
                                      pt.proj,boundary.poly,ProtectedAreas)
    buff100_list_yesSDM <- c(buff100_list_yesSDM,spp.buff.100)
    print(paste0("Created rasterized buffers for ",
                 sp_yesSDM[[spp_now]]$species[1]))
  }
  #plot(buff50_list_yesSDM[[3]])
  
  ### Set species distribution layer(s) to use in analyses
  # for species with SDM, we will use SDM and three rasterized buffer layers 
  # created above:
  reps_matrix <- data.frame(
    SpeciesDist = c("20km buffers","50km buffers","100km buffers",
                    "SDM"#,"SDM","SDM"
                    ),
    GBufferSize = c(20000,50000,100000,
                    50000#,20000,100000
                    ))
  RasterList_yesSDM <- list(
    buff20_list_yesSDM,buff50_list_yesSDM,buff100_list_yesSDM,
    sdm_list,sdm_list,sdm_list)
  
  ### Run analyses !
  # for species with SDM
  all_results_yesSDM_now <- rep_gap_analysis(
    Species_list = speciesList_yesSDM, 
    Occurrence_data = allOcc_yesSDM, 
    Raster_list = RasterList_yesSDM,
    Rasterized_buffers_list = buff100_list_yesSDM,
    Ecoregions_shp = ecoregions, 
    Pro_areas = ProtectedAreas, 
    Occurrence_data_raw = allOccRaw_yesSDM, 
    Select_database = "GBIF", 
    reps = reps_matrix)
  
  # add genus results to df with all results
  all_results_yesSDM <- rbind(all_results_yesSDM, all_results_yesSDM_now)
  all_results_yesSDM <- all_results_yesSDM %>% distinct()
  
}

tail(all_results_yesSDM)

all_results_yesSDM_ROUND1 <- all_results_yesSDM
write.csv(all_results_yesSDM_ROUND1, file.path(
  main_dir,"RasterizedBuffer-vs-SDM_results","all_results_yesSDM_ALL.csv"))

#all_results_yesSDM <- read.csv(file.path(
#  main_dir,"RasterizedBuffer-vs-SDM_results","all_results_yesSDM_ALL.csv"))

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








all_results_ex <- all_results_yesSDM %>% filter(CalcType=="Ex situ")
all_results_in <- all_results_yesSDM %>% filter(CalcType=="In situ")






SummaryHTML(
  Species_list=speciesList_yesSDM, 
  Occurrence_data=allOcc_yesSDM, 
  Raster_list=RasterList_yesSDM,
  Buffer_distance=50000, 
  Ecoregions_shp=ecoregions,
  Pro_areas=ProtectedAreas, 
  Output_Folder=main_dir, 
  writeRasters=T,
  Occurrence_data_raw=allOccRaw_yesSDM,
  Select_database="GBIF")
  








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
