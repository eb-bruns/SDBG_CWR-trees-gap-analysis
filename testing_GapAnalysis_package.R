
################################################################################
##Load package
#install.packages("GapAnalysis")
library(GapAnalysis)

##Load additional suggested packages
my.packages <- c('dplyr', 'sp', 'tmap', 'data.table', 'sf', 'methods',
  'geosphere', 'fasterize', 'rmarkdown', 'knitr', 'rgdal',
  'rgeos', 'kableExtra', 'DT',
  #EBB: additional packages for leaflet mapping and other visualizations:
  'leaflet','RColorBrewer','raster','terra','ggplot2'
  )
lapply(my.packages, require, character.only=TRUE)

## assign main working directory
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")
## set up file paths
path.pts.raw <- file.path(main_dir,occ_dir,"raw_occurrence_data")
path.pts <- file.path(main_dir,occ_dir,"standardized_occurrence_data",
                      "taxon_edited_points")
path.sdm <- file.path(main_dir,"sdm_rasters")

##Obtaining protected areas raster
# this seems to download the 10 arc minutes version instead of 2.5 arc min
#data(ProtectedAreas)
#str(ProtectedAreas)
#res(ProtectedAreas) # resolution = 0.1666667 0.1666667
# so I downloaded from the source to read that in instead...
# read in PA file from Dataverse:
#   https://dataverse.harvard.edu/dataverse/GapAnalysis
ProtectedAreas <- raster(file.path(
  main_dir,"gis_layers/wdpa_reclass.tif"))
res(ProtectedAreas) # resolution = 0.04166667 0.04166667

##Obtaining ecoregions shapefile
# this seems to download a version that is cropped to the western US ?
#data(ecoregions)
#head(ecoregions)
# so I downloaded from the source to read that in instead...
# read in ecoregions file downloaded from the source
ecoregions <- sf::st_read(file.path(
  main_dir,"gis_layers/terr-ecoregions-TNC/tnc_terr_ecoregions.shp"))
# crop to target regions, to make a little smaller
ecoregions <- ecoregions[ecoregions$WWF_REALM2 == "Nearctic" | 
                         ecoregions$WWF_REALM2 == "Neotropic",]
#plot(ecoregions)



### TESTING WITH JUGLANS CALIFORNICA, J. HINDSII & J. MAJOR

##Obtaining occurrences from example
#data(CucurbitaData)
#str(CucurbitaData)
  # read in raw occurrence data (before keeping only geolocated & removing dups)
  # and edit to match format needed for GapAnalysis
OccDataRaw <-  read.csv(file.path(path.pts.raw,
                                  "all_occurrence_data_raw_2023-01-04.csv"))
OccDataRaw <- OccDataRaw %>%
  dplyr::mutate(type = recode(database,
                                "Ex_situ" = "G",
                                .default = "H")) %>%
  dplyr::select(-species) %>%
  dplyr::rename(species = taxon_name_accepted,
                latitude = decimalLatitude,
                longitude = decimalLongitude) %>%
  dplyr::select(species,latitude,longitude,type,database)
OccDataRaw$species <- gsub(" ","_",OccDataRaw$species)
str(OccDataRaw)
  # read in filtered occurrence data (ready for spatial analyses)
  # and edit to match format needed for GapAnalysis
occ_files <- list.files(path.pts, pattern = ".csv", full.names = T)
occ_dfs <- lapply(occ_files, read.csv)
OccData <- Reduce(rbind, occ_dfs); rm(occ_files,occ_dfs)
OccData <- OccData %>%
  dplyr::mutate(database = recode(database,
                                    "Ex_situ" = "G",
                                    .default = "H")) %>%
  dplyr::filter(
    database == "Ex_situ" | #keep all ex situ points, even if flagged
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
OccData$species <- gsub(" ","_",OccData$species)
str(OccData)

##Obtaining raster_list
#data(CucurbitaRasters)
#CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#str(CucurbitaRasters)
# read in SDM output rasters and stack
sdm_files <- list.files(path.sdm, pattern = ".tif", full.names = T)
sdm_list <- lapply(sdm_files, raster)
length(sdm_list) #90

##Obtaining species names from the data
#speciesList <- unique(CucurbitaData$species)
# separate list for species with SDM and without
speciesList_all <- sort(unique(OccData$species))
speciesList_yesSDM <- list.files(path.sdm, pattern = ".tif", full.names = F)
speciesList_yesSDM <- gsub("-spdist_thrsld_median.tif","",speciesList_yesSDM)
speciesList_noSDM <- speciesList_all[!(speciesList_all %in% speciesList_yesSDM)]

##EBB: Testing creating a rasterized buffer layer to use instead of SDM
# function for creating aggregated buffers around occ points then rasterizing
rasterized.buffer <- function(df,radius,pt_proj,buff_proj,boundary,template_raster){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c("longitude", "latitude"),
                   crs=pt_proj)
  # reproject to specified projection
  proj_df <- project(spat_pts,buff_proj)
  # place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # make into simple feature object type
  buffers_clip_sf <- sf::st_as_sf(buffers_clip)
  # rasterize the buffers
  rasterized_buffer <- fasterize::fasterize(buffers_clip_sf, template_raster)
  # trim raster extent to area with values
  rasterized_buffer_trim <- raster::trim(rasterized_buffer,padding=5)
  #raster_blueprint <- raster()
  #extent(raster_blueprint) <- extent(buffers_clip_sf)
  #res(raster_blueprint) <- res(pa_raster)
  #rasterized_buffer <- rasterize(buffers_clip_sf, raster_blueprint)
  # return raster
  return(rasterized_buffer_trim)
}
# define projections
pt.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#calc.proj <- "+proj=aea + lat_1=29.5 + lat_2=45.5 + lat_0=37.5 + lon_0=-96 +x_0=0 +y_0=0 + ellps =GRS80 +datum=NAD83 + units=m +no_defs"
# create boundary for clipping buffers to land
ecoregions_sf <- as(ecoregions, "SpatVector")
ecoregions_proj <- project(ecoregions_sf,pt.proj)
boundary.poly <- aggregate(ecoregions_proj,dissolve = TRUE)

# create rasertized buffers for target taxa
  # do for multiple buffer sizes: 20km, 50km, 100km
buff20_list <- list(); buff50_list <- list(); buff100_list <- list()
sp_split <- split(OccData, as.factor(OccData$species))
for(i in 1:length(sp_split)){
  # 20km buffers
  spp.buff.20 <- rasterized.buffer(sp_split[[i]],20000,pt.proj,pt.proj,
                                   boundary.poly,ProtectedAreas)
  buff20_list <- c(buff20_list,spp.buff.20)
  # 50km buffers
  spp.buff.50 <- rasterized.buffer(sp_split[[i]],50000,pt.proj,pt.proj,
                                   boundary.poly,ProtectedAreas)
  buff50_list <- c(buff50_list,spp.buff.50)
  # 100km buffers
  spp.buff.100 <- rasterized.buffer(sp_split[[i]],100000,pt.proj,pt.proj,
                                   boundary.poly,ProtectedAreas)
  buff100_list <- c(buff100_list,spp.buff.100)
  print(paste0("Ending ",sp_split[[i]]$species[1]))
}
plot(buff50_list[[3]])


# function for running gap analysis for multiple buffer sizes and multiple
#   species distribution rasters
rep_gap_analysis <- function(Species_list, Occurrence_data, Raster_list,
                             Rasterized_buffers_list,
                             Ecoregions_shp, Pro_areas,
                             Occurrence_data_raw, Select_database,
                             buffer_sizes, spdist_type){
  all_results <- data.frame()
  for(i in 1:length(spdist_type)){
    # run all ex situ analyses
    ex_results <- data.frame()
    for(j in 1:length(buffer_sizes)){
        FCSex_df <- FCSex(Species_list=Species_list,
                          Occurrence_data=Occurrence_data,
                          Raster_list=Raster_list[[i]],
                          Buffer_distance=buffer_sizes[j],
                          Ecoregions_shp=Ecoregions_shp,
                          Occurrence_data_raw=Occurrence_data_raw,
                          Select_database=Select_database)
        FCSex_df$CalcType <- "Ex situ"
        FCSex_df$SpeciesDist <- spdist_type[i] 
        FCSex_df$BufferSizeKm <- (buffer_sizes[j]/1000)
        # rename columns to remove ex (we have column for that instead)
        colnames(FCSex_df)[2:6] <- c("SRS","GRS","ERS","FCS","FCS_class")
        # rename columns to include buffer size used (this was for wide df)
        #colnames(FCSex_df)[2:6] <- c(
        #  paste0("SRSex",spdist_type[i]),
        #  paste0("GRSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
        #  paste0("ERSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
        #  paste0("FCSex",buffer_sizes[j]/1000,"_",spdist_type[i]),
        #  paste0("FCSex_class",buffer_sizes[j]/1000,"_",spdist_type[i]))
        # add to all results dataframe
        ex_results <- rbind(ex_results,FCSex_df)
    }
    # run all in situ analyses
    in_results <- FCSin(Species_list=Species_list,
                        Occurrence_data=Occurrence_data,
                        Raster_list=Raster_list[[i]],
                        Rasterized_buffers_list=Rasterized_buffers_list,
                        Ecoregions_shp=Ecoregions_shp,
                        Pro_areas=Pro_areas,
                        Occurrence_data_raw=Occurrence_data_raw,
                        Select_database=Select_database)
    in_results$CalcType <- "In situ"
    in_results$SpeciesDist <- spdist_type[i] 
    in_results$BufferSizeKm <- "N/A"
    # rename columns to remove ex (we have column for that instead)
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

source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSex.R")
source("/Users/emily/Documents/GitHub/GapAnalysis/R/FCSin.R")

bufferSize <- c(20000,50000,100000)

# for species with SDM
spdistType <- c("SDM","A) 20km buffers","B) 50km buffers","C) 100km buffers")
RasterList <- list(sdm_list,buff20_list,buff50_list,buff100_list)
  
ga_all_results <- rep_gap_analysis(
  Species_list=speciesList, Occurrence_data=JData, Raster_list=RasterList, 
  Rasterized_buffers_list=JDataBuff100,
  Ecoregions_shp=ecoregions, Pro_areas=ProtectedAreas, 
  Occurrence_data_raw=JDataRaw, Select_database="IUCN_RedList", 
  buffer_sizes=bufferSize, spdist_type=spdistType)

# 

all_results_ex <- ga_all_results %>% filter(CalcType=="Ex situ")
all_results_in <- ga_all_results %>% filter(CalcType=="In situ")

# all SRS results
ggplot(
  ga_all_results, aes(x = SRS, y = species)) + 
  geom_point(color="black") +
  xlim(0,100) +
  facet_grid(cols = vars(CalcType))
# all FCS results
ggplot(
  ga_all_results, aes(x = FCS, y = SpeciesDist)) + 
  geom_point(color="black") +
  xlim(0,100) +
  facet_grid(rows = vars(species), cols = vars(CalcType))
# all GRS results
ggplot(
  ga_all_results, aes(x = GRS, y = CalcType)) + 
  geom_point(color="purple") +
  xlim(0,100) +
  facet_grid(rows = vars(species), cols = vars(SpeciesDist))
# all ERS results
ggplot(
  ga_all_results, aes(x = ERS, y = species)) + 
  geom_point(color="blue") +
  xlim(0,100) +
  facet_grid(rows = vars(SpeciesDist), cols = vars(BufferSizeKm))








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
