################################################################################

### GapAnalysis_additions.R
### Authors: Emily Beckman Bruns
### Funding:

### Creation date: 17 November 2022
### Last full check and update:
### R version 4.2.2

### DESCRIPTION:
  # This script uses the GapAnalysis package as a base for calculating
  # spatial in situ (protected areas) and ex situ (botanical gardens and
  # genebanks) representation of a species for conservation.
  # We test using buffers around occurrence points versus a Species
  # Distribution Model (SDM)

### INPUTS:
  #

### OUTPUTS:
  #

################################################################################
# Load libraries
################################################################################

my.packages <- c('GapAnalysis', 'dplyr', 'sp', 'tmap', 'data.table', 'sf',
  'methods', 'geosphere', 'data.table','fasterize', 'rmarkdown', 'knitr',
  'rgdal', 'rgeos', 'kableExtra', 'DT')
# install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping/"
# or use 0-1_set_workingdirectory.R script:
source("./Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-1_set_working_directory.R")

# set up file paths we'll use in this script
path.pts <- file.path(main_dir,"OLD-occurrence_points","OUTPUTS_FROM_R","taxon_edited_points")
path.sdm <- file.path(main_dir,"gis_layers","PNAS_2020_SDMs")

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
  # remove a couple target taxa that have no occurrence data
  taxon_list <- taxon_list %>%
    filter(taxon_name_accepted != "Carya x ludoviciana")

# list of target taxon names (accepted)
taxon_names <- sort(unique(taxon_list$taxon_name_accepted))
length(taxon_names) #95

## IF YOU HAVE JUST ONE OR A FEW TARGET TAXA, CREATE A LIST BY HAND:
#taxon_names <- c("Quercus havardii")
#length(taxon_names)

################################################################################
# Read in data
################################################################################

## occurrence point data

# compiled using 2-0_get_occurrence_data.R & 3-0_compile_occurrence_data.R
# read in data for each species and combine into one dataframe
all_occ <- data.frame()
for(i in 1:length(taxon_names)){
  taxon_nm <- gsub(" ","_",taxon_names[i])
  taxon_occ <- read.csv(file.path(path.pts,paste0(taxon_nm,".csv")))
  all_occ <- rbind(all_occ, taxon_occ)
}
head(all_occ)

# edit data to match format needed for GapAnalysis
all_occ <- all_occ %>%
  dplyr::mutate(database = recode(database,
                           "Ex_situ" = "G",
                           .default = "H")) %>%
  dplyr::rename(species = taxon_name_acc,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         type = database) %>%
  dplyr::select(species,latitude,longitude,type)
all_occ$species <- gsub(" ","_",all_occ$species)
str(all_occ)

## species list input for GapAnalysis
speciesList <- unique(all_occ$species)
speciesList

##Obtaining raster_list
#data(CucurbitaRasters)
#CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#str(CucurbitaRasters)
  # read in SDM output rasters and stack
JcalRaster <- raster(file.path(main_dir,path.sdm,"Juglans_californica_PNAS_2020_SDM.tif"))
JhinRaster <- raster(file.path(main_dir,path.sdm,"Juglans_hindsii_PNAS_2020_SDM.tif"))
JmajRaster <- raster(file.path(main_dir,path.sdm,"Juglans_major_PNAS_2020_SDM.tif"))
JRasters <- list(JmajRaster,JcalRaster,JhinRaster)
str(JRasters)

##Obtaining protected areas raster
  # this seems to download the 10 arc minutes version instead of 2.5 arc min
data(ProtectedAreas)
str(ProtectedAreas)
  res(ProtectedAreas) # resolution = 0.1666667 0.1666667
  # so I downloaded from the source to read that in instead...
# read in PA file from Dataverse:
#   https://dataverse.harvard.edu/dataverse/GapAnalysis
ProtectedAreas <- raster(file.path(
  main_dir,"gis_layers/wdpa_reclass.tif"))
  res(ProtectedAreas) # resolution = 0.04166667 0.04166667

##Obtaining ecoregions shapefile
data(ecoregions)
head(ecoregions)

##Running all three ex situ gap analysis steps using FCSex function
FCSex_df <- FCSex(Species_list=speciesList,
                  Occurrence_data=JData,
                  Raster_list=JRasters,
                  Buffer_distance=50000,
                  Ecoregions_shp=ecoregions
)
FCSex_df

##Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=JData,
                  Raster_list=JRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)
FCSin_df

##Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)
FCSc_mean_df

##Running Conservation indicator across taxa
indicator_df <- indicator(FCSc_mean_df)
indicator_df




##EBB: Testing creating a rasterized buffer layer to use instead of SDM
library("terra")
# function for creating aggregated buffers around occurrence points
create.buffers <- function(df,radius,pt_proj,buff_proj,boundary){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c("decimalLongitude", "decimalLatitude"),
                   crs=pt_proj)
  # reproject to specified projection
  proj_df <- project(spat_pts,buff_proj)
  # place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # make into object that can be mapped in leaflet
  buffers_clip_sf <- sf::st_as_sf(buffers_clip)
  # return buffer polygons
  return(buffers_clip_sf)
}
# define projections
pt.proj <- "+proj=longlat +datum=WGS84"
calc.proj <- "+proj=aea + lat_1=29.5 + lat_2=45.5 + lat_0=37.5 + lon_0=-96 +x_0=0 +y_0=0 + ellps =GRS80 +datum=NAD83 + units=m +no_defs"
# create boundary for clipping buffers to land
ecoregions_sf <- as(ecoregions, "SpatVector")
ecoregions_proj <- project(ecoregions_sf,pt.proj)
boundary.poly <- aggregate(ecoregions_proj,dissolve = TRUE)
## read in occurrence records (includes ex situ)
#spp.pts <- read.csv(file.path(path.pts, paste0(spp.all[i], ".csv")))
## create layer of buffers around points
spp.pts.buffer <- create.buffers(JcalData,50000,pt.proj,pt.proj,boundary.poly)
## rasterize the vector
raster_blueprint <- raster()
extent(raster_blueprint) <- extent(spp.pts.buffer)
res(rasterized_buffer) <- res(ProtectedAreas)
rasterized_buffer <- rasterize(spp.pts.buffer, raster_blueprint)
plot(rasterized_buffer)

##Run the gap analysis again with this new layer
FCSex_df_buff <- FCSex(Species_list=speciesList[2],
                  Occurrence_data=JData,
                  Raster_list=rasterized_buffer,
                  Buffer_distance=50000,
                  Ecoregions_shp=ecoregions
)
FCSex_df_buff
FCSin_df_buff <- FCSin(Species_list=speciesList[2],
                  Occurrence_data=JData,
                  Raster_list=rasterized_buffer,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)
FCSin_df_buff





##Generate summary HTML file with all result
GetDatasets()

  #####Adding test data for running summaryHTML.Rmd separately
  Sl <- speciesList
  Od <- JData
  Rl <- JRasters
  Buffer_distance <- 50000
  Ecoregions_shp <- ecoregions
  Pro_areas <- ProtectedAreas

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
