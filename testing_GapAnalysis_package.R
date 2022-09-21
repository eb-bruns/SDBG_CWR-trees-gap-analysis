

##Load package
#install.packages("GapAnalysis")
library(GapAnalysis)

##Load additional suggested packages
my.packages <- c('dplyr', 'sp', 'tmap', 'data.table', 'sf', 'methods', 'geosphere', 'data.table',
'fasterize', 'rmarkdown', 'knitr', 'rgdal', 'rgeos', 'kableExtra', 'DT')
lapply(my.packages, require, character.only=TRUE)

## assign main working directory
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis"

### TESTING WITH JUGLANS CALIFORNICA & J. HINDSII

##Obtaining occurrences from example
#data(CucurbitaData)
#str(CucurbitaData)
  # read in data
JcalData <- read.csv(file.path(
  main_dir,"In situ - H - records/outputs/taxon_edited_points/Juglans_californica.csv"))
JhinData <- read.csv(file.path(
  main_dir,"In situ - H - records/outputs/taxon_edited_points/Juglans_hindsii.csv"))
JData <- rbind(JcalData,JhinData)
str(JData)
  # edit to match format needed for GapAnalysis
JData <- JData %>%
  mutate(database = recode(database,
                           "Ex_situ" = "G",
                           .default = "H")) %>%
  rename(species = taxon_name_acc,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         type = database) %>%
  select(species,latitude,longitude,type)
JData$species <- gsub(" ","_",JData$species)
str(JData)

##Obtaining species names from the data
#speciesList <- unique(CucurbitaData$species)
speciesList <- unique(JData$species)
speciesList

##Obtaining raster_list
#data(CucurbitaRasters)
#CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#str(CucurbitaRasters)
  # read in SDM output rasters and stack
JcalRaster <- raster(file.path(
  main_dir,"In situ - H - records/inputs/PNAS_2020_SDMs/Juglans_californica_PNAS_2020_SDM.tif"))
JhinRaster <- raster(file.path(
  main_dir,"In situ - H - records/inputs/PNAS_2020_SDMs/Juglans_hindsii_PNAS_2020_SDM.tif"))
JRasters <- list(JcalRaster,JhinRaster)
str(JRasters)

##Obtaining protected areas raster
data(ProtectedAreas)
str(ProtectedAreas)

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
str(FCSex_df)

##Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=JData,
                  Raster_list=JRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)
str(FCSin_df)

##Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)
str(FCSc_mean_df)

##Running Conservation indicator across taxa
indicator_df <- indicator(FCSc_mean_df)
str(indicator_df)

##Generate summary HTML file with all result
GetDatasets()
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
