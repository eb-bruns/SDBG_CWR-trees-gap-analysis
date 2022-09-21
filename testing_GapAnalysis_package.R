

##Load package
#install.packages("GapAnalysis")
library(GapAnalysis)

##Load additional suggested packages
my.packages <- c('sp', 'tmap', 'data.table', 'sf', 'methods', 'geosphere', 'data.table',
'fasterize', 'rmarkdown', 'knitr', 'rgdal', 'rgeos', 'kableExtra', 'DT')
lapply(my.packages, require, character.only=TRUE)

##Obtaining occurrences from example
data(CucurbitaData)
str(CucurbitaData)

##Obtaining species names from the data
speciesList <- unique(CucurbitaData$species)

##Obtaining raster_list
data(CucurbitaRasters)
CucurbitaRasters <- raster::unstack(CucurbitaRasters)
str(CucurbitaRasters)

##Obtaining protected areas raster
data(ProtectedAreas)
str(ProtectedAreas)

##Obtaining ecoregions shapefile
data(ecoregions)
head(ecoregions)

#Running all three ex situ gap analysis steps using FCSex function
FCSex_df <- FCSex(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Buffer_distance=50000,
                  Ecoregions_shp=ecoregions
)
str(FCSex_df)

#Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)
str(FCSin_df)

## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)
str(FCSc_mean_df)

##Running Conservation indicator across taxa
indicator_df <- indicator(FCSc_mean_df)
str(indicator_df)

## Generate summary HTML file with all result
GetDatasets()
summaryHTML_file <- SummaryHTML(Species_list=speciesList,
                                Occurrence_data = CucurbitaData,
                                Raster_list=CucurbitaRasters,
                                Buffer_distance=50000,
                                Ecoregions_shp=ecoregions,
                                Pro_areas=ProtectedAreas,
                                Output_Folder=".",
                                writeRasters=FALSE)



##Usage with different buffer distances for ex situ gap analysis
#Buffer distances for 5, 10, and 20 km respectively

buffer_distances <- c(5000,10000,20000)

SRSex_df <- SRSex(Species_list = speciesList,
                  Occurrence_data = CucurbitaData)

FCSex_df_list <- list()


#Running all three ex situ gap analysis steps using FCSex function

#Choose if gap maps are calculated for ex situ gap analysis using diferent buffer size
Gap_Map=FALSE

for(i in 1:length(speciesList)){

  FCSex_df_list[[i]] <- FCSex(Species_list=speciesList[i],
                    Occurrence_data=CucurbitaData,
                    Raster_list=CucurbitaRasters[i],
                    Buffer_distance=buffer_distances[i],
                    Ecoregions_shp=ecoregions,
                    Gap_Map=Gap_Map)



};rm(i)

#Returning FCSex object
if(Gap_Map==TRUE){
  FCSex_df <- list(FCSex=do.call(rbind,lapply(FCSex_df_list, `[[`, 1)),
                 GRSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 2)),
                 ERSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 3))
                 )
} else {
  FCSex_df <- do.call(rbind,FCSex_df_list)
}


#Running all three in situ gap analysis steps using FCSin function

FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas,
                  Gap_Map = NULL)


## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)


##Running Conservation indicator across taxa
indicator_df  <- indicator(FCSc_mean_df)
