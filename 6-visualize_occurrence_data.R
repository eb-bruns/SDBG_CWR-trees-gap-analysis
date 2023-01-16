################################################################################

## 6-visualize_occurrence_data.R

### Authors: Emily Beckman Bruns
### Funding: Base script was funded by the Institute of Museum and Library 
# Services (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
# between the United States Botanic Garden and San Diego Botanic Garden
# (subcontracted to The Morton Arboretum), with support from
# Botanic Gardens Conservation International U.S.

### Last updated: 3 January 2023
### R version 4.2.2

### DESCRIPTION:
# Creates interactive (HTML) occurrence point map for each target species,
#   for exploring. Includes toggles that show points flagged in
#   5-refine_occurrence_data.R

### INPUTS:
# target_taxa_with_synonyms.csv
# Occurrence points from 5-refine_occurrence_data.R

### OUTPUTS:
# interactive_maps/visualize_occurrence_data folder with HTML map for each 
#   target taxa (e.g., Quercus_lobata_occurrence_map.html), which can be 
#   downloaded and opened in your browser for exploring

################################################################################
# Load libraries
################################################################################

rm(list=ls())
my.packages <- c("ggplot2","maps","leaflet","RColorBrewer","dplyr","raster",
                 "terra")
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"
  
# or use 0-set_working_directory.R script:
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")
  
################################################################################
# Use leaflet package to create interactive maps to explore (html)
################################################################################

# set up file paths
output.maps <- file.path(main_dir, "interactive_maps","visualize_occurrence_data")
path.gis <- file.path(main_dir,"gis_layers")
path.pts <- file.path(main_dir,"occurrence_data",
                       "standardized_occurrence_data","taxon_edited_points")
path.rasters <- file.path(main_dir,"sdm_rasters")

# select target taxa
taxon_list <- read.csv(file.path(main_dir,taxa_dir,
                                 "target_taxa_with_synonyms.csv"), 
                       header = T, na.strings=c("","NA"),colClasses="character")
  # add country distribution data
taxon_dist <- read.csv(file.path(main_dir,taxa_dir,
  "target_taxa_with_native_dist.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
taxon_list <- left_join(taxon_list,taxon_dist)
  ### UPDATE THIS AS NEEDED:
#no_sdm <- c("Asimina tetramera","Asimina x nashii","Carya x lecontei",
#            "Carya x ludoviciana","Castanea x neglecta",
#            "Juglans microcarpa var. stewartii","Prunus eremophila",
#            "Prunus serotina var. eximia","Prunus x orthosepala")
  # select accepted taxa and remove one that has no occurrence points
target_taxa <- taxon_list %>%
  dplyr::filter(taxon_name_status == "Accepted"
      # optionally, remove species with no SDM (list created manually above)
#        & !(taxon_name_accepted %in% no_sdm)
        )
  nrow(target_taxa) #87 // 96
spp.all <- unique(gsub(" ","_",target_taxa$taxon_name_accepted))
spp.all
  # list of native countries for each target species
countries <- target_taxa$all_native_dist_iso2
  # read in country boundaries shapefile
world_polygons <- vect(file.path(path.gis,
  "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))                             

# create folder for output maps, if not yet created
if(!dir.exists(output.maps)) dir.create(output.maps, recursive=T)

### cycle through each species file and create map
for(i in 1:length(spp.all)){

  # read in records
  spp.now <- read.csv(file.path(path.pts, paste0(gsub("\\.","",spp.all[i]), 
                                                 ".csv")))
  
  #target.iso <- unlist(strsplit(countries[i],split="; "))
  #target_countries <- world_polygons[world_polygons$IOS %in% target.iso,]

  ## palette based on database
  # set database as factor and order appropriately
  spp.now$database <- factor(spp.now$database,
    levels = c("Ex_situ","GBIF","NorthAm_herbaria","iDigBio",
               "IUCN_RedList","FIA","BIEN"))
  spp.now <- spp.now %>% arrange(desc(database))
  # create color palette
  colors <- c("#adbb3f","#819756","#5fbb9a","#6a9ebd","#7b83cc","#7264de",
              "#3c2c7a")
  database.pal <- colorFactor(palette=colors,
    levels = c("Ex_situ","GBIF","NorthAm_herbaria","iDigBio",
               "IUCN_RedList","FIA","BIEN"))

  ## read in species distribution model and select color for mapping
    ## this is the old code for using SDMs from Khoury et al 2020 (PNAS):
    #genus <- strsplit(spp.now$taxon_name_accepted[1]," ")[[1]][1]
    #taxon <- gsub(" ","%20",spp.now$taxon_name_accepted[1])
    #raster_path <- paste0(
    #  "https://raw.githubusercontent.com/dcarver1/cwr_pnas_results/main/speciesLevelData20210706/",
    #  genus,"/",taxon,"/",taxon,"__thrsld_median.tif")
    #download.file(raster_path,
    #  destfile=file.path(path.rasters,paste0(spp.all[i], "_PNAS_2020_SDM.tif")))
#  spp.raster <- raster(file.path(
#    path.rasters,paste0(spp.all[i],"-spdist_thrsld_median.tif")))
  # reclassification of 0 to NA, so we don't view the ecoregions layer
  # this was not needed for the 2020 PNAS SDMs
#  spp.raster[spp.raster==0] <- NA
  # set color palette
#  raster.pal <- colorNumeric("#e0aee6",values(spp.raster),na.color = "transparent")
  
  # create map
  try(final_map <- leaflet() %>%
    # Base layer groups
    #addProviderTiles(providers$CartoDB.PositronNoLabels,
    #  group = "CartoDB.PositronNoLabels") %>%
    addProviderTiles(providers$CartoDB.Positron,
      group = "CartoDB.Positron") %>%
    addControl(paste0("<b>",spp.all[i]), position = "topright") %>%
    addControl(
      "Toggle the checkboxes below on/off to view flagged points (colored red) in each category.</br>
      If no points turn red when box is checked, there are no points flagged in that category.</br>
      Click each point for more information about the record.",
      position = "topright") %>%
    # SDM
#    addRasterImage(spp.raster,colors=raster.pal,opacity = 0.8) %>%
	  # Native country outlines
	  #addPolygons(data = target_countries, fillColor = "transparent",
		#  weight = 2, opacity = 0.8, color = "#7a7a7a") %>%
    # Color by database
    addCircleMarkers(data = spp.now, ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      color = ~database.pal(database),radius = 4,
      fillOpacity = 0.9, stroke = T) %>%
    # Overlay groups (can toggle)
    addCircleMarkers(data = spp.now %>% filter(!.cen & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Within 500m of country/state centroid (.cen)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.urb & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "In urban area (.urb)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.inst & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Within 100m of biodiversity institution (.inst)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.con & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Not in reported country (.con)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.outl & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Geographic outlier (.outl)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.nativectry),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Outside native country (.nativectry)") %>%
    addCircleMarkers(data = spp.now %>%
      filter(basisOfRecord == "FOSSIL_SPECIMEN" |
        basisOfRecord == "LIVING_SPECIMEN"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "FOSSIL_SPECIMEN or LIVING_SPECIMEN (basisOfRecord)") %>%
    addCircleMarkers(data = spp.now %>%
      filter(establishmentMeans == "INTRODUCED" |
        establishmentMeans == "MANAGED" |
        establishmentMeans == "INVASIVE" |
        establishmentMeans == "CULTIVATED"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "INTRODUCED, MANAGED, CULTIVATED, or INVASIVE (establishmentMeans)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.yr1950 & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Recorded prior to 1950 (.yr1950)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.yr1980 & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Recorded prior to 1980 (.yr1980)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.yrna & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Year:</b> ",year,"<br/>",
        "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
        "<b>Dataset name:</b> ",datasetName,"<br/>",
        "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      radius=4,stroke=T,color="black",weight=1,fillColor="red",fillOpacity=0.8,
      group = "Year unknown (.yrna)") %>%
    # Layers control
    addLayersControl(
      #baseGroups = c("CartoDB.PositronNoLabels",
      #               "CartoDB.Positron",
      #               "Esri.WorldTopoMap",
      #               "Stamen.Watercolor"),
      overlayGroups = c("Within 500m of country/state centroid (.cen)",
                        "In urban area (.urb)",
                        "Within 100m of biodiversity institution (.inst)",
                        "Not in reported country (.con)",
                        "Geographic outlier (.outl)",
                        "Outside native country (.nativectry)",
                        "FOSSIL_SPECIMEN or LIVING_SPECIMEN (basisOfRecord)",
                        "INTRODUCED, MANAGED, CULTIVATED, or INVASIVE (establishmentMeans)",
                        "Recorded prior to 1950 (.yr1950)",
                        "Recorded prior to 1980 (.yr1980)",
                        "Year unknown (.yrna)"),
      options = layersControlOptions(collapsed = FALSE)) %>%
    ## uncomment these to hide the group (unchecked in controls on map)
    #hideGroup("Within 500m of country/state centroid (.cen)") %>%
    hideGroup("In urban area (.urb)") %>%
    #hideGroup("Within 100m of biodiversity institution (.inst)") %>%
    hideGroup("Not in reported country (.con)") %>%
    #hideGroup("Geographic outlier (.outl)") %>%
    #hideGroup("Outside native country (.nativectry)") %>%
    #hideGroup("FOSSIL_SPECIMEN or LIVING_SPECIMEN (basisOfRecord)") %>%
    #hideGroup("INTRODUCED, MANAGED, CULTIVATED, or INVASIVE (establishmentMeans)") %>%
    hideGroup("Recorded prior to 1950 (.yr1950)") %>%
    hideGroup("Recorded prior to 1980 (.yr1980)") %>%
    hideGroup("Year unknown (.yrna)") %>%
#    addLegend(labels = c("Present","Absent"),colors = c("#e0aee6","white"),
#      title = "Species Distribution Model", position = "bottomright",
#      opacity = 0.8) %>%
    addLegend(pal = database.pal, values = unique(spp.now$database),
      title = "Occurrence point</br>source database", position = "bottomright", opacity = 0.8) %>%
    addControl(
      "See https://github.com/eb-bruns/SDBG_CWR-trees-gap-analysis
      for information about data sources and flagging methodology.",
      position = "bottomleft"))
  final_map

  # save map
  try(htmlwidgets::saveWidget(final_map, file.path(output.maps,
    paste0(spp.all[i], "_occurrence_map.html"))))

  cat("\tEnding ", spp.all[i], ", ", i, " of ", length(spp.all), ".\n\n", sep="")
}
