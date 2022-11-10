################################################################################

## 4-0_plot_occurrence_points_CWRcopy.R

### Authors: Emily Beckman Bruns & Christy Rollinson
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Moderate edits were added with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 03 June 2022
### Last updated: 06 June 2022

### R version 4.1.3

### DESCRIPTION:
  # Creates interactive (HTML) occurrence point map for each target species,
  #   for exploring. Includes toggles that show points flagged in
  #   3-1_refine_occurrence_points.R
  #   Also creates two fixed basic (PNG) maps for each target species: one with
  #   all valid occurrence points (output from 3-0_compile_occurrence_points.R)
  #   and another with all flagged points removed (output from
  #   3-1_refine_occurrence_points.R)

### DATA IN:
  # Occurrence points from 3-1_refine_occurrence_points.R

### DATA OUT:
  # spp_interactive_maps folder with HTML map for each target species
  #   (e.g., Quercus_lobata_leafet_map.html), which can be downloaded and opened
  #   in your browser for exploring

################################################################################
# Load libraries
################################################################################

rm(list=ls())
my.packages <- c("ggplot2","maps","leaflet","RColorBrewer","dplyr","raster")
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/occurrence_points"
#script_dir <- "./Documents/GitHub/OccurrencePoints/scripts"
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/In situ - H - records"

# or use 0-1_set_workingdirectory.R script:
#source("./Documents/GitHub/OccurrencePoints/scripts/0-1_set_workingdirectory.R")
#source("scripts/0-1_set_workingdirectory.R")


################################################################################
################################################################################
# Use leaflet package to create interactive maps to explore (html)
################################################################################

# set up file paths
output <- file.path(main_dir, "outputs")
path.pts <- file.path(output, "taxon_edited_points")
path.figs <- file.path(output, "taxon_interactive_maps")
path.rasters <- file.path(main_dir,"inputs","PNAS_2020_SDMs")

# select target taxa
taxon_list <- read.csv(file.path(main_dir,"inputs",
  "target_taxa_with_syn.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
  # add country distribution data
taxon_dist <- read.csv(file.path(main_dir,"inputs","known_distribution",
  "target_taxa_with_native_dist.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
taxon_list <- left_join(taxon_list,taxon_dist)
no_sdm <- c("Asimina incana","Asimina longifolia","Asimina obovata",
            "Asimina parviflora","Asimina pygmaea","Asimina reticulata",
            "Asimina tetramera","Juglans jamaicensis","Persea palustris",
            "Prunus serotina") # this one throws an error because it is too large?
  # select accepted taxa and remove one that has no occurrence points
target_taxa <- taxon_list %>%
  dplyr::filter(taxon_name_status == "Accepted" &
                taxon_name_acc != "Prunus +orthosepala" &
      # optionally, remove species with no SDM from 2020 anlysis
                !grepl("\\+",taxon_name_acc) &
                !(taxon_name_acc %in% no_sdm))
  nrow(target_taxa) #78
spp.all <- unique(gsub(" ","_",target_taxa$taxon_name_acc))
spp.all
  # list of native countries for each target species
countries <- target_taxa$all_native_dist_iso2
  # load polygon data
load(file.path(poly_dir, "inputs", "gis_data", "admin_shapefiles.RData"))

# create folder for output maps, if not yet created
if(!dir.exists(path.figs)) dir.create(path.figs, recursive=T)

# create folder for rasters, if not yet created
if(!dir.exists(path.rasters)) dir.create(path.rasters, recursive=T)

### cycle through each species file and create map
for(i in 67:length(spp.all)){

  # read in records
  spp.now <- read.csv(file.path(path.pts, paste0(spp.all[i], ".csv")))

  target.iso <- unlist(strsplit(countries[i],split="; "))
  target_countries <- adm0.poly[adm0.poly@data$country.iso_a2 %in% target.iso,]

  ## palette based on database
  # set database as factor and order appropriately
  spp.now$database <- factor(spp.now$database,
    levels = c("Ex_situ",#"FIA",
               "GBIF","US_Herbaria","iDigBio",
               #"BISON","BIEN",
               "IUCN_RedList"))
  spp.now <- spp.now %>% arrange(desc(database))
  # create color palette
  colors <- c("#188562","#147053","#115e46","#0d4735","#093326")
  database.pal <- colorFactor(palette=colors,
    levels = c("Ex_situ",#"FIA",
               "GBIF","US_Herbaria","iDigBio",
               #"BISON","BIEN",
               "IUCN_RedList"))

  ## read in raster from Khoury et al 2020 (PNAS)
  genus <- strsplit(spp.now$taxon_name_acc[1]," ")[[1]][1]
  taxon <- gsub(" ","%20",spp.now$taxon_name_acc[1])
  raster_path <- paste0(
    "https://raw.githubusercontent.com/dcarver1/cwr_pnas_results/main/speciesLevelData20210706/",
    genus,"/",taxon,"/",taxon,"__thrsld_median.tif")
  download.file(raster_path,
    destfile=file.path(path.rasters,paste0(spp.all[i], "_PNAS_2020_SDM.tif")))
  spp.raster <- raster(file.path(path.rasters,paste0(spp.all[i], "_PNAS_2020_SDM.tif")))
    # select color for raster when mapped
  raster.pal <- colorNumeric("#c2a763",values(spp.raster),na.color = "transparent")

  # create map
    map <- leaflet() %>%
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
    # SDM from PNAS 2020
    addRasterImage(spp.raster,colors=raster.pal,opacity = 0.8) %>%
	  # Native country outlines
	  addPolygons(data = target_countries, fillColor = "transparent",
		  weight = 2, opacity = 0.8, color = "#7a7a7a") %>%
    # Color by database
    addCircleMarkers(data = spp.now, ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
    addCircleMarkers(data = spp.now %>% filter(!.gtsnative),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
      group = "Outside GTS native country (.gtsnative)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.rlnative),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
      group = "Outside IUCN RL native country (.rlnative)") %>%
    addCircleMarkers(data = spp.now %>%
      filter(basisOfRecord == "FOSSIL_SPECIMEN" |
        basisOfRecord == "LIVING_SPECIMEN"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        establishmentMeans == "INVASIVE"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
      group = "INTRODUCED, MANAGED, or INVASIVE (establishmentMeans)") %>%
    addCircleMarkers(data = spp.now %>% filter(!.yr1950 & database!="Ex_situ"),
      ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
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
                        "Outside GTS native country (.gtsnative)",
                        "Outside IUCN RL native country (.rlnative)",
                        "FOSSIL_SPECIMEN or LIVING_SPECIMEN (basisOfRecord)",
                        "INTRODUCED, MANAGED, or INVASIVE (establishmentMeans)",
                        "Recorded prior to 1950 (.yr1950)",
                        "Recorded prior to 1980 (.yr1980)",
                        "Year unknown (.yrna)"),
      options = layersControlOptions(collapsed = FALSE)) %>%
    #hideGroup("Within 500m of country/state centroid (.cen)") %>%
    hideGroup("In urban area (.urb)") %>%
    #hideGroup("Within 100m of biodiversity institution (.inst)") %>%
    #hideGroup("Not in reported country (.con)") %>%
    #hideGroup("Geographic outlier (.outl)") %>%
    #hideGroup("Outside GTS native country (.gtsnative)") %>%
    #hideGroup("Outside IUCN RL native country (.rlnative)") %>%
    hideGroup("In IUCN RL introduced country (.rlintroduced)") %>%
    #hideGroup("FOSSIL_SPECIMEN or LIVING_SPECIMEN (basisOfRecord)") %>%
    #hideGroup("INTRODUCED, MANAGED, or INVASIVE (establishmentMeans)") %>%
    hideGroup("Recorded prior to 1950 (.yr1950)") %>%
    hideGroup("Recorded prior to 1980 (.yr1980)") %>%
    hideGroup("Year unknown (.yrna)") %>%
    addLegend(labels = c("Present","Absent"),colors = c("#c2a763","white"),
      title = "PNAS 2020 SDM", position = "bottomright",
      opacity = 0.8) %>%
    addLegend(pal = database.pal, values = unique(spp.now$database),
      title = "Occurrence point</br>source database", position = "bottomright", opacity = 0.8) %>%
    addControl(
      "See https://github.com/MortonArb-CollectionsValue/OccurrencePoints
      for information about data sources and flagging methodology.",
      position = "bottomleft")
  map

  # save map
  htmlwidgets::saveWidget(map, file.path(path.figs,
    paste0(spp.all[i], "_leaflet_map.html")))

  cat("\tEnding ", spp.all[i], ", ", i, " of ", length(spp.all), ".\n\n", sep="")
}
