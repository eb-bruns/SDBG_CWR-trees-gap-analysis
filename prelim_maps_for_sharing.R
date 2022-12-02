################################################################################

## 4-1_prelim_maps_for_sharing.R

### Authors: Emily Beckman Bruns
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Updates/additions with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 18 October 2022
### Last updated: 18 October 2022

### R version 4.1.3

### DESCRIPTION:
  # Creates interactive (HTML) map for each target species,

### DATA IN:
  # Occurrence points from 3-1_refine_occurrence_points.R
  # Global ecoregions:
  # 	Terrestrial Ecoregions of the World, via WWF (Olson et al., 2001)
  #		https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
  # Protected areas?

### DATA OUT:
  # spp_interactive_maps folder with HTML map for each target species
  #   (e.g., Quercus_lobata_prelim_map.html), which can be downloaded and opened
  #   in your browser for exploring

################################################################################
# Load libraries
################################################################################

rm(list=ls())
my.packages <- c("leaflet","RColorBrewer","dplyr","raster","terra",
  "rnaturalearth","rnaturalearthhires","Polychrome") #"ggplot2","maps",
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

filter <- dplyr::filter
select <- dplyr::select

################################################################################
# Functions
################################################################################

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

# create map
create.full.map <- function(eco_clip,spp.raster,raster.pal,spp.pts.buffer,spp.pts,
                       ex.pts.buffer,ex.pts){
  map <- leaflet() %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
      group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    #addPolygons(data = eco_clip,
    #            fillColor = "white", fillOpacity = 0,
    #            color = "#b8bab6", weight = 1, opacity = 0.7,
    #            popup = eco_clip$ECO_NAME,
    #            group = "Ecoregions (Olson et al. 2001)") %>%
    # SDM from PNAS 2020
    addRasterImage(spp.raster, colors=raster.pal, opacity = 0.8,
      group = "Species distribution model (Khoury et al. 2020)") %>%
    # In situ buffers
  	addPolygons(data = spp.pts.buffer,
  		fillColor = "#965da7", fillOpacity = 0.45,
  		weight = 0, opacity = 0, color = "white",
  		smoothFactor = 0,
      group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
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
      color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)") %>%
    # Ex situ buffers
  	addPolygons(data = ex.pts.buffer,
  		fillColor = "#80bf30", fillOpacity = 0.7,
  		weight = 0, color = "white", opacity = 0,
  		smoothFactor = 0,
      group = "Estimated area represented by ex situ material (50km buffers around source localities)") %>%
    # Ex situ points
    addCircleMarkers(data = ex.pts, ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Collection date:</b> ",year,"<br/>",
        "<b>Provenance type:</b> ",basisOfRecord,"<br/>",
        "<b>Institution name:</b> ",datasetName,"<br/>",
        "<b>Number of individuals:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      color = "#335707", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Source localities for ex situ material (in botanic gardens/genebanks)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
      position = "topright") %>%
  	addLegend(labels =
  		c("Species distribution model (Khoury et al. 2020)",
  			"Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)"),
  		colors = c("#0e6187","#56116b","#965da7","#335707","#80bf30"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Note that occurrence points are verbatim from each source (not filtered)",
      position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(#"Ecoregions (Olson et al. 2001)",
                        "Species distribution model (Khoury et al. 2020)",
                        "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
                        "Estimated species distribution (50km buffers around occurrence points)",
                        "Source localities for ex situ material (in botanic gardens/genebanks)",
                        "Estimated area represented by ex situ material (50km buffers around source localities)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
    #hideGroup("Ecoregions (Olson et al. 2001)") %>%
    addControl(
      "Toggle the checkboxes on and off to control which layers are visible...",
      position = "bottomright") %>%
    addControl(
      "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      position = "bottomleft") %>%
    setView(-98, 40, zoom = 5)
  return(map)
}

# create map
create.noex.map <- function(eco_clip,spp.raster,raster.pal,spp.pts.buffer,spp.pts){
  map <- leaflet() %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
      group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    #addPolygons(data = eco_clip,
    #            fillColor = "white", fillOpacity = 0,
    #            color = "#b8bab6", weight = 1, opacity = 0.7,
    #            popup = eco_clip$ECO_NAME,
    #            group = "Ecoregions (Olson et al. 2001)") %>%
    # SDM from PNAS 2020
    addRasterImage(spp.raster, colors=raster.pal, opacity = 0.8,
      group = "Species distribution model (Khoury et al. 2020)") %>%
    # In situ buffers
  	addPolygons(data = spp.pts.buffer,
  		fillColor = "#965da7", fillOpacity = 0.45,
  		weight = 0, opacity = 0, color = "white",
  		smoothFactor = 0,
      group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
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
      color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
      position = "topright") %>%
  	addLegend(labels =
  		c("Species distribution model (Khoury et al. 2020)",
  			"Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)"),
  		colors = c("#0e6187","#56116b","#965da7","#335707","#80bf30"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Note that occurrence points are verbatim from each source (not filtered)",
      position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(#"Ecoregions (Olson et al. 2001)",
                        "Species distribution model (Khoury et al. 2020)",
                        "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
                        "Estimated species distribution (50km buffers around occurrence points)",
                        "Source localities for ex situ material (in botanic gardens/genebanks)",
                        "Estimated area represented by ex situ material (50km buffers around source localities)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
    #hideGroup("Ecoregions (Olson et al. 2001)") %>%
    addControl(
      "Toggle the checkboxes on and off to control which layers are visible...",
      position = "bottomright") %>%
    addControl(
      "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      position = "bottomleft") %>%
    setView(-98, 40, zoom = 5)
  return(map)
}

# create map
create.noex.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts){
  map <- leaflet() %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
      group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    #addPolygons(data = eco_clip,
    #            fillColor = "white", fillOpacity = 0,
    #            color = "#b8bab6", weight = 1, opacity = 0.7,
    #            popup = eco_clip$ECO_NAME,
    #            group = "Ecoregions (Olson et al. 2001)") %>%
    # In situ buffers
  	addPolygons(data = spp.pts.buffer,
  		fillColor = "#965da7", fillOpacity = 0.45,
  		weight = 0, opacity = 0, color = "white",
  		smoothFactor = 0,
      group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
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
      color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
      position = "topright") %>%
  	addLegend(labels =
  		c("Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)"),
  		colors = c("#56116b","#965da7","#335707","#80bf30"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Note that occurrence points are verbatim from each source (not filtered)",
      position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(#"Ecoregions (Olson et al. 2001)",
                        "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
                        "Estimated species distribution (50km buffers around occurrence points)",
                        "Source localities for ex situ material (in botanic gardens/genebanks)",
                        "Estimated area represented by ex situ material (50km buffers around source localities)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
    #hideGroup("Ecoregions (Olson et al. 2001)") %>%
    addControl(
      "Toggle the checkboxes on and off to control which layers are visible...",
      position = "bottomright") %>%
    addControl(
      "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      position = "bottomleft") %>%
    setView(-98, 40, zoom = 5)
  return(map)
}

# create map
create.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts,ex.pts.buffer,ex.pts){
  map <- leaflet() %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
      group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    #addPolygons(data = eco_clip,
    #            fillColor = "white", fillOpacity = 0,
    #            color = "#b8bab6", weight = 1, opacity = 0.7,
    #            popup = eco_clip$ECO_NAME,
    #            group = "Ecoregions (Olson et al. 2001)") %>%
    # In situ buffers
  	addPolygons(data = spp.pts.buffer,
  		fillColor = "#965da7", fillOpacity = 0.45,
  		weight = 0, opacity = 0, color = "white",
  		smoothFactor = 0,
      group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
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
      color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)") %>%
    # Ex situ buffers
  	addPolygons(data = ex.pts.buffer,
  		fillColor = "#80bf30", fillOpacity = 0.7,
  		weight = 0, color = "white", opacity = 0,
  		smoothFactor = 0,
      group = "Estimated area represented by ex situ material (50km buffers around source localities)") %>%
    # Ex situ points
    addCircleMarkers(data = ex.pts, ~decimalLongitude, ~decimalLatitude,
      popup = ~paste0(
        "<b>Accepted species name:</b> ",taxon_name_acc,"<br/>",
        "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
        "<b>Source database:</b> ",database,"<br/>",
        "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
        "<b>Collection date:</b> ",year,"<br/>",
        "<b>Provenance type:</b> ",basisOfRecord,"<br/>",
        "<b>Institution name:</b> ",datasetName,"<br/>",
        "<b>Number of individuals:</b> ",establishmentMeans,"<br/>",
        "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
        "<b>ID:</b> ",UID),
      color = "#335707", radius = 3, fillOpacity = 0.9, #stroke = T,
      group = "Source localities for ex situ material (in botanic gardens/genebanks)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
      position = "topright") %>%
  	addLegend(labels =
  		c("Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)"),
  		colors = c("#56116b","#965da7","#335707","#80bf30"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Note that occurrence points are verbatim from each source (not filtered)",
      position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(#"Ecoregions (Olson et al. 2001)",
                        "Occurrence points (GBIF, IUCN Red List, iDigBio, herbaria consortia)",
                        "Estimated species distribution (50km buffers around occurrence points)",
                        "Source localities for ex situ material (in botanic gardens/genebanks)",
                        "Estimated area represented by ex situ material (50km buffers around source localities)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
    #hideGroup("Ecoregions (Olson et al. 2001)") %>%
    addControl(
      "Toggle the checkboxes on and off to control which layers are visible...",
      position = "bottomright") %>%
    addControl(
      "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      position = "bottomleft") %>%
    setView(-98, 40, zoom = 5)
  return(map)
}

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/occurrence_points"
#script_dir <- "./Documents/GitHub/OccurrencePoints/scripts"
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"
local_dir <- "/Users/emily/Desktop/*work"

# or use 0-1_set_workingdirectory.R script:
#source("./Documents/GitHub/OccurrencePoints/scripts/0-1_set_workingdirectory.R")
#source("scripts/0-1_set_workingdirectory.R")


################################################################################
################################################################################
# Use leaflet package to create interactive maps to explore (html)
################################################################################

# define projections we'll use throughout
#		points will be WGS84
pt.proj <- "+proj=longlat +datum=WGS84"
#   for calculations we need something in meters, like Equal Earth Projection
# calc.proj <- "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
calc.proj <- "+proj=aea + lat_1=29.5 + lat_2=45.5 + lat_0=37.5 + lon_0=-96 +x_0=0 +y_0=0 + ellps =GRS80 +datum=NAD83 + units=m +no_defs"

# read in / load polygon data
  # ecoregions
ecoregions <- vect(file.path(main_dir,"gis_layers","global-ecoregions",
  "wwf_terr_ecos.shp"))
  # country boundaries
world_countries <- vect(file.path(main_dir,"gis_layers",
  "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))
  # create subset with only target countries
target_iso <- c("US","MX","CA")
target_countries_shp <- subset(world_countries,
  world_countries$ISO %in% target_iso,)
  # create polygon for clipping points later (project to pt projection)
#target_countries.pt <- project(target_countries_shp,pt.proj)
#boundary.pt <- aggregate(target_countries.pt,dissolve = TRUE)
  # create polygon for clipping buffers later
ecoregions_proj <- project(ecoregions,pt.proj)
boundary.poly <- aggregate(ecoregions_proj,dissolve = TRUE)
  # create clipped version of ecoregions
eco_clip <- crop(ecoregions_proj,target_countries_shp)
eco_clip <- sf::st_as_sf(eco_clip)

# set up file paths
path.pts <- file.path(main_dir,"occurrence_points","OUTPUTS_FROM_R","taxon_edited_points")
path.sdm <- file.path(main_dir,"gis_layers","PNAS_2020_SDMs")
path.eco <- file.path(main_dir,"gis_layers","global-ecoregions")
path.out.figs <- file.path(local_dir,"interactive_maps")

# create folder for output maps, if not yet created
if(!dir.exists(path.out.figs)) dir.create(path.outlfigs, recursive=T)

# select target taxa
taxon_list <- read.csv(file.path(main_dir,"target_taxa_with_synonyms.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
  # add country distribution data
taxon_dist <- read.csv(file.path(main_dir,"taxa_metadata",
  "target_taxa_with_native_dist.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
taxon_list <- left_join(taxon_list,taxon_dist)
target_taxa <- taxon_list %>%
  filter(taxon_name_status == "Accepted")
  # list of native countries for each target species
countries <- target_taxa$all_native_dist_iso2

  # create subsets of taxa that have data
spp.all <- list.files(path = path.pts, pattern = ".csv", full.names = F)
  spp.all <- unique(gsub(".csv","",spp.all))
yes_sdm <- list.files(path = path.sdm, pattern = ".tif", full.names = F)
  yes_sdm <- unique(gsub("_PNAS_2020_SDM.tif","",yes_sdm))

##
### cycle through each species file and create map
##
for(i in 1:length(spp.all)){

  ## read in occurrence records (includes ex situ)
  spp.pts <- read.csv(file.path(path.pts, paste0(spp.all[i], ".csv")))

  ## create layer of buffers around points
  spp.pts.buffer <- create.buffers(spp.pts,50000,pt.proj,pt.proj,boundary.poly)

  ## create subset of ecoregions within buffers (takes too long)
  #eco.sel <- intersect(spp.pts.buffer,eco_clip)
  #  # create ecoregion color palette
  #  eco.pal <- colorFactor(palette = "Greys", domain = eco.sel$ECO_ID,
  #  	reverse = F, na.color = "white")

  ## subset with just ex situ points
  ex.pts <- spp.pts %>% filter(database == "Ex_situ")
  # if there are ex situ points
  if(nrow(ex.pts)>0){
    ## create layer of buffers around ex situ points
    ex.pts.buffer <- create.buffers(ex.pts,50000,pt.proj,pt.proj,boundary.poly)
        # split ex situ data by number of individuals, to use different symbols
        #exsitu1 <- exsitu %>% arrange(establishmentMeans) %>%
        #  filter(establishmentMeans <= few_indiv)
        #exsitu2 <- exsitu %>% arrange(establishmentMeans) %>%
        #  filter(establishmentMeans > few_indiv & establishmentMeans < many_indiv)
        #exsitu3 <- exsitu %>% arrange(establishmentMeans) %>%
        #  filter(establishmentMeans >= few_indiv)
    # if there is a 2020 SDM
    if(spp.all[i] %in% yes_sdm){
      ## read in raster from Khoury et al 2020 (PNAS)
      spp.raster <- raster(file.path(path.sdm,paste0(spp.all[i],"_PNAS_2020_SDM.tif")))
      # select color for raster when mapped
      raster.pal <- colorNumeric("#0e6187",values(spp.raster),na.color = "transparent")
      ## map everything
      print("mapping everything")
      map <- create.full.map(eco_clip,spp.raster,raster.pal,spp.pts.buffer,spp.pts,
                             ex.pts.buffer,ex.pts)
    } else { # no sdm
      ## map without sdm
      print("mapping without SDM")
      map <- create.nosdm.map(eco_clip,spp.pts.buffer,spp.pts,ex.pts.buffer,ex.pts)
    }
  } else { # no ex situ
    if(spp.all[i] %in% yes_sdm){
      ## read in raster from Khoury et al 2020 (PNAS)
      spp.raster <- raster(file.path(path.sdm,paste0(spp.all[i],"_PNAS_2020_SDM.tif")))
      # select color for raster when mapped
      raster.pal <- colorNumeric("#0e6187",values(spp.raster),na.color = "transparent")
      ## map without ex situ
      print("mapping without ex situ")
      map <- create.noex.map(eco_clip,spp.raster,raster.pal,spp.pts.buffer,spp.pts)
    } else { # no ex situ and no sdm
      ## map without ex situ or sdm
      print("mapping without SDM & ex situ")
      map <- create.noex.nosdm.map(eco_clip,spp.pts.buffer,spp.pts)
    }
  }

  # save map
  htmlwidgets::saveWidget(map, file.path(path.out.figs,
    paste0(spp.all[i], "__prelim_map.html")))

  cat("\tEnding ", spp.all[i], ", ", i, " of ", length(spp.all), ".\n\n", sep="")
}
