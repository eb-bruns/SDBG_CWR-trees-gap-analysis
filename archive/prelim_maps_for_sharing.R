################################################################################

## prelim_maps_for_sharing.R

### Authors: Emily Beckman Bruns
### Funding:
# Base script was funded by the Institude of Museum and Library Services
#   (IMLS MFA program grant MA-30-18-0273-18 to The Morton Arboretum).
# Significant updates/additions with funding from a cooperative agreement
#   between the United States Botanic Garden and San Diego Botanic Garden
#   (subcontracted to The Morton Arboretum), with support from
#   Botanic Gardens Conservation International U.S.

### Creation date: 18 October 2022
### Last updated: 18 October 2022
### R version 4.2.2

### DESCRIPTION:
  # Creates interactive (HTML) map for each target taxon

### DATA IN:
  # Occurrence points from 5-refine_occurrence_points.R
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

################################################################################
# Functions
################################################################################

source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/mapping_functions.R")

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/Gap-Analysis-Mapping"

# or use 0-set_working-directory.R script:
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")

# set up file paths
path.pts <<- file.path(main_dir,occ_dir,"standardized_occurrence_data","taxon_edited_points")
path.sdm <<- file.path(main_dir,gis_dir,"PNAS_2020_SDMs")
path.out.figs <<- file.path(main_dir,"interactive_maps")

################################################################################
################################################################################
# Use leaflet package to create interactive maps to explore (html)
################################################################################

# define projections we'll use throughout
#		points will be WGS84
#pt.proj <- "+proj=longlat +datum=WGS84"
pt.proj <<- "EPSG:4326"
#   for calculations we need something in meters, like Equal Earth Projection
# calc.proj <- "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
calc.proj <<- "+proj=aea + lat_1=29.5 + lat_2=45.5 + lat_0=37.5 + lon_0=-96 +x_0=0 +y_0=0 + ellps =GRS80 +datum=NAD83 + units=m +no_defs"

# read in / load polygon data
  # ecoregions
ecoregions <- vect(file.path(main_dir,gis_dir,"terr-ecoregions-TNC", 
                             "tnc_terr_ecoregions.shp"))
# country boundaries
world_countries <- vect(file.path(main_dir,gis_dir,
  "UIA_World_Countries_Boundaries/World_Countries__Generalized_.shp"))
  # create subset with only target countries
#target_iso <- c("US","MX","CA")
#target_countries_shp <- subset(world_countries,
#  world_countries$ISO %in% target_iso,)
target_countries_shp <<- world_countries
  # create polygon for clipping points later (project to pt projection)
#target_countries.pt <- project(target_countries_shp,pt.proj)
#boundary.pt <- aggregate(target_countries.pt,dissolve = TRUE)
  # create polygon for clipping buffers later
ecoregions_proj <- project(ecoregions,pt.proj)
boundary.poly <<- aggregate(ecoregions_proj,dissolve = TRUE)
  # create clipped version of ecoregions
eco_clip <- crop(ecoregions_proj,target_countries_shp)
eco_clip_sf <<- sf::st_as_sf(eco_clip)

# select target taxa
taxon_list <- read.csv(file.path(main_dir,taxa_dir,"target_taxa_with_synonyms.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
  # add country distribution data
taxon_dist <- read.csv(file.path(main_dir,taxa_dir,
  "target_taxa_with_native_dist.csv"), header = T, na.strings=c("","NA"),
  colClasses="character")
taxon_list <- left_join(taxon_list,taxon_dist)
target_taxa <- taxon_list %>%
  filter(taxon_name_status == "Accepted")
  # list of native countries for each target species
countries <- target_taxa$all_native_dist_iso2

  # create subsets of taxa that have data
spp.all <- list.files(path = path.pts, pattern = ".csv", full.names = F)
spp.all <<- unique(gsub(".csv","",spp.all))
yes_sdm <- list.files(path = path.sdm, pattern = ".tif", full.names = F)
yes_sdm <<- unique(gsub("_PNAS_2020_SDM.tif","",yes_sdm))

# recreate for a few target species
#spp.all <- c("Diospyros_texana","Prunus_andersonii","Prunus_fremontii",
#             "Prunus_fasciculata","Carya_aquatica","Carya_illinoinensis",
#             "Carya_laciniosa","Carya_myristiciformis","Carya_ovalis",
#             "Carya_texana","Juglans_californica","Prunus_eremophila",
#             "Prunus_fasciculata_var_punctata")

##
### cycle through each species file and create map
##
for(i in 1:length(spp.all)){

  cat("starting ",spp.all[i],", ",i," of ",length(spp.all),".\n", sep="")
  
  ## read in occurrence records (includes ex situ)
  spp.pts <- read.csv(file.path(path.pts, paste0(spp.all[i], ".csv")))

  ## filter records based on flagging columns from 5-refine_occurrence_data.R
  spp.pts <- spp.pts %>%
    mutate(
      .cen = as.logical(.cen),
      #.urb = as.logical(.urb),
      .inst = as.logical(.inst),
      #.con = as.logical(.con),
      .outl = as.logical(.outl),
      #.gtsnative = as.logical(.gtsnative),
      #.rlnative = as.logical(.rlnative),
      .yr1950 = as.logical(.yr1950),
      .yr1980 = as.logical(.yr1980),
      .yrna = as.logical(.yrna)
    ) %>%
    # select or deselect these filters as desired:
    filter(
      database == "Ex_situ" |
      (.cen & .inst & .outl &
         #.con & 
         #.urb & .yr1950 & .yr1980 & .yrna &
         #(.gtsnative | is.na(.gtsnative)) &
         #(.rlnative  | is.na(.rlnative)) &
         #(.rlintroduced | is.na(.rlintroduced)) &
         basisOfRecord != "FOSSIL_SPECIMEN" & 
         basisOfRecord != "LIVING_SPECIMEN" &
         establishmentMeans != "INTRODUCED" & 
         establishmentMeans != "MANAGED" &
         establishmentMeans != "CULTIVATED" &
         latlong_countryCode %in% c("US","CA","MX")))

  ## create layer of buffers around points
  spp.pts.buffer <- create.buffers(spp.pts,50000,pt.proj,pt.proj,boundary.poly)

  ## create subset of ecoregions within buffers (takes too long)
  #eco.sel <- intersect(spp.pts.buffer,eco_clip_sf)
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
      map <- create.full.map(eco_clip_sf,spp.raster,raster.pal,spp.pts.buffer,spp.pts,
                             ex.pts.buffer,ex.pts)
    } else { # no sdm
      ## map without sdm
      print("mapping without SDM")
      map <- create.nosdm.map(eco_clip_sf,spp.pts.buffer,spp.pts,ex.pts.buffer,ex.pts)
    }
  } else { # no ex situ
    if(spp.all[i] %in% yes_sdm){
      ## read in raster from Khoury et al 2020 (PNAS)
      spp.raster <- raster(file.path(path.sdm,paste0(spp.all[i],"_PNAS_2020_SDM.tif")))
      # select color for raster when mapped
      raster.pal <- colorNumeric("#0e6187",values(spp.raster),na.color = "transparent")
      ## map without ex situ
      print("mapping without ex situ")
      map <- create.noex.map(eco_clip_sf,spp.raster,raster.pal,spp.pts.buffer,spp.pts)
    } else { # no ex situ and no sdm
      ## map without ex situ or sdm
      print("mapping without SDM & ex situ")
      map <- create.noex.nosdm.map(eco_clip_sf,spp.pts.buffer,spp.pts)
    }
  }

  # save map
  htmlwidgets::saveWidget(map, file.path(path.out.figs,
    paste0(spp.all[i], "__prelim_map_v2.html")))

}
