################################################################################

## 2-2_map_exsitu_versus_SDM.R

### Author: Emily Beckman Bruns
### Funding:
# Cooperative agreement between the United States Botanic Garden and
#  San Diego Botanic Garden (subcontracted to The Morton Arboretum),
#  with support from Botanic Gardens Conservation International U.S.

### Creation date: 10 June 2022
### Last updated: 10 June 2022

### R version 4.1.3

### DESCRIPTION:
  #

### DATA IN:
  #

### DATA OUT:
  #

################################################################################
# Load libraries
################################################################################

#rm(list=ls())
my.packages <- c('plyr','tidyverse','leaflet','rnaturalearth','rnaturalearthdata','raster')
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

select <- dplyr::select

################################################################################
# Set working directory
################################################################################

# Google Drive folder with input ex situ & raster data
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/CWR North America Gap Analysis/In situ - H - records"
# Google Drive folder with polygon data for mapping
#poly_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"

################################################################################
# Load functions
################################################################################

# map individual species with raster data and ex situ points + buffers
map_sp <- function(wild_dist,state_bound,gb_pts,gb_color,bg_pts,bg_color){
    # create palette
  pal <- colorNumeric(c("#329940"), values(wild_dist), na.color = "transparent") #57945f
    # get CRS
  wgs.proj <- wild_dist@crs
    # create map
  map <- leaflet() %>%
  	## background
  	addProviderTiles("CartoDB.PositronNoLabels",#"Esri.WorldGrayCanvas","Esri.WorldShadedRelief","Esri.WorldTerrain",
  		options = providerTileOptions(maxZoom = 10)) %>%
    ## taxon name label
    addControl(paste0("<b>",bg_pts$taxon_name_acc[1]), position = "topright") %>%
  	## EPA Level IV ecoregions
  	#addPolygons(
  	#	data = ecol4_sel,
  	#	fillColor = ~ecol4_pal(ecol4_sel@data$US_L4NAME), fillOpacity = 0.6,
  	#	color = ~ecol4_pal(ecol4_sel@data$US_L4NAME), weight = 0.5, opacity = 1) %>%
  	## in situ distribution (raster)
    addRasterImage(wild_dist, colors = pal, opacity = 0.7) %>%
  	## state boundaries
  	addPolygons(data = state_bound,
  		fillOpacity = 0, color = "#969696", weight = 1.5, opacity = 1) %>%
    #
    ##### POLYGONS
    #
    ### GENEBANK DATA
    ## ex situ buffers
  	addPolygons(data = create.buffers(gb_pts,50000,wgs.proj,wgs.proj),
  		smoothFactor = 0.5,	weight = 1.5, opacity = 1, color = gb_color,
  		fillOpacity = 0.3) %>%
    ### BOTANIC GARDEN DATA
    ## ex situ buffers
  	addPolygons(data = create.buffers(bg_pts,50000,wgs.proj,wgs.proj),
  		smoothFactor = 0.5,	weight = 1.5, opacity = 1, color = bg_color,
  		fillOpacity = 0.3) %>%
    #
    ##### POINTS
    #
    ### GENEBANK DATA
  	## ex situ points
  	addCircleMarkers(data = gb_pts,
  		lng = ~long_dd, lat = ~lat_dd,
  		popup = ~UID,
  		radius = 4, fillOpacity = 1, stroke = F, color = gb_color) %>%
    ### BOTANIC GARDEN DATA
  	## ex situ points
  	addCircleMarkers(data = bg_pts,
  		lng = ~long_dd, lat = ~lat_dd,
  		popup = ~UID,
  		radius = 4, fillOpacity = 1, stroke = F, color = bg_color) %>%
    #
    ## title
  	#addControl("In situ distribution and representation in ex situ collections"),
  	#	position = "topright") %>%
  	## legend for points
  	addLegend(labels =
  		c("Genebank holds the germplasm","Botanic garden holds the germplasm"),
  		colors = c(gb_color,bg_color),
      title = "Wild collection sites, plus 50km buffers",
  		position = "bottomleft", opacity = 0.8) %>%
  	## legend for raster
  	addLegend(labels = "",
  		colors = "#6ca874", title = "Predicted wild distribution (PNAS 2020)",
  		position = "bottomleft", opacity = 0.8) %>%
  	## scale bar
  	addScaleBar(position = "bottomright",
  		options = scaleBarOptions(maxWidth = 150)) %>%
  	## pick coords for center of frame and zoom level when map first appears
  	setView(-118, 35, zoom = 7)
  return(map)
}

################################################################################
################################################################################
# Maps for presentations
################################################################################

##### SUMMER 2022 PGOC MEETING

## read in ex situ point data

# ex situ data from 2020 PNAS analysis (genebank only)
#pnas_exsitu <- read.csv(file.path(main_dir,"inputs","PNAS_2020_Exsitu",
#  "Juglans2020-07-30.csv"), header = T)
#pnas_exsitu <- pnas_exsitu %>%
#  rename(taxon_name_acc = taxon,
#         lat_dd = latitude,
#         long_dd = longitude) %>%
#  filter(type == "G" & lat_dd != "NULL")
#pnas_exsitu$lat_dd <- as.numeric(pnas_exsitu$lat_dd)
#pnas_exsitu$long_dd <- as.numeric(pnas_exsitu$long_dd)
#nrow(pnas_exsitu) #76

# read in ex situ data from 2022 analysis (BG and genebank)
exsitu <- read.csv(file.path(main_dir,"inputs","Ex situ",
  "ExSitu_Compiled_Post-Geolocation_2022-06-07.csv"),
  header = T)

# read in target taxon list
taxon_list <- read.csv(file.path(main_dir, "inputs", "target_taxa_with_syn.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
str(taxon_list) #266

# number of institutions who responded to survey
length(sort(unique(exsitu$inst_short[which(!grepl("[0-9]",exsitu$inst_short))]))) #193
  # number in the US/Canada
NA_inst <- exsitu %>% distinct(inst_short,inst_country) %>%
  filter(!grepl("[0-9]",inst_short)) %>%
  filter(grepl("US|CA",inst_country)) %>%
  filter(inst_short != "NatlPlantGermplasmSystem")
  nrow(NA_inst) #106
# total number (including genebanks)
length(sort(unique(exsitu$inst_short))) #231

# map institutions
inst_exsitu <- exsitu %>%
  select(inst_short,inst_lat,inst_long) %>% distinct() %>%
  filter(inst_short != "NatlPlantGermplasmSystem")

  map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = inst_exsitu,
      lng = ~inst_long, lat = ~inst_lat,
      radius = 9,color = "purple",fillOpacity = 0.6,stroke = F,
      popup = ~inst_short) #%>%
    #addLegend(
    #  labels = "Botanic gardens that provided ex situ accessions data",
    #  colors = "purple",position = "bottomright",opacity = 1)
  map

## accession-level stats
  # by taxon and provenance type
exsitu_acc <- exsitu %>%
  mutate(prov_type = recode(prov_type,
    "N" = "Wild",
    "H?" = "Hort./Unknown",
    "W" = "Wild",
    "Z" = "Wild",
    "NG" = "Hort./Unknown",
    "H" = "Hort./Unknown",
    "U" = "Hort./Unknown")) %>%
  group_by(taxon_name_acc,prov_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
  # add taxon list, to see which have no accessions
exsitu_names <- data.frame(taxon_name_acc = sort(unique(taxon_list$taxon_name_acc)))
exsitu_names
exsitu_acc <- Reduce(full_join,list(exsitu_acc,exsitu_names))
exsitu_acc$num_records[is.na(exsitu_acc$num_records)] <- 0
exsitu_acc$prov_type[is.na(exsitu_acc$prov_type)] <- "Wild"
exsitu_acc$taxon_name_acc <- gsub("\\+","x ",exsitu_acc$taxon_name_acc)
  # plot
exsitu_acc %>%
  mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=prov_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.8)),
      axis.title.x = element_text(size = rel(2)),
      axis.title.y = element_text(size = rel(2)),
      legend.text = element_text(size = rel(1.8)),
      legend.title = element_text(size = rel(2))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target taxa (North American fruit & nut tree CWR)",
         fill = "Provenance Type")
  #ggsave(file.path(main_dir,"outputs","exsitu_prov_barchart.png")#,
  #width = NA, height = NA)
  # look at top species
top_sp <- exsitu_acc %>%
  group_by(taxon_name_acc) %>%
  summarize(num_records = sum(num_records)) %>%
  arrange(desc(num_records))
top_sp <- as.data.frame(top_sp)[1:5,1]
for(i in 1:length(top_sp)){
  print(top_sp[i])
  top <- exsitu %>%
    filter(taxon_name_acc == top_sp[i])
  print(nrow(top))
  top <- top %>%
    group_by(inst_short) %>% count() %>% ungroup() %>%
    rename(num_records = n) %>% arrange(desc(num_records))
  print(head(as.data.frame(top)))
}
## create table to compare to PNAS results
pnas <- exsitu %>%
  group_by(taxon_name_acc) %>%
  dplyr::count(coord_na = is.na(lat_dd)) %>%
  ungroup()
pnas <- rbind.fill(pnas,as.data.frame(unique(taxon_list$taxon_name_acc)))
pnas_wide <- pnas %>%
  spread("coord_na","n")
as.data.frame(pnas_wide)
## by taxon and institution type
exsitu_inst <- exsitu %>%
  mutate(inst_type = recode(inst_type,
    "Gene/Seed Bank" = "Genebank",
    "Botanic Garden; Gene/Seed Bank" = "Unknown")) %>%
  group_by(taxon_name_acc,inst_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
exsitu_names
exsitu_inst <- Reduce(full_join,list(exsitu_inst,exsitu_names))
exsitu_inst$num_records[is.na(exsitu_inst$num_records)] <- 0
exsitu_inst$inst_type[is.na(exsitu_inst$inst_type)] <- "Unknown"
exsitu_inst$taxon_name_acc <- gsub("\\+","x ",exsitu_inst$taxon_name_acc)
exsitu_inst <- exsitu_inst %>% separate("taxon_name_acc","genus",sep=" ",remove=F)
  # plot
exsitu_inst %>%
  #mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=inst_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    #facet_wrap(vars(genus), scales="free_y") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target taxa (North American fruit & nut tree CWR)",
         fill = "Institution Type")
## map all points
exsitu_pts <- exsitu %>%
  distinct(taxon_name_acc,lat_dd,long_dd,.keep_all=T) %>%
  filter(UID != "MEX208~1031700021~NG~Prunus serotina")
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = exsitu_pts,
      lng = ~long_dd, lat = ~lat_dd,
      radius = 6,color = "green",fillOpacity = 0.6,stroke = F,
      popup = exsitu_pts$UID) %>%
    addControl("North American exsitu tree CWR: All wild collection locations",
      position = "topright")
  map

## map specific species as examples

  # keep only records with lat-long
exsitu_ll <- exsitu %>%
  filter(!is.na(lat_dd))
  # separate genebank and bg data
bg_exsitu <- exsitu_ll %>%
  filter(data_source == "ex_situ_BG_survey" &
         inst_short != "NatlPlantGermplasmSystem")
  nrow(bg_exsitu) #833
gene_exsitu <- exsitu_ll %>%
  filter(data_source != "ex_situ_BG_survey" |
         inst_short == "NatlPlantGermplasmSystem")
  nrow(gene_exsitu) #469

# get state boundary data
state_bound <- ne_states(country="united states of america")

### Juglans californica
  # select species pts of interest
  gene_pts <- gene_exsitu %>% filter(taxon_name_acc == "Juglans californica")
    nrow(gene_pts) #15
  bg_pts <- bg_exsitu %>% filter(taxon_name_acc == "Juglans californica")
    nrow(bg_pts) #17
  # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"inputs","PNAS_2020_SDMs",
    "Juglans_californica_PNAS_2020_SDM.tif"))
  # create maps
  map_sp(wild_dist,state_bound,
         gene_pts,"#3e146e",
         bg_pts,"transparent")
  map_sp(wild_dist,state_bound,
         gene_pts,"#3e146e",
         bg_pts,"#d615a9") ##ff7875
### Juglans hindsii
  # select species pts of interest
  gene_pts <- gene_exsitu %>% filter(taxon_name_acc == "Juglans hindsii")
    nrow(gene_pts) #13
  bg_pts <- bg_exsitu %>% filter(taxon_name_acc == "Juglans hindsii") %>%
    filter(UID != "TrompenburgGArb~22913~W~Juglans hindsii" &
           UID != "TrompenburgGArb~22896~W~Juglans hindsii" &
           UID != "BedgeburyNatlPinetum~2015/130~W~Juglans hindsii")
    nrow(bg_pts) #5
  # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"inputs","PNAS_2020_SDMs",
    "Juglans_hindsii_PNAS_2020_SDM.tif"))
  # create maps
  map_sp(wild_dist,state_bound,
         gene_pts,"#3e146e",
         bg_pts,"transparent")
  map_sp(wild_dist,state_bound,
         gene_pts,"#3e146e",
         bg_pts,"#d615a9")



##### SPRING 2020 WORKING GROUP MEETING

##
##### just fruit accessions
##

fruit <- data_sel %>% filter(fruit_nut == "Fruit"); nrow(fruit) #5424
sort(unique(fruit$taxon_name_acc))

# number of institutions who responded to survey
length(sort(unique(fruit$inst_short[which(!grepl("[0-9]",fruit$inst_short))]))) #175
  # number in the US/Canada
NA_inst <- fruit %>% distinct(inst_short,inst_country) %>%
  filter(!grepl("[0-9]",inst_short)) %>% filter(grepl("US|CA",inst_country))
  nrow(NA_inst) #96
# total number (including genebanks)
length(sort(unique(fruit$inst_short))) #205

# map institutions
inst_fruit <- fruit %>%
  select(inst_short,inst_lat,inst_long) %>% distinct()

  map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = inst_fruit,
      lng = ~inst_long, lat = ~inst_lat,
      radius = 9,color = "purple",fillOpacity = 0.6,stroke = F) %>%
    addLegend(
      labels = "Institutions that provided fruit data",
      colors = "purple",position = "bottomright",opacity = 1)
  map

## accession-level stats
  # by taxon and provenance type
fruit_acc <- fruit %>%
  mutate(prov_type = recode(prov_type,
    "N" = "Wild",
    "H?" = "Hort./Unknown",
    "W" = "Wild",
    "Z" = "Wild",
    "NG" = "Hort./Unknown",
    "H" = "Hort./Unknown",
    "U" = "Hort./Unknown")) %>%
  group_by(taxon_name_acc,prov_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
  # add taxon list, to see which have no accessions
fruit_names <- data.frame(taxon_name_acc = sort(unique(taxon_list[taxon_list$fruit_nut=="Fruit",]$taxon_name_acc)))
fruit_names
fruit_acc <- Reduce(full_join,list(fruit_acc,fruit_names))
fruit_acc$num_records[is.na(fruit_acc$num_records)] <- 0
fruit_acc$prov_type[is.na(fruit_acc$prov_type)] <- "Wild"
fruit_acc$taxon_name_acc <- gsub("\\+","x ",fruit_acc$taxon_name_acc)
  # plot
fruit_acc %>%
  mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=prov_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR fruit taxa",
         fill = "Provenance Type")
  #ggsave(file.path(main_dir,"outputs","fruit_prov_barchart.png")#,
  #width = NA, height = NA)
  # look at top species
top_sp <- fruit_acc %>%
  group_by(taxon_name_acc) %>%
  summarize(num_records = sum(num_records)) %>%
  arrange(desc(num_records))
top_sp <- as.data.frame(top_sp)[1:5,1]
for(i in 1:length(top_sp)){
  print(top_sp[i])
  top <- fruit %>%
    filter(taxon_name_acc == top_sp[i])
  print(nrow(top))
  top <- top %>%
    group_by(inst_short) %>% count() %>% ungroup() %>%
    rename(num_records = n) %>% arrange(desc(num_records))
  print(head(as.data.frame(top)))
}
## create table to compare to PNAS results
pnas <- fruit %>%
  group_by(taxon_name_acc) %>%
  count(coord_na = is.na(lat_dd)) %>%
  ungroup()
pnas_wide <- pnas %>%
  spread("coord_na","n")
as.data.frame(pnas_wide)
## by taxon and institution type
fruit_inst <- fruit %>%
  mutate(inst_type = recode(inst_type,
    "Gene/Seed Bank - Genesys/WEIWS" = "Gene/Seed Bank",
    "Botanic Garden; Gene/Seed Bank" = "Unknown")) %>%
  group_by(taxon_name_acc,inst_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
fruit_names
fruit_inst <- Reduce(full_join,list(fruit_inst,fruit_names))
fruit_inst$num_records[is.na(fruit_inst$num_records)] <- 0
fruit_inst$inst_type[is.na(fruit_inst$inst_type)] <- "Unknown"
fruit_inst$taxon_name_acc <- gsub("\\+","x ",fruit_inst$taxon_name_acc)
fruit_inst <- fruit_inst %>% separate("taxon_name_acc","genus",sep=" ",remove=F)
  # plot
fruit_inst %>%
  #mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=inst_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    #facet_wrap(vars(genus), scales="free_y") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR fruit taxa",
         fill = "Institution Type")
## map all points
fruit_pts <- fruit %>%
  distinct(taxon_name_acc,lat_dd,long_dd,.keep_all=T) %>%
  filter(UID != "MEX208~1031700021~NG~Prunus serotina")
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = fruit_pts,
      lng = ~long_dd, lat = ~lat_dd,
      radius = 6,color = "green",fillOpacity = 0.6,stroke = F,
      popup = fruit_pts$UID) %>%
    addControl("North American fruit tree CWR: All wild collection locations",
      position = "topright")
  map
## map specific species as examples
  # get state boundary data
state_bound <- ne_states(country="united states of america")
### Asimina triloba
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Asimina triloba" & !is.na(lat_dd))
  nrow(sp_pts) #29
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Asimina triloba__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Diospyros texana
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Diospyros texana" & !is.na(lat_dd))
  nrow(sp_pts) #12
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Diospyros texana__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Persea borbonia
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Persea borbonia" & !is.na(lat_dd))
  nrow(sp_pts) #4
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Persea borbonia__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Prunus virginiana
    # select species pts of interest
  sp_pts <- fruit_pts %>%
    filter(taxon_name_acc == "Prunus virginiana" & !is.na(lat_dd))
  nrow(sp_pts) #34
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Prunus virginiana__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)


##
##### just nut accessions
##

nut <- data_sel %>% filter(fruit_nut == "Nut"); nrow(nut) #3829

# number of institutions who responded to survey
length(sort(unique(nut$inst_short[which(!grepl("[0-9]",nut$inst_short))]))) #171
  # number in the US/Canada
NA_inst <- nut %>% distinct(inst_short,inst_country) %>%
  filter(!grepl("[0-9]",inst_short)) %>% filter(grepl("US|CA",inst_country))
  nrow(NA_inst) #93
# total number (including genebanks)
length(sort(unique(nut$inst_short))) #187

# map institutions
inst_nut <- nut %>%
  select(inst_short,inst_lat,inst_long) %>% distinct()

  map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = inst_nut,
      lng = ~inst_long, lat = ~inst_lat,
      radius = 9,color = "purple",fillOpacity = 0.6,stroke = F) %>%
    addLegend(
      labels = "Institutions that provided nut data",
      colors = "purple",position = "bottomright",opacity = 1)
map

  # by taxon and provenance type
nut_acc <- nut %>%
  mutate(prov_type = recode(prov_type,
    "N" = "Wild",
    "H?" = "Hort./Unknown",
    "W" = "Wild",
    "Z" = "Wild",
    "NG" = "Hort./Unknown",
    "H" = "Hort./Unknown",
    "U" = "Hort./Unknown")) %>%
  group_by(taxon_name_acc,prov_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
  # add taxon list, to see which have no accessions
nut_names <- data.frame(taxon_name_acc = sort(unique(taxon_list[taxon_list$fruit_nut=="Nut",]$taxon_name_acc)))
nut_names
nut_acc <- Reduce(full_join,list(nut_acc,nut_names))
nut_acc$num_records[is.na(nut_acc$num_records)] <- 0
nut_acc$prov_type[is.na(nut_acc$prov_type)] <- "Wild"
nut_acc$taxon_name_acc <- gsub("\\+","x ",nut_acc$taxon_name_acc)
nut_inst <- nut_inst %>% separate("taxon_name_acc","genus",sep=" ",remove=F)
  # plot
nut_acc %>%
  mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=prov_type, y=num_records, x=taxon_name_acc)) +
    #geom_bar(position="stack", stat="identity") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR nut taxa",
         fill = "Provenance Type")
  # look at top species
top_sp <- nut_acc %>%
  group_by(taxon_name_acc) %>%
  summarize(num_records = sum(num_records)) %>%
  arrange(desc(num_records))
top_sp <- as.data.frame(top_sp)[1:5,1]
for(i in 1:length(top_sp)){
  print(top_sp[i])
  top <- nut %>%
    filter(taxon_name_acc == top_sp[i])
  print(nrow(top))
  top <- top %>%
    group_by(inst_short) %>% count() %>% ungroup() %>%
    rename(num_records = n) %>% arrange(desc(num_records))
  print(head(as.data.frame(top)))
}
  # create table to compare to PNAS results
pnas <- nut %>%
  group_by(taxon_name_acc) %>%
  count(coord_na = is.na(lat_dd)) %>%
  ungroup()
pnas_wide <- pnas %>%
  spread("coord_na","n")
as.data.frame(pnas_wide)
  # by taxon and instition type
nut_inst <- nut %>%
  mutate(inst_type = recode(inst_type,
    "Gene/Seed Bank - Genesys/WEIWS" = "Gene/Seed Bank",
    "Botanic Garden; Gene/Seed Bank" = "Unknown")) %>%
  group_by(taxon_name_acc,inst_type) %>%
  count() %>%
  ungroup() %>%
  rename(num_records = n)
nut_names
nut_inst <- Reduce(full_join,list(nut_inst,nut_names))
nut_inst$num_records[is.na(nut_inst$num_records)] <- 0
nut_inst$inst_type[is.na(nut_inst$inst_type)] <- "Unknown"
nut_inst$taxon_name_acc <- gsub("\\+","x ",nut_inst$taxon_name_acc)
nut_inst <- nut_inst %>% separate()
  # plot
nut_inst %>%
  #mutate(taxon_name_acc = fct_reorder(taxon_name_acc, desc(num_records))) %>%
  ggplot(aes(fill=inst_type, y=num_records, x=taxon_name_acc)) +
    geom_bar(position="stack", stat="identity") +
    #facet_grid(vars(genus), scales="free_y") +
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.text = element_text(size = rel(1.6)),
      legend.title = element_text(size = rel(1.8))) +
    #scale_color_manual(values = c("yellow", "green")) +
    labs(y = "Number of ex situ records",
         x = "Target CWR nut taxa",
         fill = "Institution Type")
  # map all points
nut_pts <- nut %>%
  distinct(taxon_name_acc,lat_dd,long_dd,.keep_all=T) #%>%
  #filter(UID != "MEX208~1031700021~NG~Prunus serotina")
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    addCircleMarkers(
      data = nut_pts,
      lng = ~long_dd, lat = ~lat_dd,
      radius = 6,color = "green",fillOpacity = 0.6,stroke = F,
      popup = nut_pts$UID)
  map
## map specific species as examples
  # get state boundary data
#state_bound <- ne_states(country="united states of america")
### Juglans microcarpa
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Juglans microcarpa" & !is.na(lat_dd))
  nrow(sp_pts) #10
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Juglans microcarpa__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Corylus cornuta
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Corylus cornuta" & !is.na(lat_dd))
  nrow(sp_pts) #18
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Corylus cornuta__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Castanea pumila
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Castanea pumila" & !is.na(lat_dd))
  nrow(sp_pts) #9
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Castanea pumila__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
### Carya glabra
    # select species pts of interest
  sp_pts <- nut_pts %>%
    filter(taxon_name_acc == "Carya glabra" & !is.na(lat_dd))
  nrow(sp_pts) #26
    # read in PNAS raster data
  wild_dist <- raster(file.path(main_dir,"Carya glabra__thrsld_median.tif"))
    # create map
  map_sp(wild_dist,state_bound,sp_pts)
