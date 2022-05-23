### Author: Emily Beckman Bruns  ###  Date: 01/26/2022
### Funding:

### DESCRIPTION:
  #

### INPUTS:
  #

### OUTPUTS:
  #


#################
### LIBRARIES ###
#################

my.packages <- c('plyr','xlsx','tidyverse','stringr','textclean','data.table')
#'readxl','leaflet','raster','sp','rgeos','rgdal','RColorBrewer','ggplot2','maptools')
lapply(my.packages, require, character.only=TRUE)
#  install.packages(my.packages) #Turn on to install current versions

select <- dplyr::select

#################
### FUNCTIONS ###
#################


##############
### SCRIPT ###
##############

main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/BGCI-US/CWR North America (2021-2022)/Ex situ survey"

###############################################################################
# Load data
###############################################################################

# read in target species list
target_sp <- read.csv(file.path(main_dir,
  "North America CWR taxa list with synonyms.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
head(target_sp); nrow(target_sp) #266
length(unique(target_sp$Taxon.accepted)) #90

# read in PlantSearch data for High Priority taxa
ps <- read.csv(file.path(main_dir, "CWR PS match 01.13.22.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
head(ps); nrow(ps) #3576

# read in all GardenSearch data
gs <- read.csv(file.path(main_dir, "All-GSgardens_Jan25-2022.csv"),
  header = T, na.strings=c("","NA"), colClasses="character")
head(gs); nrow(gs) #3736

###############################################################################
# Create summary table with taxa at each instittuion (for mail merge)
###############################################################################

# replace strange characters in taxon names
ps$Taxon.including.synonyms <- gsub(".var\\.."," var. ",ps$Taxon.including.synonyms)
ps$Taxon.including.synonyms <- gsub(".subsp\\.."," subsp. ",ps$Taxon.including.synonyms)
ps$Taxon.including.synonyms <- gsub(".f\\.."," f. ",ps$Taxon.including.synonyms)

# summary table with species at each institution
ps_pivot <- ps %>%
  arrange(Taxon.including.synonyms) %>%
  distinct(Institution,Taxon.including.synonyms,.keep_all=T) %>%
  group_by(Institution) %>%
  mutate(target_taxa = paste0(Taxon.including.synonyms,collapse = ", ")) %>%
  mutate(accepted_taxa = paste0(Taxon.accepted,collapse = ", ")) %>%
  add_tally() %>%
  ungroup() %>%
  distinct(Institution,.keep_all=T) %>%
  select(-Taxon.accepted,-Taxon.including.synonyms,
    -PlantID,-GenHyb,-Gen,-SpecHyb,-Spec,-infRank,-InfraEPI,-Cultivar) %>%
  arrange(Institution)
head(as.data.frame(ps_pivot)); nrow(ps_pivot) #422

# fix institution names
ps_pivot <- ps_pivot %>%
  separate("Institution","Institution",sep=" / \\?|, The",remove=T)
unique(ps_pivot$Institution)

# remove titles after contact names
ps_pivot <- ps_pivot %>%
  separate("DirectorName","DirectorName",sep=", |~|\\(",remove=T) %>%
  separate("CuratorName","CuratorName",sep=", |~|\\(",remove=T) %>%
  separate("PlantRecordsOfficerName","PlantRecordsOfficerName",sep=", |~|\\(",remove=T)

# clean up contact name and email columns
  # if curator and plant records emails are the same, remove plant records email
ps_pivot[which(ps_pivot$PlantRecordsOfficerEmail ==
  ps_pivot$CuratorEmail),]$PlantRecordsOfficerEmail <- NA
  # if primary and director emails are the same, remove primary email
ps_pivot[which(ps_pivot$DirectorEmail ==
  ps_pivot$PrimaryEmail),]$PrimaryEmail <- NA
  # if director name but no email, remove name
ps_pivot[which(!is.na(ps_pivot$DirectorName) &
  is.na(ps_pivot$DirectorEmail)),]$DirectorName <- NA
  # if curator email but no name, add "Curator" as name
ps_pivot[which(!is.na(ps_pivot$CuratorEmail) &
  is.na(ps_pivot$CuratorName)),]$CuratorName <- "Curator"
  # if curator name but no email, remove name
ps_pivot[which(!is.na(ps_pivot$CuratorName) &
  is.na(ps_pivot$CuratorEmail)),]$CuratorName <- NA
  # if plant records email but no name, add "Plant Records Officer" as name
ps_pivot[which(!is.na(ps_pivot$PlantRecordsOfficerEmail) &
  is.na(ps_pivot$PlantRecordsOfficerName)),]$PlantRecordsOfficerName <- "Plant Records Officer"
  # if plant records name but no email, remove name
ps_pivot[which(!is.na(ps_pivot$PlantRecordsOfficerName) &
  is.na(ps_pivot$PlantRecordsOfficerEmail)),]$PlantRecordsOfficerName <- NA

# fill contact email and name columns

  ## Curator and Plant Records Officer names
ps_pivot$PlantRecordsOfficerName[which(is.na(ps_pivot$PlantRecordsOfficerName))] <- ""
ps_pivot$CuratorName[which(is.na(ps_pivot$CuratorName))] <- ""
ps_pivot$contact_name <- NA
for(i in 1:nrow(ps_pivot)){
  if(ps_pivot$PlantRecordsOfficerName[i] == ""){
    ps_pivot$contact_name[i] <- str_trim(str_to_title(ps_pivot$CuratorName[i]))
  } else if(ps_pivot$CuratorName[i] == ""){
    ps_pivot$contact_name[i] <- str_trim(str_to_title(ps_pivot$PlantRecordsOfficerName[i]))
  } else if (str_trim(str_to_title(ps_pivot$PlantRecordsOfficerName[i])) !=
      str_trim(str_to_title(ps_pivot$CuratorName[i]))){
    ps_pivot$contact_name[i] <- paste0(
      str_trim(str_to_title(ps_pivot$CuratorName[i]))," (",ps_pivot$CuratorEmail[i],") and ",
      str_trim(str_to_title(ps_pivot$PlantRecordsOfficerName[i]))," (",ps_pivot$PlantRecordsOfficerEmail[i],")")
  } else {
    ps_pivot$contact_name[i] <- str_trim(str_to_title(ps_pivot$CuratorName[i]))
  }
}
ps_pivot$contact_name[which(ps_pivot$contact_name == "")] <- NA
  # check out the results
sort(unique(ps_pivot$contact_name))
  # replace anything strange
ps_pivot$contact_name <- gsub("Secretary   : ","",ps_pivot$contact_name,)
str(ps_pivot)
  # create email column
plantrecords <- ps_pivot
plantrecords$contact_email <- plantrecords$PlantRecordsOfficerEmail
plantrecords <- plantrecords[which(!is.na(plantrecords$contact_email)),]
plantrecords$contact_type <- "Plant Records Officer"
nrow(plantrecords) #186
curator <- ps_pivot
curator$contact_email <- curator$CuratorEmail
curator <- curator[which(!is.na(curator$contact_email)),]
curator$contact_type <- "Curator"
nrow(curator) #325
ps_join <- full_join(plantrecords,curator)
nrow(ps_join) #511

  # only add director if no other contact, primary as last resort
add <- ps_pivot[which(is.na(ps_pivot$contact_name)),]
  ## Director names
add$contact_name <- str_trim(str_to_title(add$DirectorName))
add$contact_name[which(is.na(add$contact_name))] <- "Colleague"
add$contact_name <- gsub(" And "," and ",add$contact_name)
add$contact_name <- gsub("^Xx$","Colleague",add$contact_name)
add$contact_name <- gsub("^-$","Colleague",add$contact_name)
add$contact_type <- "Director"
unique(add$contact_name)
  # create email column
add$contact_email <- add$DirectorEmail
add[which(is.na(add$DirectorEmail)),]$contact_email <-
  add[which(is.na(add$DirectorEmail)),]$PrimaryEmail
add[which(is.na(add$DirectorEmail) & !is.na(add$contact_email)),]$contact_type <- "Primary"
add[which(is.na(add$contact_email)),]
  # join everything together
ps_join2 <- full_join(ps_join,add)
ps_join2[which(is.na(ps_join2$contact_email)),]$contact_type <- "X NO CONTACT INFO"
  # arrange by number of species
ps_join2 <- ps_join2 %>% arrange(desc(n),contact_type)
nrow(ps_join2) #566
str(ps_join2)

# final edits

### !! GARDENS WITH NO CONTACT INFO & NUMBER OF TARGET TAXA AT EACH !!
ps_join2[which(is.na(ps_join2$contact_email)),c(1,15)]
  # add contacts found online
ps_join2[which(ps_join2$Institution == "Forstbotanischer Garten Eberswalde"),16:18] <-
 list("Mr. Prof. Dr. Matthias Barth","buero.praesident@hnee.de","Director")
ps_join2[which(ps_join2$Institution == "Forstbotanischer Garten und Arboretum"),16:18] <-
 list("Prof. Dr. Andrea Polle","apolle@gwdg.de","Curator")
ps_join2[which(ps_join2$Institution == "Botanischer Garten der Universität Freiburg"),16:18] <-
 list("Colleague","botanischer.garten@biologie.uni-freiburg.de","Primary")

### !! NETWORKS !!
ps_join2[which(ps_join2$InstitutionType == "Network"),1:2]
  # remove netwoks that overlap with other institutions we are contacting
ps_final <- ps_join2[which(ps_join2$Institution != "Center for Plant Conservation (USA)" &
                           ps_join2$Institution != "Bishop Museum - Checklist of Cultivated Plants of Hawai'i"),]
nrow(ps_final) #563

# set up dummy row as first row, with longest inputs in each column,
#   to be sure nothing gets trimmed during mail merge
instl <- ps_final[nchar(ps_final$Institution)==max(nchar(ps_final$Institution)),]$Institution[1]
namel <- ps_final[nchar(ps_final$contact_name)==max(nchar(ps_final$contact_name)),]$contact_name[1]
taxal <- ps_final[nchar(ps_final$target_taxa)==max(nchar(ps_final$target_taxa)),]$target_taxa[1]
dummy_row <- data.frame(Institution=instl, contact_name=namel, target_taxa=taxal,
  contact_email="ebeckman@mortonarb.org")
ps_final <- full_join(dummy_row,ps_final)

# write file
write.xlsx(ps_final, file.path(main_dir,"Mail Merge",
  paste0("NA-CWR_ExSituSurvey_DirectEmail_Contacts_",Sys.Date(), ".xlsx")),
  col.names = T, row.names=F, append = F)

###############################################################################
# Create list of U.S. and Canada institutions not already contacted (for secondary mail merge)
###############################################################################

# remove insitutions already contacted
gs_nam <- gs %>% filter(CountryCode == "US" | CountryCode == "CA")
nrow(gs_nam) #1156
direct_contact <- unique(ps$Institution)
gs_contact <- gs_nam %>% filter(!(Institution %in% direct_contact))
nrow(gs_contact) #942
unique(gs_contact$InstitutionType)

# look at networks
gs_contact[which(gs_contact$InstitutionType == "Network"),]
# remove networks
gs_contact <- gs_contact %>% filter(InstitutionType != "Network")
nrow(gs_contact) #935

# fix institution names
gs_contact <- gs_contact %>%
  separate("Institution","Institution",sep=" / \\?|, The",remove=T)

# remove titles after contact names
gs_contact <- gs_contact %>%
  separate("DirectorName","DirectorName",sep=", |~|\\(",remove=T) %>%
  separate("CuratorName","CuratorName",sep=", |~|\\(",remove=T) %>%
  separate("PlantRecordsOfficerName","PlantRecordsOfficerName",sep=", |~|\\(",remove=T)

# clean up contact name and email columns
  # if curator and plant records emails are the same, remove plant records email
gs_contact[which(gs_contact$PlantRecordsOfficerEmail ==
  gs_contact$CuratorEmail),]$PlantRecordsOfficerEmail <- NA
  # if primary and director emails are the same, remove primary email
gs_contact[which(gs_contact$DirectorEmail ==
  gs_contact$PrimaryEmail),]$PrimaryEmail <- NA
  # if director name but no email, remove name
gs_contact[which(!is.na(gs_contact$DirectorName) &
  is.na(gs_contact$DirectorEmail)),]$DirectorName <- NA
  # if curator email but no name, add "Curator" as name
gs_contact[which(!is.na(gs_contact$CuratorEmail) &
  is.na(gs_contact$CuratorName)),]$CuratorName <- "Curator"
  # if curator name but no email, remove name
gs_contact[which(!is.na(gs_contact$CuratorName) &
  is.na(gs_contact$CuratorEmail)),]$CuratorName <- NA
  # if plant records email but no name, add "Plant Records Officer" as name
gs_contact[which(!is.na(gs_contact$PlantRecordsOfficerEmail) &
  is.na(gs_contact$PlantRecordsOfficerName)),]$PlantRecordsOfficerName <- "Plant Records Officer"
  # if plant records name but no email, remove name
gs_contact[which(!is.na(gs_contact$PlantRecordsOfficerName) &
  is.na(gs_contact$PlantRecordsOfficerEmail)),]$PlantRecordsOfficerName <- NA

# fill contact email and name columns

  ## Curator and Plant Records Officer names
gs_contact$PlantRecordsOfficerName[which(is.na(gs_contact$PlantRecordsOfficerName))] <- ""
gs_contact$CuratorName[which(is.na(gs_contact$CuratorName))] <- ""
gs_contact$contact_name <- NA
for(i in 1:nrow(gs_contact)){
  if(gs_contact$PlantRecordsOfficerName[i] == ""){
    gs_contact$contact_name[i] <- str_trim(str_to_title(gs_contact$CuratorName[i]))
  } else if(gs_contact$CuratorName[i] == ""){
    gs_contact$contact_name[i] <- str_trim(str_to_title(gs_contact$PlantRecordsOfficerName[i]))
  } else if (str_trim(str_to_title(gs_contact$PlantRecordsOfficerName[i])) !=
      str_trim(str_to_title(gs_contact$CuratorName[i]))){
    gs_contact$contact_name[i] <- paste0(
      str_trim(str_to_title(gs_contact$CuratorName[i]))," (",gs_contact$CuratorEmail[i],") and ",
      str_trim(str_to_title(gs_contact$PlantRecordsOfficerName[i]))," (",gs_contact$PlantRecordsOfficerEmail[i],")")
  } else {
    gs_contact$contact_name[i] <- str_trim(str_to_title(gs_contact$CuratorName[i]))
  }
}
gs_contact$contact_name[which(gs_contact$contact_name == "")] <- NA
  # check out the results
sort(unique(gs_contact$contact_name))
  # create email column
plantrecords <- gs_contact
plantrecords$contact_email <- plantrecords$PlantRecordsOfficerEmail
plantrecords <- plantrecords[which(!is.na(plantrecords$contact_email)),]
plantrecords$contact_type <- "Plant Records Officer"
nrow(plantrecords) #96
curator <- gs_contact
curator$contact_email <- curator$CuratorEmail
curator <- curator[which(!is.na(curator$contact_email)),]
curator$contact_type <- "Curator"
nrow(curator) #398
gs_join <- full_join(plantrecords,curator)
nrow(gs_join) #494

  # only add director if no other contact, primary as last resort
add <- gs_contact[which(is.na(gs_contact$contact_name)),]
  ## Director names
add$contact_name <- str_trim(str_to_title(add$DirectorName))
add$contact_name[which(is.na(add$contact_name))] <- "Colleague"
add$contact_name <- gsub(" And "," and ",add$contact_name)
add$contact_name <- gsub("^Xx$","Colleague",add$contact_name)
add$contact_name <- gsub("^-$","Colleague",add$contact_name)
add$contact_type <- "Director"
unique(add$contact_name)
  # create email column
add$contact_email <- add$DirectorEmail
add[which(is.na(add$DirectorEmail)),]$contact_email <-
  add[which(is.na(add$DirectorEmail)),]$PrimaryEmail
add[which(is.na(add$DirectorEmail) & !is.na(add$contact_email)),]$contact_type <- "Primary"
add[which(is.na(add$contact_email)),]
  # join everything together
gs_join2 <- full_join(gs_join,add)
gs_join2[which(is.na(gs_join2$contact_email)),]$contact_type <- "X NO CONTACT INFO"
nrow(gs_join2) #1011
str(gs_join2)

# final edits

### !! GARDENS WITH NO CONTACT INFO !!
gs_join2[which(is.na(gs_join2$contact_email)),2]

# set up dummy row as first row, with longest inputs in each column,
#   to be sure nothing gets trimmed during mail merge
instl <- gs_join2[nchar(gs_join2$Institution)==max(nchar(gs_join2$Institution)),]$Institution[1]
namel <- gs_join2[nchar(gs_join2$contact_name)==max(nchar(gs_join2$contact_name)),]$contact_name[1]
dummy_row <- data.frame(Institution=instl, contact_name=namel,
  contact_email="ebeckman@mortonarb.org")
gs_final <- full_join(dummy_row,gs_join2)

# add entry that sends to Abby, so she can forward along
abby_row <- data.frame(Institution="your institution", contact_name="Colleague",
  contact_email="abby.meyer@bgci.org")
gs_final <- full_join(gs_final,abby_row)

# write file
write.xlsx(gs_final, file.path(main_dir,"Mail Merge",
  paste0("NA-CWR_ExSituSurvey_BroaderEmail_Contacts_",Sys.Date(), ".xlsx")),
  col.names = T, row.names=F, append = F)













################################################################################
# Map state richness
################################################################################

#setwd("./../../../../..")
local_dir <- "Desktop/work"
drive_dir <- "/Volumes/GoogleDrive/My Drive/BGCI-US SEPCA Project 2021/Underlying data for analyses/State species richness and garden locations interactive maps"

# READ IN POLYGONS

# read in shapefiles
	# https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
states_shp <- readOGR(file.path(local_dir,"cb_2018_us_state_20m/cb_2018_us_state_20m.shp"))
	# https://www.arcgis.com/home/item.html?id=2ca75003ef9d477fb22db19832c9554f
world_shp <- readOGR(file.path(local_dir,"countries_shp/countries.shp"))

# READ IN SPECIES RICHNESS DATA

# read in distribution information
state_richness <- read.csv(file.path(local_dir,"SE_PCA_richness - SE_PCA_richness.csv"),
  as.is=T, na.strings=c("","NA"))
state_richness <- read.csv(file.path(local_dir,"SE_PCA_richness_No.PS.or.survey.inst.csv"),
  as.is=T, na.strings=c("","NA"))

# JOIN STATE/COUNTRY POLYGONS AND RICHNESS DATA

states_joined <- merge(states_shp,state_richness)
# select states with taxa
#filled_states <- state_richness[which(state_richness$num_taxa != "0"),]$STUSPS
#states_filled <- states_joined[states_joined@data$STUSPS %in% filled_states,]
#head(states_filled)
# select states with no taxa
#empty_states <- target_states <- state_richness[which(state_richness$num_taxa == "0"),]$STUSPS
#states_empty <- states_joined[states_joined@data$STUSPS %in% empty_states,]
#head(states_empty)

# select target states
target_states <- state_richness[which(state_richness$target_state == "Y"),]$STUSPS
states_target <- states_joined[states_joined@data$STUSPS %in% target_states,]
head(states_target)
# select non-target states
empty_states <- state_richness[which(state_richness$target_state == "N"),]$STUSPS
states_empty <- states_joined[states_joined@data$STUSPS %in% empty_states,]
head(states_empty)

# MAP

# function to create color palettes
create_pal <- function(spdf_col,palette,rev){
	pal <- colorBin(palette = palette, bins = bins,
		domain = spdf_col, reverse = rev, na.color = "transparent")
	return(pal)
}
display.brewer.all()

# create color bins and labels, then color palette
bins <- c(1,20,50,100,200,Inf)
labels <- c("1-19","20-49","50-99","100-199","200+")
pal <- create_pal(states_target@data$num_taxa,"YlOrBr",F)

# function to create maps
map.richness <- function(boundary_poly,empty_boundary,#target_boundary,
    boundary_richness,my_palette,labels,legend_txt){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
		#addPolygons(data = empty_boundary,
		#	label = ~NAME,
		#	color = "#404040",
		#	weight = 1,
		#	opacity = 1,
		#	fillColor = "white",
		#	fillOpacity = 0) %>%
	  addPolygons(data = boundary_poly,
			label = ~NAME,
			color = "#404040",
			weight = 5,
			opacity = 1,
			fillColor = ~my_palette(boundary_richness),
			fillOpacity = 0.8) %>%
		#addPolygons(data = target_boundary,
		#	label = ~NAME,
		#	color = "black",
		#	weight = 5,
		#	opacity = 1,
		#	fillColor = "white",
		#	fillOpacity = 0) %>%
		leafem::addStaticLabels(
			.,
			data = boundary_poly,
			label = boundary_poly@data$num_taxa,
			style = list("font-weight"="bold","font-size"="20px","color"="black")) %>%
		addLegend(pal = my_palette,
			values = boundary_richness,
			opacity = 0.8,
			title = legend_txt,
			labFormat = function(type, cuts, p) {paste0(labels)},
			position = "bottomright") %>%
    setView(-98, 40, zoom = 5)
	return(map)
}

richness_map <- map.richness(states_target,states_empty,
  states_target@data$num_taxa,pal,labels,"Number of high priority taxa")
richness_map
htmlwidgets::saveWidget(richness_map, file = file.path(main_dir,"State_Taxon_Richness_SE_PCA.html"))

################################################################################
# Map institutions
################################################################################

# read in summary PlantSearch data

# Y/N high priority taxa
#inst <- read.csv(file.path(drive_dir,"All gardens by Y:N high priority taxa.csv"),
#  as.is=T, na.strings=c("","NA"))
#pre_2010 <- inst %>% filter(UpdateYearCategory == "pre-2010-")
#priority_taxa <- inst %>% filter(UpdateYearCategory != "pre-2010-" & high_priority_taxa. == "Y")
#no_priority_taxa <- inst %>% filter(UpdateYearCategory != "pre-2010-" & high_priority_taxa. == "N")

# Y/N wild material
inst <- read.csv(file.path(local_dir,"All gardens wild_yn.csv"),
  as.is=T, na.strings=c("","NA"))
inst$Latitude <- as.numeric(inst$Latitude)
inst$Longitude <- as.numeric(inst$Longitude)
inst <- inst %>% filter(!is.na(Latitude) & !is.na(Longitude))

pre_2010 <- inst %>% filter(UpdateYearCategory == "pre-2010-")
wild_material <- inst %>% filter(UpdateYearCategory != "pre-2010-" & wild_yn == "Y")
no_wild_material <- inst %>% filter(UpdateYearCategory != "pre-2010-" & wild_yn == "N")
no_target_taxa <- inst %>% filter(UpdateYearCategory != "pre-2010-" & wild_yn == "A-NONE")
# combine bottom two categories
  no_target_taxa <- rbind(no_target_taxa,pre_2010)

map.inst <- function(boundary_poly,pre_2010,wild,not_wild,no_target){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    #addCircleMarkers(
    #  data = pre_2010,
    #  lng = ~Longitude, lat = ~Latitude,
    #  radius = 2,color = "#cbd1f5",fillOpacity = 1,stroke = F,
    #  popup = pre_2010$Institution) %>%
    addCircleMarkers(
      data = no_target,
      lng = ~Longitude, lat = ~Latitude,
      radius = 2,color = "#f9c4c5",fillOpacity = 1,stroke = F, #b89ce6
      popup = no_target$Institution) %>%
    addCircleMarkers(
      data = not_wild,
      lng = ~Longitude, lat = ~Latitude,
      radius = 3.5,color = "#f689b5",fillOpacity = 1,stroke = F, #c752c3
      popup = not_wild$Institution) %>%
    addCircleMarkers(
      data = wild,
      lng = ~Longitude, lat = ~Latitude,
      radius = 4,color = "#A80874",fillOpacity = 1,stroke = F, #700348
      popup = wild$Institution) %>%
    addLegend(
      labels = c("wild-origin material for Southeast U.S. high priority taxa",
                 "Southeast U.S. high priority taxa but no wild material",
                 "no known material for Southeast U.S. high priority taxa"#,
                 #"no collections data available"
                ),
      colors = c("#A80874","#f689b5","#f9c4c5"),
      title = "Ex situ living collections with:",
      position = "bottomright",
      opacity = 1)
	return(map)
}

inst_map <- map.inst(states_target,pre_2010,wild_material,no_wild_material,no_target_taxa)
#htmlwidgets::saveWidget(inst_map, file = file.path(main_dir,"PlantSearch_Institutions_SE_PCA.html"))
inst_map

# version with state lines
map.state <- function(boundary_poly,pre_2010,wild,not_wild,no_target){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
    #addCircleMarkers(
    #  data = pre_2010,
    #  lng = ~Longitude, lat = ~Latitude,
    #  radius = 4,color = "#cbd1f5",fillOpacity = 1,stroke = F,
    #  popup = pre_2010$Institution) %>%
    addCircleMarkers(
      data = no_target,
      lng = ~Longitude, lat = ~Latitude,
      radius = 4,color = "#f9c4c5",fillOpacity = 1,stroke = F,
      popup = no_target$Institution) %>%
    addCircleMarkers(
      data = not_wild,
      lng = ~Longitude, lat = ~Latitude,
      radius = 5.5,color = "#f689b5",fillOpacity = 1,stroke = F,
      popup = not_wild$Institution) %>%
    addCircleMarkers(
      data = wild,
      lng = ~Longitude, lat = ~Latitude,
      radius = 6,color = "#A80874",fillOpacity = 1,stroke = F,
      popup = wild$Institution) %>%
		addPolygons(data = boundary_poly,
			color = "#545454",weight = 1,	opacity = 1,
		  fillColor = "#e3e3d8",fillOpacity = 0.2) #%>%
    #addLegend(
    #  labels = c("wild-origin material for Southeast U.S. high priority taxa",
    #             "Southeast U.S. high priority taxa but no wild material",
    #             "no Southeast U.S. high priority taxa",
    #             "no collections data available"),
    #  colors = c("#700348","#a733ab","#a191cf","#cbd1f5"),
    #  title = "Ex situ living collections with:",
    #  position = "bottomright",
    #  opacity = 1)
	return(map)
}

state_map <- map.state(states_target,pre_2010,wild_material,no_wild_material,no_target_taxa)
#htmlwidgets::saveWidget(inst_map, file = file.path(main_dir,"PlantSearch_Institutions_SE_PCA.html"))
state_map

################################################################################
# Case study mapping
################################################################################

#insitu <- read.csv(file.path(drive_dir,"insitu_pts_case_studies.csv"),
#  as.is=T, na.strings=c("","NA"))
exsitu <- read.csv(file.path(drive_dir,"exsitu_compiled_standardized_2021-04-23_Wild_ex-wild_origin_records.csv"),
  as.is=T, na.strings=c("","NA"))
exsitu <- exsitu %>% filter(!is.na(lat_dd) & !is.na(long_dd))

spp.all <- c("Tsuga caroliniana","Platanthera integrilabia","Helianthus verticillatus",
  "Macbridea caroliniana","Ribes echinellum","Schwalbea americana","Arabis georgiana")
spp.ptsize <- c(12,12,12,6,12,8,12)
spp.zoom <- c(7,7,7,6,7,5,7)

#i <- 1

# cycle through each species file and create map
for(i in 1:length(spp.all)){

  # select records
  #spp.now_i <- insitu %>% filter(species == spp.all[i]) %>%
  #  rename(taxon_name_acc = species,localityi = locality,UID = gbifID)
  #  spp.now_i$database <- "GBIF"
  spp.now_e <- exsitu %>% filter(taxon_name_acc == spp.all[i]) %>%
    rename(decimalLatitude = lat_dd, decimalLongitude = long_dd)
  #  spp.now_e$database <- "Ex_situ"
  #spp.now <- rbind.fill(spp.now_i,spp.now_e)

  ## palette based on database
  # set database as factor and order appropriately
  #spp.now$database <- factor(spp.now$database,
  #  levels = c("Ex_situ","GBIF"))
  #spp.now <- spp.now %>% arrange(desc(database))
  # create color palette
  #colors <- c("#cf8d5f","#c1c46c")
  #database.pal <- colorFactor(palette=colors,
  #  levels=c("Ex_situ","GBIF"))

  # create map
  map <- leaflet() %>%
    # Base layer
    addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
    addControl(paste0("<b>",spp.all[i]), position = "topright") %>%
    # Add ex situ points
    addCircleMarkers(data = spp.now_e, ~decimalLongitude, ~decimalLatitude,
      popup = paste0(
        "<b>UID:</b> ",spp.now_e$UID,"<br/>",
        "<b>Institution:</b> ",spp.now_e$inst_short,"<br/>",
        "<b>Collection Year:</b> ",spp.now_e$coll_year,"<br/>",
        "<b>Original source:</b> ",spp.now_e$orig_source,"<br/>",
        "<b>Germplasm type:</b> ",spp.now_e$germ_type,"<br/>",
        "<b>Geolocation method:</b> ",spp.now_e$latlong_flag,"<br/>",
        "<b>Geolocation notes:</b> ",spp.now_e$geolocation_notes),
      color = "#b80656",
      radius = spp.ptsize[i],fillOpacity = 0.6,stroke = F) %>%
    # Add legend
	  addLegend(labels = "Source locality for wild-origin ex situ accession",
		  colors = "#b80656",title = "Legend",position = "bottomright",opacity = 0.6) %>%
    setView(-84, 34, zoom = spp.zoom[i])
  map

  # save map
  htmlwidgets::saveWidget(map, file.path(local_dir,
    paste0(spp.all[i], "_ExsituWildSources_InteractiveMap.html")))

  cat("\tEnding ", spp.all[i], ", ", i, " of ", length(spp.all), ".\n\n", sep="")
}





























### OLD ###
## remove duplicates for each institution with multiple PS IDs
# Altanta BG
a <- ps_all[which(ps_all$Institution == "CPC: Atlanta BG "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "1147" |
                  ps_all$GardenID == "5638" |
                  ps_all$GardenID == "5639" |
                  ps_all$GardenID == "5640"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove1 <- ps_all[which(ps_all$Institution == "CPC: Atlanta BG " & !(ps_all$accepted_SCIENTIFIC_NAME == keep)),]
# Bok Tower
a <- ps_all[which(ps_all$Institution == "CPC: Bok Tower "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "4437" |
                  ps_all$GardenID == "4459" |
                  ps_all$GardenID == "1150"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove2 <- ps_all[which(ps_all$Institution == "CPC: Bok Tower " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# Cincinnati Zoo
a <- ps_all[which(ps_all$Institution == "CPC: Cincinnati Zoo "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "4503" |
                  ps_all$GardenID == "4502"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove3 <- ps_all[which(ps_all$Institution == "CPC: Cincinnati Zoo " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# Desert BG
a <- ps_all[which(ps_all$Institution == "CPC: Desert BG "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "69" |
                  ps_all$GardenID == "4508"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove4 <- ps_all[which(ps_all$Institution == "CPC: Desert BG " & !(ps_all$accepted_SCIENTIFIC_NAME == keep)),]
# Fairchild Tropical BG
a <- ps_all[which(ps_all$Institution == "CPC: Fairchild Tropical BG "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "87" |
                  ps_all$GardenID == "5641" |
                  ps_all$GardenID == "5642"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove5 <- ps_all[which(ps_all$Institution == "CPC: Fairchild Tropical BG " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# LBJ Wildflower Center
a <- ps_all[which(ps_all$Institution == "CPC: LBJ Wildflower Center "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "1163" |
                  ps_all$GardenID == "4757"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove6 <- ps_all[which(ps_all$Institution == "CPC: LBJ Wildflower Center " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# Missouri BG
a <- ps_all[which(ps_all$Institution == "CPC: Missouri BG "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "4515" |
                  ps_all$GardenID == "771"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove7 <- ps_all[which(ps_all$Institution == "CPC: Missouri BG " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# North Carolina BG
a <- ps_all[which(ps_all$Institution == "CPC: North Carolina BG "),]$accepted_SCIENTIFIC_NAME
b <- ps_all[which(ps_all$GardenID == "639" |
                  ps_all$GardenID == "4484"),]$accepted_SCIENTIFIC_NAME
keep <- setdiff(a,b)
remove8 <- ps_all[which(ps_all$Institution == "CPC: North Carolina BG " & !(ps_all$accepted_SCIENTIFIC_NAME %in% keep)),]
# JOIN ALL TOGETHER FOR REMOVAL
remove_a <- list(remove1,remove2,remove3,remove4,remove5,remove6,remove7,remove8)
remove_all <- Reduce(rbind,remove_a)
ps_all2 <- anti_join(ps_all,remove_all)
