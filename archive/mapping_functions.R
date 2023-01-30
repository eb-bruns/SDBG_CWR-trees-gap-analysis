

create.buffers <- function(df,radius,pt_proj,buff_proj,boundary,
                           lat.col,long.col){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c(long.col, lat.col), crs=pt_proj)
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
create.full.map <- function(spp.raster,raster.pal,spp.pts.buffer,
                            spp.pts,ex.pts.buffer,ex.pts,taxon.nm
#                            ,pro_areas,eco_clip
                            ){
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
                     group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
#    addPolygons(data = eco_clip,
#                fillColor = "white", fillOpacity = 0,
#                color = "#dec171", weight = 2, opacity = 0.7,
#                popup = eco_clip$ECO_NAME,
#                group = "Ecoregions (Olson et al. 2001)") %>%
    # Protected areas
#    addRasterImage(pro_areas, colors="#753a10", opacity = 0.7,
#                   group = "Protected Areas (WDPA)") %>%
    # Species Distribution Model
    addRasterImage(spp.raster, colors=raster.pal, opacity = 0.8,
                   group = "Species distribution model") %>%
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~longitude, ~latitude,
                     #popup = ~paste0(
                       #"<b>Accepted taxon name:</b> ",taxon_name_accepted,"<br/>",
                       #"<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       #"<b>Source database:</b> ",database,"<br/>",
                       #"<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
                       #"<b>Year:</b> ",year,"<br/>",
                       #"<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       #"<b>Dataset name:</b> ",datasetName,"<br/>",
                       #"<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       #"<b>Latitude:</b> ",latitude,"<br/>",
                       #"<b>Longitude:</b> ",longitude,"<br/>",
                       #"<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       #"<b>ID:</b> ",UID),
                     color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)") %>%
    # Ex situ buffers
    addPolygons(data = ex.pts.buffer,
                fillColor = "#80bf30", fillOpacity = 0.7,
                weight = 0, color = "white", opacity = 0,
                smoothFactor = 0,
                group = "Estimated area represented by ex situ material (50km buffers around source localities)") %>%
    # Ex situ points
    addCircleMarkers(data = ex.pts, ~longitude, ~latitude,
                     #popup = ~paste0(
                       #"<b>Accepted taxon name:</b> ",taxon_name_accepted,"<br/>",
                       #"<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       #"<b>Source database:</b> ",database,"<br/>",
                       #"<b>Institution name:</b> ",datasetName,"<br/>",
                       #"<b>Accession number:</b> ",references,"<br/>",
                       #"<b>Collection date:</b> ",year,"<br/>",
                       #"<b>Number of individuals:</b> ",individualCount,"<br/>",
                       #"<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       #"<b>ID:</b> ",UID),
                     color = "#335707", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Source localities for ex situ material (in botanic gardens/genebanks)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",taxon.nm)),
               position = "topright") %>%
    addLegend(labels = c(
      "Species distribution model",
      "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
      "Estimated species distribution (50km buffers around occurrence points)",
      "Source localities for ex situ material (in botanic gardens/genebanks)",
      "Estimated area represented by ex situ material (50km buffers around source localities)"
#      ,"Protected Areas (WDPA)",
#      "Ecoregions (Olson et al. 2021)"
       ),
      colors = c("#0e6187","#56116b","#965da7","#335707","#80bf30"
#                ,"#753a10","#dec171"
                 ),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(
        "Species distribution model",
        "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)"
#        ,"Protected Areas (WDPA)",
#        "Ecoregions (Olson et al. 2001)"
         ),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
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
create.noex.map <- function(eco_clip,spp.raster,raster.pal,spp.pts.buffer,
                            spp.pts,pro_areas){
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
                     group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    addPolygons(data = eco_clip,
                fillColor = "white", fillOpacity = 0,
                color = "#dec171", weight = 2, opacity = 0.7,
                popup = eco_clip$ECO_NAME,
                group = "Ecoregions (Olson et al. 2001)") %>%
    # Protected areas
    addRasterImage(pro_areas, colors="#753a10", opacity = 0.7,
                   group = "Protected Areas (WDPA)") %>%
    # Species Distribution Model
    addRasterImage(spp.raster, colors=raster.pal, opacity = 0.8,
                   group = "Species distribution model") %>%
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",taxon_name_accepted,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",decimalLatitude,"<br/>",
                       "<b>Longitude:</b> ",decimalLongitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
               position = "topright") %>%
    addLegend(labels = c(
      "Species distribution model",
      "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
      "Estimated species distribution (50km buffers around occurrence points)",
      "Protected Areas (WDPA)",
      "Ecoregions (Olson et al. 2021)"),
      colors = c("#0e6187","#56116b","#965da7","#753a10","#dec171"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(
        "Species distribution model",
        "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Protected Areas (WDPA)",
        "Ecoregions (Olson et al. 2001)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
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
create.noex.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts,pro_areas){
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
                     group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    addPolygons(data = eco_clip,
                fillColor = "white", fillOpacity = 0,
                color = "#dec171", weight = 2, opacity = 0.7,
                popup = eco_clip$ECO_NAME,
                group = "Ecoregions (Olson et al. 2001)") %>%
    # Protected areas
    addRasterImage(pro_areas, colors="#753a10", opacity = 0.7,
                   group = "Protected Areas (WDPA)") %>%
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",taxon_name_accepted,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",decimalLatitude,"<br/>",
                       "<b>Longitude:</b> ",decimalLongitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
               position = "topright") %>%
    addLegend(labels = c(
      "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
      "Estimated species distribution (50km buffers around occurrence points)",
      "Protected Areas (WDPA)",
      "Ecoregions (Olson et al. 2021)"),
      colors = c("#56116b","#965da7","#753a10","#dec171"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(
        "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Protected Areas (WDPA)",
        "Ecoregions (Olson et al. 2001)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
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
create.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts,ex.pts.buffer,
                             ex.pts,pro_areas){
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$CartoDB.Positron,
                     group = "CartoDB.Positron") %>%
    ### Overlay groups (can toggle)
    # Ecoregions
    addPolygons(data = eco_clip,
                fillColor = "white", fillOpacity = 0,
                color = "#dec171", weight = 2, opacity = 0.7,
                popup = eco_clip$ECO_NAME,
                group = "Ecoregions (Olson et al. 2001)") %>%
    # Protected areas
    addRasterImage(pro_areas, colors="#753a10", opacity = 0.7,
                   group = "Protected Areas (WDPA)") %>%
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated species distribution (50km buffers around occurrence points)") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~decimalLongitude, ~decimalLatitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",taxon_name_accepted,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with duplicate record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",decimalLatitude,"<br/>",
                       "<b>Longitude:</b> ",decimalLongitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)") %>%
    # Ex situ buffers
    addPolygons(data = ex.pts.buffer,
                fillColor = "#80bf30", fillOpacity = 0.7,
                weight = 0, color = "white", opacity = 0,
                smoothFactor = 0,
                group = "Estimated area represented by ex situ material (50km buffers around source localities)") %>%
    # Ex situ points
    addCircleMarkers(data = ex.pts, ~decimalLongitude, ~decimalLatitude,
                     popup = ~paste0(
                       "<b>Accepted species name:</b> ",taxon_name_accepted,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>Institution name:</b> ",datasetName,"<br/>",
                       "<b>Accession number:</b> ",references,"<br/>",
                       "<b>Collection date:</b> ",year,"<br/>",
                       "<b>Number of individuals:</b> ",individualCount,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#335707", radius = 3, fillOpacity = 0.9, #stroke = T,
                     group = "Source localities for ex situ material (in botanic gardens/genebanks)") %>%
    ### Legend and notes
    addControl(paste0("Preliminary data for ","<b>", gsub("_"," ",spp.all[i])),
               position = "topright") %>%
    addLegend(labels = c(
      "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
      "Estimated species distribution (50km buffers around occurrence points)",
      "Source localities for ex situ material (in botanic gardens/genebanks)",
      "Estimated area represented by ex situ material (50km buffers around source localities)",
      "Protected Areas (WDPA)",
      "Ecoregions (Olson et al. 2021)"),
      colors = c("#56116b","#965da7","#335707","#80bf30","#753a10","#dec171"),
      title = "Legend", position = "topright", opacity = 0.8) %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    # Layers control
    addLayersControl(
      overlayGroups = c(
        "Occurrence points (GBIF, iDigBio, IUCN Red List, herbaria consortia, FIA, BIEN)",
        "Estimated species distribution (50km buffers around occurrence points)",
        "Source localities for ex situ material (in botanic gardens/genebanks)",
        "Estimated area represented by ex situ material (50km buffers around source localities)",
        "Protected Areas (WDPA)",
        "Ecoregions (Olson et al. 2001)"),
      options = layersControlOptions(collapsed = FALSE),
      position = "bottomright") %>%
    addControl(
      "Toggle the checkboxes on and off to control which layers are visible...",
      position = "bottomright") %>%
    addControl(
      "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      position = "bottomleft") %>%
    setView(-98, 40, zoom = 5)
  return(map)
}
