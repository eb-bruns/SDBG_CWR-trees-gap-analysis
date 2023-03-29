

filter.points <- function(occ_data, manual_edits=NULL){
  
  ## first make sure all the T/F columns are logical type
  occ_data <- occ_data %>%
    dplyr::mutate(.cen = as.logical(.cen),
                  .inst = as.logical(.inst),
                  .outl = as.logical(.outl),
                  .con = as.logical(.con),
                  .urb = as.logical(.urb),
                  .yr1950 = as.logical(.yr1950),
                  .yr1980 = as.logical(.yr1980),
                  .yrna = as.logical(.yrna),
                  .nativectry = as.logical(.nativectry))
  
  ## filter occurrence data based on filter columns created in
  ##  5-refine-occurrence-data.R
  occ_fltr <- occ_data %>%
    dplyr::filter(
      # choose to keep ex situ points even if flagged (and manually remove
      #   any that are bad) or comment out next line and filter along with
      #   points from other databases
      database == "Ex_situ" |
        # choose which filters you want to use (comment out those you don't want)
        (.cen & .inst & .outl &
           #.con & .urb & .yr1950 & .yr1980 & .yrna &
           (.nativectry | is.na(.nativectry)) &
           basisOfRecord != "FOSSIL_SPECIMEN" &
           basisOfRecord != "LIVING_SPECIMEN" &
           establishmentMeans != "INTRODUCED" &
           establishmentMeans != "MANAGED" &
           establishmentMeans != "INVASIVE" &
           establishmentMeans != "CULTIVATED"
        ))
  print(paste0("Removed ", nrow(occ_data)-nrow(occ_fltr)," of ", nrow(occ_data),
               " points based on flagging colums"))
  
  ## (optionally) add additional edits for specific taxa
  # for example, here we add the geographic outliers back in for Carya ovata:
  add <- occ_data %>%
    dplyr::filter(taxon_name_accepted == "Carya ovata") %>%
    dplyr::filter(database == "Ex_situ" |
                    (.cen & .inst &
                       (.nativectry | is.na(.nativectry)) &
                       basisOfRecord != "FOSSIL_SPECIMEN" &
                       basisOfRecord != "LIVING_SPECIMEN" &
                       establishmentMeans != "INTRODUCED" &
                       establishmentMeans != "MANAGED" &
                       establishmentMeans != "INVASIVE" &
                       establishmentMeans != "CULTIVATED"))
  occ_fltr <- occ_fltr %>% dplyr::filter(taxon_name_accepted != "Carya ovata")
  occ_fltr <- suppressMessages(full_join(occ_fltr,add))
  
  # another example, we remove anything that is Prunus glandulosa (a bad
  #  synonym):
  occ_fltr <- occ_fltr %>%
    dplyr::filter(taxon_name != "Prunus glandulosa")
  
  ## (optionally) check document with manual point edits to see if anything
  ##  needs to be removed or added back in
  if(!missing(manual_edits)){
    for(i in 1:length(unique(occ_data$taxon_name_accepted))){
      target_taxon <- unique(occ_data$taxon_name_accepted)[i]
      print(paste0("Checking manual edits for ",target_taxon))
      taxon_edits <- manual_edits[which(
        manual_edits$taxon_name_accepted == target_taxon),]
      # remove if in bounding box
      # bounding box is in format: lat-max, long-min, lat-min, long-max
      #  (i.e. the SE US would be something like 34, -101, 25, -77)
      if(!is.na(taxon_edits$bounding_box)){
        bounds <- unlist(strsplit(taxon_edits$bounding_box,"; "))
        for(j in 1:length(bounds)){
          within <- unlist(strsplit(bounds[j],", "))
          remove <- occ_fltr %>%
            dplyr::filter(
              taxon_name_accepted == target_taxon &
                (decimalLatitude < as.numeric(within[1]) &
                   decimalLongitude > as.numeric(within[2]) &
                   decimalLatitude > as.numeric(within[3]) &
                   decimalLongitude < as.numeric(within[4])))
          occ_fltr <- occ_fltr %>%
            dplyr::filter(!(UID %in% unique(remove$UID)))
          print(paste0("--Removed ", nrow(remove),
                       " points based on bounding box"))
        }
      }
      # remove if ID listed in remove_id
      if(!is.na(taxon_edits$remove_id)){
        remove <- unlist(strsplit(taxon_edits$remove_id,"; "))
        occ_fltr <- occ_fltr %>% dplyr::filter(!(UID %in% remove))
        print(paste0("--Removed ", length(remove),
                     " points based on IDs to remove"))
      }
      # add back if ID listed in keep_id
      if(!is.na(taxon_edits$keep)){
        keep <- unlist(strsplit(taxon_edits$keep,"; "))
        add <- occ_fltr %>% dplyr::filter(UID %in% keep)
        occ_fltr <- suppressMessages(full_join(occ_fltr,add))
        print(paste0("--Added back ", length(keep),
                     " points based on IDs to keep"))
      }
    }
  }
  
  return(occ_fltr)
  
}
