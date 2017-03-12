ColourPicker <- function(variable){
  
  if(any(levels(variable) %in% c("Broadleaf deciduous forest", "Herbaceous"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Broadleaf deciduous forest","Broadleaf evergreen forest","Cropland", 
                                 "Cropland/Other vegetation mosaic",
                                 "Herbaceous","Herbaceous with spare tree/shrub", "Needleleaf evergreen forest",
                                 "Paddy field", "Sparse vegetation", "Shrub", "Wetland", "Unknown/Other"),
                     colour =c("4CAF50", "1B5E20","E65100", "FF9800","6D4C41","A1887F","A5D6A7","9C27B0",
                                "F06292","26C6DA","0288D1","607D8B"))
    colos <- tp$colour[which(levels(variable) == tp$habitat)]
    return(colos)
    }
  
  if(any(levels(variable) %in% c("Pasture", "Primary vegetation", "Production - Arable"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Pasture","Primary vegetation", "Production - Arable","Production - Crop plantations",
                                 "Production - Wood plantation", "Secondary vegetation", "Unknown"),
                     colour =c("FFEB3B", "1B5E20", "E65100", "FF9800", "795548", "4CAF50", "607D8B"))
    colos <- tp$colour[which(levels(variable) == tp$habitat)]
    return(colos)
    
  }
  
  
  
}




