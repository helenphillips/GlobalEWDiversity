ColourPicker <- function(variable){
  
  
  
  if(any(levels(variable) %in% c("Broadleaf deciduous forest", "Herbaceous"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Production - Herbaceous", "Production - Plantation", "Cropland/Other vegetation mosaic", 
                                "Urban", "Bare area (consolidated", "Bare area (unconsolidated", 
                                "Water bodies", "Mixed forest", "Tree open", "Herbaceous with spare tree/shrub", 
                                "Shrub", "Herbaceous", "Sparse vegetation", "Broadleaf evergreen forest", 
                                "Broadleaf deciduous forest", "Needleleaf evergreen forest", "Needleleaf deciduous forest"),
                    colour = c("ffff64", "ffff00", "c8c864", "c31400", "dcdcdc", 'fff5d7',
                                "0046c8", "788200", "8ca000", "be9600", "966400", "ffb432", "ffebaf",
                                "286400", "00a000", "003c00", "285000"))

    }
  
  
  if(any(levels(variable) %in% c("Pasture", "Primary vegetation", "Production - Arable"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Pasture","Primary vegetation", "Production - Arable","Production - Crop plantations",
                                 "Production - Wood plantation", "Secondary vegetation", "Urban", "Unknown"),
                     colour =c("FFEB3B", "1B5E20", "E65100", "FF9800", "795548", "4CAF50", "788084", "607D8B"))
    
  }
  
  if(any(levels(variable) %in% c("Annual crop", "Perennial crops","Integrated systems"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Primary vegetation","Secondary vegetation","Annual crop",            
                            "Perennial crops","Integrated systems","Tree plantations",       
                                "Pastures (grazed lands)","Urban", "Unknown"),
                     colour = c("006400", "008B00", "EEEE00",
                                "8B8B00", "8B8970", "8B2323", 
                                "CD8500","8B8378", "292421"))
    }
  
  colos <- tp$colour[match(levels(variable), tp$habitat)]
  
  return(colos)
}




