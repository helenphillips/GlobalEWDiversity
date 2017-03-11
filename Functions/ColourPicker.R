
ColourPicker <- function(variable){
  
  if(any(levels(variable) %in% c("Broadleaf deciduous forest", "Herbaceous"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Broadleaf deciduous forest","Broadleaf evergreen forest","Cropland",                        
                                 "Herbaceous","Herbaceous with spare tree/shrub", "Needleleaf evergreen forest"),
                     colour =c("696346", "FFD300", "68089D", "103B9D", "0041D8", "1B1B1B"))
                      # colour = c("red", "green", "blue", "grey", "purple", "black"))
    colos <- tp$colour[which(levels(variable) == tp$habitat)]
    return(colos)
    }
  
  if(any(levels(variable) %in% c("Pasture", "Primary vegetation", "Production - Arable"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Pasture","Primary vegetation", "Production - Arable","Production - Crop plantations",
                                 "Production - Wood plantation", "Secondary vegetation"),
                     colour =c("696346", "FFD300", "68089D", "103B9D", "0041D8", "1B1B1B"))
    # colour = c("red", "green", "blue", "grey", "purple", "black"))
    colos <- tp$colour[which(levels(variable) == tp$habitat)]
    return(colos)
    
  }
  
  
  
}




