## This function takes in two the dataframes from individual studies
## and formats the columns etc.

formatSites <- function(sites){
  
  names(sites) <- gsub("\\s+","_", names(sites))
  names(sites) <- gsub("\\(","_", names(sites))
  names(sites) <- gsub("\\)","", names(sites))
  names(sites) <- gsub("\\?","", names(sites))
  names(sites) <- gsub("\\%","percent", names(sites))
  names(sites)[which(names(sites) == "C/N_ratio")] <- "C.N_ratio"
  names(sites)[which(names(sites) == "WRB/FAO_SoilType")] <- "WRB_FAO_SoilType"
  

  sites$Study_Name <- as.factor(sites$Study_Name)
  sites$Site_Name <- as.factor(sites$Site_Name)
  sites$Observational <- as.factor(sites$Observational)
  sites$Latitude__decimal_degrees <- as.numeric(sites$Latitude__decimal_degrees)
  sites$Longitude__Decimal_Degrees <- as.numeric(sites$Longitude__Decimal_Degrees)
  sites$Altitude__m <- as.numeric(sites$Altitude__m)
  sites$Country <- as.factor(sites$Country)
  sites$Sample_StartDate_Month <- as.integer(sites$Sample_StartDate_Month)
  sites$Sample_StartDate_Year <- as.integer(sites$Sample_StartDate_Year)
  sites$Sample_EndDate_Month <- as.integer(sites$Sample_EndDate_Month)
  sites$Sample_EndDate_Year <- as.integer(sites$Sample_EndDate_Year)
  sites$ExtractionMethod <- as.factor(sites$ExtractionMethod)
  sites$Sampled_Area <- as.integer(sites$Sampled_Area)
  sites$Sampled_Area_Unit <- as.factor(sites$Sampled_Area_Unit)
  sites$Sample_Effort <- as.numeric(sites$Sample_Effort)
  sites$PH <- as.numeric(sites$PH)
  sites$PH_Collection_Method <- as.factor(sites$PH_Collection_Method)
  sites$CEC <- as.numeric(sites$CEC)
  sites$Base_Saturation_percent <- as.numeric(sites$Base_Saturation_percent)
  sites$Organic_Carbon__percent <- as.numeric(sites$Organic_Carbon__percent)
  sites$Soil_Organic_Matter__percent <- as.numeric(sites$Soil_Organic_Matter__percent)
  sites$C.N_ratio <- as.numeric(sites$C.N_ratio)                                 
  sites$Sand__percent <- as.numeric(sites$Sand__percent)                                  
  sites$Silt__percent <- as.numeric(sites$Silt__percent)
  sites$Clay__percent <- as.numeric(sites$Clay__percent)
  sites$USDA_SoilTexture <- as.factor(sites$USDA_SoilTexture)
  sites$Soil_Moisture_percent <- as.numeric(sites$Soil_Moisture_percent)
  sites$WRB_FAO_SoilType <- as.factor(sites$WRB_FAO_SoilType)                                 
  sites$LandUse <- as.factor(sites$LandUse)                              
  sites$HabitatCover <- as.factor(sites$HabitatCover)
  sites$IPBES_Habitat_Units <- as.factor(sites$IPBES_Habitat_Units)
  sites$Management_System <- as.factor(sites$Management_System)
  sites$Tillage <- as.integer(sites$Tillage)
  sites$Pesticide<- as.integer(sites$Pesticide)
  sites$Fertilizer <- as.integer(sites$Fertilizer)
  sites$Selectively_harvested <- as.integer(sites$Selectively_harvested)
  sites$Clear_cut <- as.integer(sites$Clear_cut)
  sites$Fire <- as.integer(sites$Fire)
  sites$Stocking_rate <- as.integer(sites$Stocking_rate)
  sites$Grazing_all_year <- as.integer(sites$Grazing_all_year)
  sites$Rotation <- as.integer(sites$Rotation)
  sites$Monoculture <- as.integer(sites$Monoculture)
  sites$Planted <- as.integer(sites$Planted)
  
  sites$SpeciesRichness <- as.numeric(sites$SpeciesRichness)                          
  sites$SpeciesRichnessUnit <- as.factor(sites$SpeciesRichnessUnit)
  sites$Habitat_as_described <- as.factor(sites$Habitat_as_described)              
  sites$file <- as.factor(sites$file)
  sites$WetBiomass <- as.numeric(sites$WetBiomass)
  sites$WetBiomassUnits <- as.factor(sites$WetBiomassUnits)
  sites$Abundance <- as.numeric(sites$Abundance)
  sites$Abundance_Units <- as.factor(sites$Abundance_Units)
  
  sites$Study_site <- as.factor(paste(sites$Study_Name, sites$Site_Name))
  
  names(sites)[which(names(sites) == "WetBiomass")] <- "Site_WetBiomass"
  names(sites)[which(names(sites) == "WetBiomassUnits")] <- "Site_WetBiomassUnits"

  names(sites)[which(names(sites) == "Abundance")] <- "Site_Abundance"
  names(sites)[which(names(sites) == "Abundance_Units")] <- "Site_AbundanceUnits"

  return(sites)
}

formatSpecies <- function(species){
  
  names(species) <- gsub("\\s+","_", names(species))
  names(species) <- gsub("\\(","_", names(species))
  names(species) <- gsub("\\)","", names(species))
  names(species)[which(names(species) == "Native/Non-native")] <- "Native.Nonnative"

  species$Study_ID <- as.factor(species$Study_ID)
  species$Site_Name <- as.factor(species$Site_Name)
  species$SpeciesBinomial <- as.factor(species$SpeciesBinomial)
  species$MorphospeciesID <- as.factor(species$MorphospeciesID)
  species$Genus <- as.factor(species$Genus)
  species$Family <- as.factor(species$Family)
  species$LifeStage <- as.factor(species$LifeStage)
  species$Native.Nonnative <- as.factor(species$Native.Nonnative)
  species$Functional_Type <- as.factor(species$Functional_Type)
  species$Abundance <- as.numeric(species$Abundance)
  species$Abundance_Unit <- as.factor(species$Abundance_Unit)
  species$file <- as.factor(species$file)
  species$WetBiomass <- as.numeric(species$WetBiomass)
  species$WetBiomassUnits <- as.factor(species$WetBiomassUnits)
  
  species$Study_site <- as.factor(paste(species$Study_ID, species$Site_Name))
  
  return(species)
}

SiteLevels <- function(sites){
  
  sites$Management_System <- factor(sites$Management_System, levels = c("None",  "Annual crop", "Perennial crops", "Integrated systems",                   
                                                                        "Tree plantations", "Pastures (grazed lands)","Unknown"))
  
  sites$LandUse <- factor(sites$LandUse, levels = c("Primary vegetation", "Secondary vegetation", "Pasture" ,
                                                    "Production - Arable", "Production - Crop plantations", 
                                                    "Production - Wood plantation",
                                                    "Urban", "Unknown"))
  
  sites$HabitatCover <- factor(sites$HabitatCover, levels = c("Broadleaf deciduous forest", "Broadleaf evergreen forest",
                                                              "Needleleaf deciduous forest","Needleleaf evergreen forest",
                                                              "Mixed forest", "Tree open",
                                                              "Cropland","Cropland/Other vegetation mosaic",
                                                              "Herbaceous", "Herbaceous with spare tree/shrub",
                                                              "Shrub", "Sparse vegetation",
                                                              "Urban","Bare area (consolidated",
                                                              "Paddy field","Wetland", "Water bodies", "Unknown"))
  
  sites$LU_Mgmt <- factor(sites$LU_Mgmt, levels = c( "Primary vegetation", "Secondary vegetation", "Pastures (grazed lands)", "Annual crop", 
                                                     "Perennial crops", "Integrated systems", "Tree plantations", 
                                                     "Urban", "Unknown" ))
  
  if("intensity" %in% names(sites)){sites$intensity <- as.factor(sites$intensity)}
  
  return(sites)
}
