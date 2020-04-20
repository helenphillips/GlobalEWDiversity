## THIS SCRIPT MAKES THE DATASET THAT WILL BECOME OPEN-ACCESS


## VARIABLES --------------
data_in <-"9_Data"


if(!dir.exists("22_Data")){
  dir.create("22_Data")
}
data_out <- "22_Data"



## WORKING DIRECTORY ------

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "TSGIS02"){
  setwd("C:/sWorm/EarthwormAnalysis")
}

## LOAD DATA ------------
dat <- read.csv(file.path(data_in, "sWorm_CompleteDataSet_correction.csv"))


## REMOVE COLUMNS -------

notNeeded <- c("ID",                                              
"Observational","Latitude__decimal_degrees",   
"Longitude__Decimal_Degrees","Altitude__m",               
"Country","Sample_StartDate_Month",   
"Sample_StartDate_Year","Sample_EndDate_Month",     
"Sample_EndDate_Year","ExtractionMethod", "CEC","CEC_unit", "CEC_mean",      "Soil_Organic_Matter__percent",
"SOM_mean",
"Sampled_Area","Sampled_Area_Unit",       
"Sample_Effort","Base_Saturation_percent",    
"BaseSaturation_mean","WRB_FAO_SoilType",         
"LandUse",     "biome",             "Soil_Moisture_percent",                   
"IPBES_Habitat_Units","Management_System",   "USDA_SoilTexture",    
"Tillage","Pesticide",               
"Fertilizer","Selectively_harvested",     
"Clear_cut","Fire",                 
"Stocking_rate","Grazing_all_year",        
"Rotation","Monoculture",      "HabitatCover",      
"Planted","Habitat_as_described",     
"bio10_2","bio10_3",                   
"bio10_5",               
"bio10_6", 
"bio10_8","bio10_9",           
"bio10_10","bio10_11",              
"bio10_13",         
"bio10_14",      "Study_site", 
"bio10_16","bio10_17",                  
"bio10_18","bio10_19",                 
"TAXNWRB_1",    "ScalePET","ScalePETSD",     "scaleAridity",             
"country","LU_Mgmt",                  
"scalePH")


dat <- dat[,-(which(names(dat) %in% notNeeded))]

## RE-ORDER COLUMNS -------


toOrder <- c("file","Study_Name",                  
"Site_Name",                  
"PH",                          
"PH_Collection_Method","PH_mean",                     
"Organic_Carbon__percent",     
"OC_mean","C.N_ratio",                   
"CN_mean","Sand__percent",               
"Silt__percent","Clay__percent",               
"sand_silt_clay_mean" ,                    
"ph_new",
"PHIHOX","CLYPPT",                      
"SLTPPT","SNDPPT",                      
"CECSOL","ORCDRC",                      
"phFinal","ClayFinal",                   
"SandFinal","SiltFinal",                   
"OCFinal",
"ESA",
"bio10_1",                     
"bio10_4","bio10_7",                     
"bio10_12","bio10_15",                    
"elevation","SnowMonths",  "SnowMonths_cat",                
"Aridity",            
"PETyr",                    
"PET_SD",                
"SpeciesRichness","SpeciesRichnessUnit",         
"Site_WetBiomass","Site_WetBiomassUnits",        
"Site_Biomassm2","logBiomass",
"Site_Abundance","Site_AbundanceUnits" ,        
"Sites_Abundancem2","logAbundance")

dat <- dat[,match(toOrder, names(dat))]

write.csv(dat, file = file.path(data_out, "sWormModelData_correction.csv"))


## CREATE META-DATA ------------


df1 <- data.frame(colName = names(dat))

write.csv(df1, file = file.path(data_out, "sWormModelData_meta-data_correction.csv"), row.names = FALSE)

