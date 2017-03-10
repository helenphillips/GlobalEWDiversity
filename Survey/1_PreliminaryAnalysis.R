########################################################
# 0. Set Working directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\Survey\\")
}

########################################################
# 1. Load libraries
########################################################

library(reshape)

########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

if(!dir.exists("EdittedData")){
  dir.create("EdittedData")
}
Data_out <- "EdittedData"

########################################################
# 3. Extracting most recent CSV
########################################################

files <- list.files(file.path("Data"))
file_dates <- sapply(strsplit(files, "_"), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
all <- read.csv(file.path(file.path("Data"), paste(date, "_Soil biodiversity survey short-report.csv", sep="")))

########################################################
# 4. Establishing relevant columns
########################################################

taxa <- which(names(all) %in% c("Bacteria", "Fungi", "Archaea", "Protists", "Rotifera", "Nematoda", "Mollusca", "Annelida","Tardigrada", "Arthropoda",
                 "Chelicerata.Arachnida", "Myriapoda", "Crustacea.Malacostracea", "Hexapoda.Ectognatha.Insecta", "Hexapoda.Entognatha"))

numSpp <- grep("How.many.species", names(all))                                                                                                                                                                
contribution <- grep("Would.you.be.willing.to.contribute.your.data.to.a.global", names(all))
geog <- grep("geographical.area", names(all))
collection <- which(names(all) %in% c("a.museum.collection", "a.personal.dataset", "data.compiled.from.publications.or.books"))
########################################################
# 5. Basic Stats
########################################################
nrow(all) ## 176 responses
table(all[,contribution]) 
#No answer   Maybe    No   Yes 
#15           39     8   114 

########################################################
## What sort of data
########################################################

collects <- all[,c(1, collection)]
c <- melt(collects, id = "X.")
c <- droplevels(c[which(c$value != ""),])
table(c$value)
# a museum collection                       a personal dataset 
# 43                                      131 
#data compiled from publications or books 
# 62 

########################################################
## What databases do you know of
########################################################

dbs <- all$What.databases.are.you.aware.of.that.contain.information.on.soil.fauna.diversity.
dbs <- dbs[dbs != ""]
dbs <- as.vector(dbs)
x <- strsplit(dbs, ",")
x <- unlist(x)
x <- strsplit(x, "\n")
x <- unlist(x)
x<- as.factor(x)

levels(x)[agrep("gbif", levels(x), ignore.case = TRUE)] <- "GBIF"
levels(x)[agrep("Edapobase", levels(x), ignore.case = TRUE)] <- "Edapobase"
levels(x)[agrep("senckenberg", levels(x), ignore.case = TRUE)] <- "None"
levels(x)[agrep("drilobase", levels(x), ignore.case = TRUE)] <- "Drilobase"
levels(x)[agrep("Macrofauna", levels(x), ignore.case = TRUE)] <- "Macrofauna"
levels(x)[agrep("LTER", levels(x), ignore.case = TRUE)] <- "LTER"
levels(x)[agrep("BETSI", levels(x), ignore.case = TRUE)] <- "BETSI"
levels(x)[agrep("www.annelida.net", levels(x), ignore.case = TRUE)] <- "annelida.net"
levels(x)[agrep("earthmicrobiome", levels(x), ignore.case = TRUE)] <- "Earth Microbiome project"
levels(x)[grep("Bold", levels(x), ignore.case = TRUE)] <- "Bold/NCBI/Genbank"
levels(x)[grep("NCBI", levels(x), ignore.case = TRUE)] <- "Bold/NCBI/Genbank"
levels(x)[agrep("GenBank", levels(x), ignore.case = TRUE)] <- "Bold/NCBI/Genbank"
levels(x)[agrep("predicts", levels(x), ignore.case = TRUE)] <- "PREDICTS project"
levels(x)[grep("fao", levels(x), ignore.case = TRUE)] <- "FAO"
levels(x)[grep("Non", levels(x), ignore.case = TRUE)] <- "None"
levels(x)[grep("not aware", levels(x), ignore.case = TRUE)] <- "None"
levels(x)[grep("do not know", levels(x), ignore.case = TRUE)] <- "None"
levels(x)[grep("don't know", levels(x), ignore.case = TRUE)] <- "None"

levels(x)[grep("I have a", levels(x), ignore.case = TRUE)] <- "Personal dataset"

x_sorted <- sort(table(x),decreasing=T)
x_sorted <- x_sorted[x_sorted > 1]

x_sorted <- x_sorted[names(x_sorted) %in% c("None","GBIF","Edapobase","Bold/NCBI/Genbank", "Drilobase",                                                
"BETSI", "Macrofauna", "Earth Microbiome project", "PREDICTS project", "FAO", "annelida.net")]


pdf(file = file.path(figures, "Databases.pdf")) ## This then was converted online to be used in teh presentation
par(mar=c(10.5, 3, 2, 1))
barplot(x_sorted, las = 2)
dev.off()
########################################################
## Which taxa and where
########################################################
taxa_recorded <- all[,c(1, taxa, geog)]
m <- (melt(taxa_recorded, id = c("X.", "What.geographical.area.is.your.data.from..e.g...country.or.region.within.a.country..")))
levels(m$variable)
levels(m$value) <- gsub("\\/",".", levels(m$value))

m <- droplevels(m[which(m$value != ""),])
names(m)[2] <- "Region"
n_taxa <- data.frame(table(m$value))

write.csv(m, file= file.path(Data_out, "TaxaByRegion_toEdit.csv"))



unique(m$Region)
levels(m$Region)[levels(m$Region) == "nl"] <- "Netherlands"                                                                                                                                                                           
levels(m$Region)[levels(m$Region) == "AMERICAS and Caribbean"] <- "Caribbean"                                                                                                                                                      
levels(m$Region)[levels(m$Region) ==  "france"] <- "France"                                                                                                                                                                        
levels(m$Region)[levels(m$Region) == "Cape Town- South Africa"] <- "Cape Town, South Africa"                                                                                                                                                      
levels(m$Region)[levels(m$Region) == "UK, Ireland"] <- "UK"                                                                                                                                                               
levels(m$Region)[levels(m$Region) == "Brazil; Amazon rainforest."] <- "Amazon Rainforest"                                                                                                                                                    
levels(m$Region)[levels(m$Region) =="New York City metropolitan area"] <- "New York City"                                                                                                                                              
levels(m$Region)[levels(m$Region) == "Veluwe, province of Gelderland, Netherlands"] <- "Gelderland"                                                                                                                                   
levels(m$Region)[levels(m$Region) == "mostly within the Miami-Dade County, FL region"] <- "Miami-Dade County"                                                                                                                         
levels(m$Region)[levels(m$Region) == "63"] <- "NA"                                                                                                                                                                       
levels(m$Region)[levels(m$Region) == "tropical (India)"] <- "India"                                                                                                                                                             
levels(m$Region)[levels(m$Region) == "Switzerland, glacier forelands"] <- "Switzerland"                                                                                                                                              
levels(m$Region)[levels(m$Region) == "Agroecosystem"] <- "NA"                                                                                                                                                               
levels(m$Region)[levels(m$Region) == "All regions in Switzerland."] <- "Switzerland"                                                                                                                                                   
levels(m$Region)[levels(m$Region) == "Northern Alaska, USA"] <- "Northern Alaska"                                                                                                                                                          
levels(m$Region)[levels(m$Region) == "eastwestern China"] <- "China"                                                                                                                                                        
levels(m$Region)[levels(m$Region) == "UK, and a few records from Ireland (Tom Bolger of UCD in ireland is compiling new species lists and recording data for Irish acari - worth contacting him!)."] <- "UK"                 
levels(m$Region)[levels(m$Region) == "NE US"] <- "USA"                                                                                                                                                                        
levels(m$Region)[levels(m$Region) == "United States of America, there are 47 sites distributed across the U.S.A."] <- "USA"                                                                                                    
levels(m$Region)[levels(m$Region) == "United states of America-Eastern California-Great Basin floristic province"] <- "Great Basin"                                                                                                   
levels(m$Region)[levels(m$Region) == "Semiarid zones. Spain."] <- "Spain"                                                                                                                                                                                                                                                                                                               
levels(m$Region)[levels(m$Region) == "High alpine Subarctic Sweden"] <- "Sweden"                                                                                                                                               
levels(m$Region)[levels(m$Region) == "Note that the arthropod data is from soil extractions. Most of this data comes from the Southeastern US."] <- "USA"                                                                      
levels(m$Region)[levels(m$Region) == "Puerto Rico, forests and pastures"] <- "Puerto Rico"                                                                                                                                            
levels(m$Region)[levels(m$Region) == "South Africa and limited African records"] <- "South Africa"                                                                                                                                      
levels(m$Region)[levels(m$Region) == "These are all data from the Netherlands, most of them from the 1980s up till the present date"] <- "Netherlands"                                                                                 
levels(m$Region)[levels(m$Region) == "Scotland/UK"] <- "Scotland"                                                                                                                                                                   
levels(m$Region)[levels(m$Region) == "mostly Wales but includes samples from other areas of Europe"] <- "Wales"                                                                                                                  
levels(m$Region)[levels(m$Region) == "france, région parisienne"] <- "Paris, France"                                                                                                                                                 
levels(m$Region)[levels(m$Region) == "global"] <- "Global"                                                                                                                                                                        
                                                                                                                      
levels(m$Region)[levels(m$Region) == "grassland ecosystem of northern China"] <- "Northern China"                                                                                                                                        
levels(m$Region)[levels(m$Region) =="East-Algeria (Algeria)"] <- "Algeria"                                                                                                                                                       
levels(m$Region)[levels(m$Region) =="North americal boreal and temperate forest, peatlands"] <- "USA"                                                                                                                        
[66] "Nematodes - UK wide; Earthworms - Scotland"                                                                                                                                    
levels(m$Region)[levels(m$Region) =="Worldwide"] <- "Global"                                                                                                                                                                    
levels(m$Region)[levels(m$Region) =="primarily Europe, focus on Germany; however some data worldwide"] <- "Germany"                                                                                                               
levels(m$Region)[levels(m$Region) =="Mostly Greece. The collection extends to the Palearctic but contains many specimens from the tropics through donations and acquisitions."] <- "Greece"                                     
levels(m$Region)[levels(m$Region) =="world wide"] <- "Global"                                                                                                                                                                  
                                                                                                                                                             
levels(m$Region)[levels(m$Region) =="region within a country (France)"] <- "France"                                                                                                                                              
levels(m$Region)[levels(m$Region) =="For my  personal dataset: France. For the Betsi db, please see our web site"] <- "France"                                                                                                   
levels(m$Region)[levels(m$Region) =="southeastern US, southwestern US"] <- "USA"                                                                                                                                              
levels(m$Region)[levels(m$Region) =="world-wide, greatest strenght is eastern North America"] <- "Global"                                                                                                                        
levels(m$Region)[levels(m$Region) =="agroecosystem in Tarragona, Catalonia (SE Spain)"] <- "Spain"                                                                                                                              
levels(m$Region)[levels(m$Region) =="Brazil Eastern Amazon, Maranhão and Pará States"] <- "Brazil"                                                                                                                             
levels(m$Region)[levels(m$Region) == "Western arid land of India"] <- "India"                                                                                                                                               
[87] "Stubai Valley, Tyrol, Austria (Central Alps), Matsch Valley, South Tyrol, Italy (Central Alps), Madritsch Valley, South Tyrol, Italy (Southern Alps)"                          
[88] "Bolivia, California, Colombia, Colorado, Ecuador, El Salvador, Nicaragua, Peru"                                                                                                
[89] "More French records but also few records abroad (Italy, Spain, Austria, Romania...)"                                                                                           
levels(m$Region)[levels(m$Region) =="Germany (North)"] <- "Germany"                                                                                                                                                              
[91] "Vinschgau (Northern Italy), Tyrol (Austria)"                                                                                                                                   
[92] "Russia, Komi Republic, Nenets Autonomous District (european north-east)"                                                                                                       
levels(m$Region)[levels(m$Region) == "continental france"] <- "France"                                                                                                                                                            
levels(m$Region)[levels(m$Region) =="Finland, from southern parts of the country to the Arctic Circle."] <- "Finland"                                                                                                             
levels(m$Region)[levels(m$Region) =="Russia: Urals and Siberia (mainly)"] <- "Russia"                                                                                                                                            
[96]  "North Caucasus"                                                                                                                                                              
levels(m$Region)[levels(m$Region) =="Mostly Brazil (especially Southern and Southeastern regions), although I have a few records of worms from Bolivia, Peru, Venezuela, Argentina, Uruguay, Bolivia, and Slovenia."] <- "Brazil"
levels(m$Region)[levels(m$Region) == "Pacific Northwest, WA, ID, OR"] <- "Pacific Northwest"                                                                                                                                                 
levels(m$Region)[levels(m$Region) =="Coshocton OH, Piketon OH, Columbus OH (all USA)"] <- "Ohio"                                                                                                                               
levels(m$Region)[levels(m$Region) == "UK, north Yorkshire"] <- "Yorkshire"                                                                                                                                                           
levels(m$Region)[levels(m$Region) =="SE Kansas, USA"] <- "Kansas"                                                                                                                                                                
[104] "Badlands located in Granada, Spain (Mediterranean)"                                                                                                                            
levels(m$Region)[levels(m$Region) == "mainly around Avignon (SE of France)"] <- "Avignon"                                                                                                                                         
[107] "Russia, Republic of Komi, village Vodny"                                                                                                                                       
levels(m$Region)[levels(m$Region) == "most from Denmark, but also Holland, Wales"] <- "Denmark"                                                                                                                                    
[110] "Iran, Central Zagros and Central Alborz regions."                                                                                                                              
levels(m$Region)[levels(m$Region) == "N-America"] <- "USA"                                                                                                                                                                     
levels(m$Region)[levels(m$Region) == "SE USA.  Mostly Texas.  Also Oklahoma and Mississippi."] <- "Texas"                                                                                                                        
levels(m$Region)[levels(m$Region) == "All coninents except Antarctica"] <- "Global"                                                                                                                                               
[114] "Europe, south America, south east Asia"                                                                                                                                        
levels(m$Region)[levels(m$Region) == "none"] <- "NA"                                                                                                                                                                          
levels(m$Region)[levels(m$Region) == "The species we collected were mainly in South China."] <- "South China"                                                                                                                          
levels(m$Region)[levels(m$Region) =="Australia, mostly southern regions, mostly agricultural soils, but also urban"] <- "Australia"                                                                                                 
levels(m$Region)[levels(m$Region) == "France, �Zle-de-France"] <- "France"                                                                                                                                                        
[121] "Canada, USA, Costa Rica, some from Australasia"                                                                                                                                
levels(m$Region)[levels(m$Region) == "korea south"] <- "South Korea"                                                                                                                                                                   
[123] "Scandinavia, Australia, New Zealand, Antarctica"                                                                                                                               
[124] "Russia, Komi Republic, european North-East of Russia, the Urals mountains"                                                                                                     
[125] "Georgia, Caucasus"                                                                                                                                                             
[126] "Mediterranean basin (Spain), Great plains (Colorado, USA), Dry Tropics (Nicaragua)"                                                                                            
[128] "Canada and Alaska mainly, followed by rest of USA, Costa Rica, and a scattering from other parts of theworld"                                                                  
[129] "Australia, Solomon Islands"                                                                                                                                                    
[130] "Eastern United States"                                                                                                                                                         
[131] "Neotropical (mostly Costa Rica, Mexico, Guatemala, Honduras, Cuba, Nicaragua"                                                                                                  
levels(m$Region)[levels(m$Region) == "Israel - The Negev desert and the Arava Valley. including semi-arid, arid and hyper-arid"] <- "Israel"                                                                                      
[133] "Former Soviet Union"                                                                                                                                                          
levels(m$Region)[levels(m$Region) =="Brazil, Northeast region."] <- "Brazil"                                                                                                                                                    
levels(m$Region)[levels(m$Region) =="Russia, European North-East of Russia, especially Komi Republic"] <- "Russia"                                                                                                               
levels(m$Region)[levels(m$Region) == "All Brazil"] <- "Brazil"                                                                                                                                                                    
levels(m$Region)[levels(m$Region) == "UK, West Midlands"] <- "West Midlands"                                                                                                                                                             
[139] "New York State"                                                                                                                                                                
[140] "Europe and French Guyana"                                                                                                                                                      
[141] "SE Spain"                                                                                                                                                                      
[142] "Svalbard, Holland and Sweden"                                                                                                                                                  
[143] "Scandinavia (Norway and Sweden), Svalbard, Australia, New Zealand, Antarctica, Falkland Islands"                                                                               
[144] "Canada: boreal forest (Quebec, Ontario), subarctic taiga (Quebec), boreal peatlands (Ontario), coastal temperate and montane forests (British Columbia)"                       
levels(m$Region)[levels(m$Region) =="Netherlands only"] <- "Netherlands"                                                                                                                                                              
levels(m$Region)[levels(m$Region) == "The Netherlands, salt marsh"] <- "Netherlands"                                                                                                                                                   
levels(m$Region)[levels(m$Region) =="My data is mostly from Sweden (Southeen Sweden and area around Uppsala), the UK (transect study) and Svalbard"] <- "Sweden"                                                                
[149] "Southern India ( western ghats)" 