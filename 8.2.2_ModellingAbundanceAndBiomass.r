## SET WORKING DIRECTORY ---------------------------
setwd("C:/restore2/hp39wasi/sWorm/EarthwormAnalysis")


## VARIABLES AND LIBRARIES -------------------------
data_in <- "8_Data"

# load(file=file.path(data_in, "AbundanceBiomasDT.RData"))
library('LaF')
library(biglm)
library(stringr)

huge_file <- file.path(data_in, "AbundanceBiomasDT.csv")


## CONNECTION --------------------------------------
#First detect a data model for your file:
model <- detect_dm_csv(huge_file, sep=" ", header=TRUE)

df.laf <- laf_open(model)

## BIGLM VARS-----------------------------------------

nlines <- determine_nlines(huge_file)

## How many blocks would their be?
## So I can have a counter on my loop

nRows <- 1e6

counter <- ceiling(nlines / nRows)
counterseq <- 2:counter # because first chunk is outside loop

whichLine <- 1


## START PROCESS ---------------------------------
goto(df.laf, whichLine) # Start at the beginning

div <- next_block(df.laf,nrows=1e6)
div <- as.data.frame(str_split_fixed(div$abundance.biomass, ",", 2))
div$abundance <- as.numeric(as.character(as.vector(div[,1])))
div$biomass <- as.numeric(as.character(as.vector(div[,2])))
div[,c(1:2)] <- NULL


ff<-abundance ~ biomass

lm1 <- biglm(ff, div)


## LOOP -------------------------------------------
whichLine <- whichLine + nRows
goto(df.laf, whichLine)

for(c in counterseq){
  print(c)
  div <- next_block(df.laf,nrows=1e6)
 
  div <- as.data.frame(str_split_fixed(div$abundance.biomass, ",", 2))
  div$abundance <- as.numeric(as.character(as.vector(div[,1])))
  div$biomass <- as.numeric(as.character(as.vector(div[,2])))
  div[,c(1:2)] <- NULL
  
  
  lm1 <- update(lm1,div)
  print(summary(lm1))
  whichLine <- whichLine + nRows
  goto(df.laf, whichLine)
}

save(lm1, file =file.path(data_in, "AbundanceByBiomassBIGLM.rds"))
