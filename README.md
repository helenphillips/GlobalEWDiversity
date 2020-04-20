<<<<<<< HEAD

# Global distribution of earthworm diversity

*Phillips et al., 2019 Science*

This work almalgamates data from various sources (i.e., from data published in the literature or that has been sent by 'data providers').

I have tried to write this code so that it works on anyone's machine....but at various points, this will not work very well/at all. Sorry. 

In theory, the input for one script is the output from the previous (or a previous) script. All scripts create an output folder, for data, and place all figures into a separate folder.

Please feel free to raise any issues.

## Files

### 0_GetData.R
This script accesses my personal Google Drive to download all the files that contain the data. This will not work for anyone.

Once all the data has been downloaded, the first stage of cleaning is undertaken (i.e., ensuring all columns are the correct class). Three files are created: a bibliography, a site-level dataframe and a species-level dataframe.

For the site-level dataframe, the site-level metrics (species richness, abundance and biomass) are calculated. 

### 1_AddCHELSA.R
Based on the coordinates of each site in the site-level metrics, the relevant CHELSA data (climatologiesl LINK) are appended. The CHELSA data had been downloaded previously.

### 2_AddSoilGrids.R
Based on the coordinates of each site in the site-level metrics, the relevant SoilGrids data (LINK) are appended. The SoilGrids data had been downloaded previously.

### 3_AddOtherVariables.R
Based on the coordinates of each site in the site-level metrics, data from other global data layers are appended. For now, please see the manuscript for details on the other datalayers.

### 4_SpeciesNames.R 
This creates a dataframe of every species with an ID, along with some additional data, such as who provided the data, so that our earthworm experts could harmonise the names

### 5_DataCheck.R
Most importantly, this script converts the biomass and abundance values to a common unit. Renames factor levels, and does some basic checks

### 6_DataExploration.R
Basic exploration of the data. Also converts pH values to a common unit, snow cover to a categorical, and log's response variables.

### 7_MeasuredversusSoilGrids.R
This script investigates whether it is appropriate to use a model that contains a mixture of site-level sampled soil properties and the soil properties from SoilGrids.

### 8.1_ModellingRichness.R  & 8.1_submitScript_ModellingRichness.sh
Script for the site-level model of species richness. As this was run on the HPC, the submit script is also present.

### 8_Modelling.R
Script for all the site-level models (species richness, abundance and biomas).

### 9_MainModelsFigures.R
Script that produces some figures from the main site-level models

### 10_DataForGlobalLayers.R
Creates a dataframe of the just the predictive data that was used in the models. 

### 10.x_MapCoefficients_xxx.R & 10.x_submitScript_MapCoefficients_XX.sh
Scripts that take the relevant site-level model, and creates a raster of the predicted values globally (run on the HPC, so submit scripts included). For this, Carlos had split the underlying global data layers into regions, so the scripts/models run for a single region. 

### 11_CreateMap.R
Using the predicted regions of the globe, the final maps are assembled and saved

### 11.1_ValuesFromMap.R
From the predicted values, calculating the means, SDs etc.

### 12_MainModelsVariableImportance.R
Using randomForest models, the three main models are reconstructed, and the most important variables are identified

### 13_CrossValidationMainModels.R
For the three main models, 10-fold cross validation is performed

### 14_CrossValidationSoilGrids.R
The three site-level models (species richness, abundance and biomass) are reconstructed, but using only the SoilGrids soil properties data. 10-fold cross validation is then performed on these models

### 14.1_SoilversusSoilGrids_richness.R &  14.1_submitScript_SoilversusSoilGrids_richness.sh
The species richness SoilGrids model and cross validation is performed on the HPC.

### 15_CreateFunctionalDataset.R
Using the species-level dataset, the harmonised species names and functional groups are appended. This is a messy script! But also uses a harmonised names dataframe, which have not been made publically available (the data that has been made open-access, has the harmonised names in).

### 16_FunctionalGroupsDataset.R
After some initial cleaning, tor each site the richness, abundance and biomass of each of the functional groups present is caluclated.

### 17_LatitudinalDiversityGradient.R
Based on the species names, calculating the number of species within each equal sized band across the globe

### 17.1_LDG_fixedNumber.R
Based on the species names, caluclating the number of species within each zone of the globe, when a zone contains the same number of sites.

### 18_FunctionalGroupAnalysis.R
Models based on the functional groups at each site.

### 19_FunctionalGroupFigures.R
Figures for the models based on the fuctional groups.

### 20_Authorships.R
Determining which data we used in this analysis, in order to contact the relevant people for authorship. 

### 21_ReEmailingAuthors.R
We made initital contact (the 6 'sent' dataframes), then wanted to re-email those who had not responded.

### 22_NematodePipelineData.R
Getting the data together to send to JvdH, DR and TGC, so they can create the uncertainty maps.

### 23_OpenData.R
Script that cleans the data ready to make it open access.

### AuthorsandInstitutes.R 
Creates a csv of all the co-authors and their insitutes, ready for the manuscript

### AuthorsandInstitutes_part2.R 
I manually correct the names and address of the authors, then this script assigns the relevant number to them, so the output can be copied into the word document (and final format done)

### DownloadSoilGrids.R
Script to automatically download the SoilGrids data.

### FundingAndAcknowledgements.R
For all the funding and acknowledgements that were given by the co-authors, this script formats the information and create a text string that can be copied into teh word document.

### Meta-data.R
This script calcualtes all the numbers in the manuscript (i.e., how many sites, how many countries)

=======
# sWorm
All work relating to sWorm project

Contains lots of code.
>>>>>>> correction
