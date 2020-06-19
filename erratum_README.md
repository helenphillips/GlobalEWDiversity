# Erratum to 'Global distribution of earthworm diversity', *Phillips et al., 2019 Science*

## Explanation
Whilst preparing another manuscript that used the earthworm data, a bug was found in the data-preparation code. Very specifically, the bug was on `line 129` (but repeated on line `138` and `147`) of the original `0_GetData.R` file. Unfortunately the default value for a function was not changed, and thus specific sites were lost from specific (but not all) datasets. These sites were sites with no earthworm diversity.

Since identifying the error, the bug has been corrected (as well as any code that no longer worked as a result of the bug-fix), and the modelling approached altered (to account for the additional zeros). All underlying data has been extensively checked, and the new data has been made open access with the erratum.

File `0.1_DataCheck.R` were used to check the underlying data once the bug had been fixed. File `erratumFigures.R` was used to create the figure for the erratum, as well as produce the numbers of datasets affected, etc.

All code available on the GitHub repository have fixed the bug. Older code is still available as git commits.

We apologise for this error and any confusion that is has caused. 


## Erratum text
Since publication of the report “Global distribution of earthworm diversity” (Phillips et al., 2019), we have discovered a bug in our data-preparation code. This bug only affected some datasets, but removed from certain datasets all sites where no earthworms were recorded (“zero-sites”). Unaffected datasets still retained zeros. This bug affected all analyses in the report, generally causing our previous models to have overestimated earthworm richness, abundance and biomass, which was reflected in the predicted overall means decreasing in the revised analyses. Following re-analysis, there were also changes in the inferred geographic patterns, but only minor effects on significant environmental drivers.

In total, 40% of studies (51 datasets, 64 studies) in the richness model, 33% of studies (56 datasets, 69 studies) in the abundance model, and  28% of studies (30 datasets, 36 studies) in the biomass model were affected (Figure 1). Of those affected, most datasets only had a relatively small number of zero-sites. However, 13 affected studies used in the richness model, 14 affected studies used in the abundance model, and 7 affected studies in the biomass model more than doubled in size when the zero-sites were re-added to the underlying datasets. Four datasets were also removed during the re-analysis, as upon detailed inspection they did not meet the criteria for inclusion in relation to number of sites (see Supplementary Materials of Phillips et al., 2019).

In this Erratum, all figures, text, and supplementary materials have been updated to reflect the corrected data and analyses. The modelling approach was also changed to account for the zero-inflation in the data. The results of the re-analysis are qualitatively similar to the original paper. The general biodiversity patterns in the three maps have stayed relatively consistent, though two regions -- North America and Central Asia – deviate considerably from the previous patterns of species richness and abundance. As North American studies previously had a large number of zero-sites removed, the model now fits better in the northern USA/Canada region. The predictions in certain parts of Europe and Asia are likely affected by the North American data, and thus, the predictions of both regions are now less extrapolated across environmental space. The patterns in the global map of earthworm biomass are consistent in both versions.

Incorporation of the zero-sites has reduced our global estimates considerably. The mean global estimates are now, for richness 0.95 species (SD = 0.96; previously 2.42 species; SD = 2.19), for abundance 20.59 individuals per m2 (SD = 24.84; previously 77.89 individuals per m2; SD = 98.94), and for biomass 51.18 g per m2 (SD = 4881.2; median, 6.16; previously 2772.8 g per m2; SD = 1,312,782; median, 6.69). Climate and habitat cover are still more important drivers of earthworm communities than soil properties, though there are some changes to the importance of the different themes. Most noticeably, habitat cover has increased in importance in both the abundance and biomass model. The importance of temperature has also changed in the abundance and biomass models (decreasing and increasing, respectively). The text and figures now reflect all these changes.

The error did not quantitively affect any of the overall conclusions drawn in the original paper. We apologize for the mistake and provide a full description of the re-analysis and all changes to the text on GitHub (https://github.com/helenphillips/GlobalEWDiversity). The updated data is now available on the iDiv Data portal (https://doi.org/10.25829/idiv.1818-13-3001). 
This erratum also adds an additional author to the author list.
