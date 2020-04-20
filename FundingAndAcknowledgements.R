# I manual altered the list
# Remove one duplicate name
# And changing capitalised names to lower case
# Removing punctuation in the middle initials
# Adding a space in the middle initials

list <- read.csv("C:\\Google Drive\\sWorm\\sWorm_WS\\EarthwormAnalysis\\Manuscripts\\SpatialAnalysis\\Authorlist\\Co-Authors on the first analysis paper_2019-03-31_editted.csv")

## List of acknowledgements
acknowle <- list$Name.of.anyone.you.need.to.acknowledge..If.multiple.people..please.separate.with.a.comma......

#remove empty values
acknowle <- acknowle[which(acknowle != "")]
acknowle <- as.character(acknowle)
acknowle <- paste0(acknowle, ",")
acknowle <- do.call(paste, as.list(acknowle, sep = ","))
## Can copy this into Word now


### FUNDING
which(names(list) == "Name.of.the.funding.agency")
which(names(list) == "Funding.code.Grant.Number")

funding <- list[,31:32]
funding <- funding[which(funding$Name.of.the.funding.agency != "" | funding$Funding.code.Grant.Number != ""),]

funding <- paste0(funding$Name.of.the.funding.agency, " (", funding$Funding.code.Grant.Number, ")")

funding <- gsub("\\(\\)", "",  funding)
funding <- trimws(funding, which = "right")
funding <- paste0(funding, ",")
funding <- do.call(paste, as.list(funding))
