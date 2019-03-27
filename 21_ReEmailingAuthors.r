setwd("C:/restore2/hp39wasi/sWorm/EarthwormAnalysis")

sent <- read.csv("C:\\Google Drive\\sWorm\\sWorm_WS\\EarthwormAnalysis\\Manuscripts\\SpatialAnalysis\\DataProviders\\members_DataProvider_Take1_sent_Mar_26_2019.csv")
sent2 <- read.csv("C:\\Google Drive\\sWorm\\sWorm_WS\\EarthwormAnalysis\\Manuscripts\\SpatialAnalysis\\DataProviders\\members_DataProvider_Take2_sent_Mar_26_2019.csv")

sent$source <- NULL

sent <- rbind(sent, sent2)

form <- read.csv("C:\\Google Drive\\sWorm\\sWorm_WS\\EarthwormAnalysis\\Manuscripts\\SpatialAnalysis\\DataProviders\\Co-Authors_Responses_2019-03-26.csv")

form_emails <- form$Email.address
form_surnames <- form$Please.enter.your.surname..family.name.

all(form_emails %in% sent$Email.Address)
which(form_emails %in% sent$Email.Address)

all(form_surnames %in% sent$Last.Name)
which(!(form_surnames %in% sent$Last.Name))

## So remove from 'sent' rows where email address are already in the 'form'

tosend <- sent[-(which(sent$Email.Address %in% form_emails)),]


## Some people used different email addresses to sign up
## So a list of surnames

tosend[which(tosend$Last.Name %in% form_surnames),]

tosend <- tosend[-(which(tosend$Last.Name %in% form_surnames)),]


## To this file, I will do a manual check (becuase of formatting of anmes)
## and append on previously bounced emails
## Then resend to the list

write.csv(tosend, file = file.path("Authorship", "toResend_2019-03-26.csv"), row.names = FALSE)
