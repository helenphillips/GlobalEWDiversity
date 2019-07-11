## WORKING DIRECTORY --------------
# setwd("C:/restore2/hp39wasi/sWorm/EarthwormAnalysis")
rm(list=ls())
## VARIABLES -----------------------

data_in <- file.path("Home/data/Phillips")


library(Hmisc)
library(dplyr)

## LOAD LIST ---------------

authors <- read.csv("data/Phillips/Co-Authors on the first analysis paper (Responses)_2019-03-31 - Form Responses 3.csv")

## CHECKING ----------------

table(authors$I.consent.to.my.data.being.used.in.this.paper)

authors <- droplevels(authors[authors$I.would.like.to.be.an.author.on.this.paper == "Yes",])

table(authors$Do.you.have.more.institutions..1) ## Can remove the columns for third institution
notneededCols <- c("Timestamp",  "Email.address",                                                                                            
 "I.consent.to.my.data.being.used.in.this.paper",                                                                         
 "I.would.like.to.be.an.author.on.this.paper",                                                                         
"Do.you.have.more.institutions..1",                                                                         
 "Department.Name.2",                                                                              
 "Research.Centre.Name.2",                                                                                                
"University.Name.2",                                                                              
 "Street.name.and.number.2",                                                                                              
"City.2",                                                                                    
 "Postcode.Zipcode..please.put.quotation.marks.around.it.e.g....90210...2",                                                
 "Country.2", 
"Name.of.the.funding.agency",
 "Funding.code.Grant.Number",
"Name.of.anyone.you.need.to.acknowledge..If.multiple.people..please.separate.with.a.comma......",                         
 "Name.of.the.funding.agency.1",
 "Funding.code.Grant.Number.1",
 "Name.of.anyone.you.need.to.acknowledge..If.multiple.people..please.separate.with.a.comma.......1",
 "Please.enter.your.full.name",
 "Please.enter.an.email.address.that.we.can.contact.you.on",
"If.you.feel.comfortable.doing.so..please.indicate.why.you.no.longer.want.us.to.use.the.data")

authors <- authors[-(which(names(authors) %in% notneededCols))]


## SORT OUT NAMES FIRST ---------------------
names(authors)[which(names(authors) == "Please.enter.your.first..given..name")] <- "firstName"
names(authors)[which(names(authors) == "Please.enter.your.surname..family.name.")] <- "lastName"
names(authors)[which(names(authors) == "Please.enter.any.middle.initials..with.no.punctuation..but.spaces.between.each.initials...Leave.blank.if.no.middle.initials")] <- "middleInits"

# remove upper cases
simpleCap <- function(x) {
        s <- strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1,1)), substring(s, 2),
              sep="", collapse=" ")
}



authors$firstName <- tolower(authors$firstName)
authors$firstName <- sapply(authors$firstName, simpleCap)


authors$lastName <- tolower(authors$lastName)
authors$lastName <- sapply(authors$lastName, simpleCap)

## MIDDLE INITIALS ---------------

authors$middleInits <- paste0(authors$middleInits, ".")
authors$middleInits <- ifelse(authors$middleInits == ".", "", authors$middleInits)
authors$middleInits[authors$middleInits == "Victorovna."] <- "V."
authors$middleInits <- gsub(" ", ". ", authors$middleInits)
authors$middleInits <- gsub("\\.\\.", "\\.", authors$middleInits)

# Where people two initials with no space
authors$middleInits <- gsub("(\\S)(\\S)(\\.)", "\\1\\.\\ \\2\\3" , authors$middleInits)

authors$FullName <- paste(authors$firstName, authors$middleInits, authors$lastName)
authors$FullName <- gsub("  ", " ", authors$FullName)

# There's one author who signed up but didn't need to
authors <- authors[-(grep("Gongalsky", authors$lastName)), ]

# And there's someone who signed up twice

authors[duplicated(authors$FullName),'FullName']
grep(authors[duplicated(authors$FullName),'FullName'], authors$FullName)
# Just chekcing

authors <- authors[-(duplicated(authors$FullName)),]

authors <- authors[-(grep("Juan B.", authors$FullName)[1]),]

## Someone got their country wrong

authors$Country.1[authors$Country.1 == "Brazil" & authors$City.1 == "Berlin"] <- "Germany" 



## COUNTRY ------------------
## If we fix country, then look for duplicates within each

levels(authors$Country)[levels(authors$Country) == "België"] <- "Belgium"
levels(authors$Country)[levels(authors$Country) == "Brasil"] <- "Brazil"
levels(authors$Country)[levels(authors$Country) == "España"] <- "Spain"
levels(authors$Country)[levels(authors$Country) == "México"] <- "Mexico"
levels(authors$Country)[levels(authors$Country) == "The Netherkands"] <- "the Netherlands"
levels(authors$Country)[levels(authors$Country) == "United States"] <- "USA"
levels(authors$Country)[levels(authors$Country) == "US"] <- "USA"
levels(authors$Country)[levels(authors$Country) == "U.S.A."] <- "USA"
levels(authors$Country)[levels(authors$Country) == "Uunited States"] <- "USA"
levels(authors$Country)[levels(authors$Country) == "United States of America"] <- "USA"
authors$Country <- as.character(authors$Country)
authors$Country[authors$Country == "the Netherlands"] <- "The Netherlands" 
authors$Country[authors$Country == "Netherlands"] <- "The Netherlands" 

authors$Country <- as.factor(authors$Country)

## RENAME ----------

names(authors)[which(names(authors) == "Department.Name")] <- "Dept"
names(authors)[which(names(authors) == "Research.Centre.Name")] <- "Centr"
names(authors)[which(names(authors) == "University.Name")] <- "Uni"
names(authors)[which(names(authors) == "Street.name.and.number")] <- "Str"
names(authors)[which(names(authors) == "Postcode.Zipcode..please.put.quotation.marks.around.it.e.g....90210..")] <- "Postcode"


## REMOVE WHITE SPACE ---------


whitespace <- c("firstName", "lastName", "middleInits", "Dept", "^Centr", "Uni$", "Str$", "City$", "Postcode$")
for(c in 1:length(whitespace)){
        whichCol <- which(names(authors) == whitespace[c])        
        authors[,c] <- trimws(authors[,c], which = "both")
}


## CITY ---------------

unique(authors$City)[order(unique(authors$City))]
authors$City <- as.character(authors$City)
authors$City[authors$City == "Muencheberg"] <- "Müncheberg"
authors$City[authors$City == "Rio Cuarto"] <- "Río Cuarto"
authors$City[authors$City == "Fort Collins"] <- "Fort Collins, CO"
authors$City[authors$City == "DAVAO CITY"] <- "Davao City"
authors$City[authors$City == "Gainesvill"] <- "Gainesville, FL"
authors$City[authors$City == "Saint Paul"] <- "St Paul, MN"
authors$City[authors$City == "ATHENS"] <- "Athens"

authors$City <-sapply(authors$City, simpleCap)

## ADD COMMAS ---------------------------


addComma <- c("Dept", "^Centr", "Uni$", "Str$", "City$", "Postcode$")
for(c in 1:length(addComma)){
        whichCol <- grep(addComma[c], names(authors), perl = TRUE)
        print(whichCol)
        authors[,whichCol] <- paste0(authors[,whichCol], ", ")
        authors[,whichCol] <- ifelse(authors[,whichCol] == ", ", "", authors[,whichCol])
}

authors$Address1 <- paste0(authors$Dept, authors$Centr, authors$Uni, authors$Str, authors$City, authors$Country)



## POSTCODE -----------------------------

authors$Postcode <- gsub('\"', '', authors$Postcode)

## LOOKING WITHIN EACH COUNTRY ----------------------

## Make a dataframe
Alist <- data.frame(author = NA, dummyInst = 0)


countries <- unique(authors$Country)


for(c in 1:length(countries)) {
        country <- authors[which(authors$Country == countries[c]), ]
        country <- country[order(country$City),]
        if (nrow(country) == 1) {
                print("Only one address in this country, adding")
                ## If there's only one author in the country then we can add no problem
                ## Add a row, which had the author name, and then the same number of zeros
                # as remaining columns
                
                cols <- ncol(Alist) # The first is the author,
                
                Alist[nrow(Alist) + 1, ] <-
                        c(country$FullName, rep(0, length = (cols - 1)))
                Alist[, ncol(Alist) + 1] <-
                        c(rep(0, length = (nrow(Alist) - 1)), 1)
                names(Alist)[ncol(Alist)] <- country$Address1
        } else{
                print(country$Address1)
                anySameCountry <- readline("Are ANY of these the same?(y/n)")
                if (anySameCountry == "n") {
                        ## ADD ALL TO THE DATAFRAME
                        print("Adding all to the dataframe")
                        
                        for(r in 1:nrow(country)){
                                cols <- ncol(Alist) # The first is the author,
                                
                                Alist[nrow(Alist) + 1, ] <-
                                        c(country$FullName[r], rep(0, length = (cols - 1)))
                                Alist[, ncol(Alist) + 1] <-
                                        c(rep(0, length = (nrow(Alist) - 1)), 1)
                                names(Alist)[ncol(Alist)] <- country$Address1[r]
                        }
                        
                } else{
                        allSameCountry <- readline("Are they ALL the same?(y/n)")
                        if (allSameCountry == "y") {
                                whichAddress <- readline("Which address to use?")
                                ## Add them with the specified address
                                print("Adding all with the specified address")
                                
                                whichAddress <- as.numeric(as.character(whichAddress))
                                
                                Alist[, ncol(Alist) + 1] <- 0
                                names(Alist)[ncol(Alist)] <- country$Address1[whichAddress]
                                
                                
                                cols <- ncol(Alist) # The first is the author,
                                
                                for(r in 1:nrow(country)){
  
                                        Alist[nrow(Alist) + 1, ] <-
                                                c(country$FullName[r], rep(0, length = (cols - 2)), 1)
                                        
                                }
                                
                        } else{
                                cities <- unique(country$City)
                                
                                for(cc in 1:length(cities)){
                                        city <- country[which(country$City == cities[cc]),]
                                        if (nrow(city) == 1) {
                                                print("Only one address in this city, adding")
                                                ## If there's only one author in the city then we can add no problem
                                                ## Add a row, which had the author name, and then the same number of zeros
                                                # as remaining columns
                                                
                                                cols <- ncol(Alist) # The first is the author,
                                                
                                                Alist[nrow(Alist) + 1, ] <-
                                                        c(city$FullName, rep(0, length = (cols - 1)))
                                                Alist[, ncol(Alist) + 1] <-
                                                        c(rep(0, length = (nrow(Alist) - 1)), 1)
                                                names(Alist)[ncol(Alist)] <- city$Address1
                                        } else{
                                                print(city$Address1)
                                                anySameCity <- readline("Are ANY of these the same?(y/n)")
                                                if(anySameCity == "n"){
                                                        # Add all to data
                                                        print("Adding all the dataframe")
                                                        
                                                        
                                                        
                                                        
                                                        for(r in 1:nrow(city)){
                                                                cols <- ncol(Alist) # The first is the author,
                                                                
                                                                Alist[nrow(Alist) + 1, ] <-
                                                                        c(city$FullName[r], rep(0, length = (cols - 1)))
                                                                Alist[, ncol(Alist) + 1] <-
                                                                        c(rep(0, length = (nrow(Alist) - 1)), 1)
                                                                names(Alist)[ncol(Alist)] <- city$Address1[r]
                                                        }
                                                }else{
                                                        howMany <- readline("How many address sets are there?")
                                                        howMany <-as.numeric(as.character(howMany))
                                                        ## Add that set to the dataframe
                                                        for (h in 1:howMany) {
                                                                whichSameCity <- readline("Which are the same?(Choose another set)[numbers separated by commas]")
                                                                
                                                                whichSameCity<- strsplit(whichSameCity, split = ",")
                                                                whichSameCity <- unlist(whichSameCity)
                                                                whichSameCity <- as.numeric(as.character(whichSameCity))
                                                                
                                                                set <- city[whichSameCity,]
                                                                
                                                                print(set$Address1)
                                                                whichCity <- readline("Which address to use?")
                                                                whichCity <- as.numeric(as.character(whichCity))
                                                                ## Add them with the specified address
                                                                
                                                               
                                                                print("Adding the set with the specified address")
                                                                
                                                                
                                                                
                                                                Alist[,ncol(Alist) + 1] <- 0 # Add a new column
                                                                names(Alist)[ncol(Alist)] <- set$Address1[whichCity] # The new last column
                                                                
                                                                cols <- ncol(Alist)
                                                                
                                                                for(r in 1:nrow(set)){
                                                                        
                                                                        
                                                                        Alist[nrow(Alist) + 1, ] <-
                                                                                c(set$FullName[r], rep(0, length = (cols - 2)), 1) # The last column
                                                                        
                                                                }
                                                                
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}

write.csv(Alist, file = "data/Phillips/part1.csv", row.names = FALSE)
Alist <- read.csv("data/Phillips/part1.csv")
Blist <- Alist
## SECOND ADDRESS -------
table(authors$Do.you.have.more.institutions.)
second <- authors[which(authors$Do.you.have.more.institutions. == "Yes"),]

## REMOVE COLS
remove <- c("Dept","Centr","Uni","Str",                                                             
            "City","Postcode","Country","Do.you.have.more.institutions.")

second <- second[,which(!(names(second) %in% remove))]

names(second)[which(names(second) == "Department.Name.1")] <- "Dept"
names(second)[which(names(second) == "Research.Centre.Name.1")] <- "Centr"
names(second)[which(names(second) == "University.Name.1")] <- "Uni"
names(second)[which(names(second) == "Street.name.and.number.1")] <- "Str"
names(second)[which(names(second) == "Postcode.Zipcode..please.put.quotation.marks.around.it.e.g....90210...1")] <- "Postcode"
names(second)[which(names(second) == "Country.1")] <- "Country"
names(second)[which(names(second) == "City.1")] <- "City"


## RENAME COUNTRY --------------------
levels(second$Country)[levels(second$Country) == "México"] <- "Mexico"
levels(second$Country)[levels(second$Country) == "Brasil"] <- "Brazil"
levels(second$Country)[levels(second$Country) == "United States"] <- "USA"

## SOMEONE FILLED OUT THE FORM WRONG -------------------

second <- second[which(second$Country != ""),]

## REMOVE WHITE SPACE ---------


whitespace <- c("Dept", "^Centr", "Uni$", "Str$", "City$", "Postcode$")
for(c in 1:length(whitespace)){
        whichCol <- which(names(second) == whitespace[c])        
        second[,c] <- trimws(second[,c], which = "both")
}

## ADD COMMAS ---------------------------


addComma <- c("Dept", "^Centr", "Uni$", "Str$", "City$", "Postcode$")
for(c in 1:length(addComma)){
        whichCol <- grep(addComma[c], names(second), perl = TRUE)
        print(whichCol)
        second[,whichCol] <- paste0(second[,whichCol], ", ")
        second[,whichCol] <- ifelse(second[,whichCol] == ", ", "", second[,whichCol])
}

second$Address2 <- paste0(second$Dept, second$Centr, second$Uni, second$Str, second$City, second$Country)


## LOOP THROUGH EACH INDIVIDUALLY -------

for(r in 1:nrow(second)){
        print(r)
        print(second$Address2[r])
        country <- second$Country[r]
        country <- as.character(country)
        authors$Country <- as.character(authors$Country)
        temp <- authors[authors$Country == country,]
        
        cat("----------\n")
        print(temp$Address1)
        inOne <- readline("Is the first address in any of these?(y/n)")
        
        if(inOne == "y"){
                # Find author and instiute in original dataframe and add 1
                
                whichOne <- readline("Which one is it the same as?")
                whichOne <- as.numeric(as.character(whichOne))
                
                inst <- temp$Address1[whichOne]
                
                print("Appending to original dataframe")
                
                which(Alist[1,] == inst)
        }else{
                # Add to the dataframe
                print("adding to the dataframe")
                Alist[,ncol(Alist) + 1] <- 0
                names(Alist)[ncol(Alist)] <- second$Address2[r]
                
                AuthorRow <- grep(second$FullName[r], Alist[,1])
                
                Alist[AuthorRow,ncol(Alist)] <- 1
                
        }
        
}

Clist = Alist

## BRING IN SWORM AUTHORS ----------------------

sworm <- read.csv("data/Phillips/sWormPeople.csv",encoding = "latin1")

uniqueSworm <- unique(sworm$Institute)

for(s in 1:length(uniqueSworm)){
        # print(as.character(uniqueSworm[s]))
        
        country <- sub('.*, ', '', uniqueSworm[s])
         print(country)
        
        if(country %in% c("81280-330", "83411-000")){
                country <- "Brazil"
        }
        
        whichCol <- grep(country, names(Alist))
        temp <- names(Alist)[whichCol]
        
        print(as.character(uniqueSworm[s]))
        cat("-----------\n")
        print(temp)
        whichName <- readline("Is the first address in the second?(y/n")
        
        if(whichName == "y"){
                whichRow <- readline("Which row should be replaced by the first?")
                whichRow <- as.numeric(as.character(whichRow))
                names(Alist)[whichRow] <- as.character(uniqueSworm[s])
        }else{
                print("Then we can ignore for now")
        }
        
        
        
}

Dlist = Alist

justNames <- sworm[!(duplicated(sworm[,1])),]
fAuthor <- justNames[1:(nrow(justNames)-2),1]
lAuthor <- justNames[(nrow(justNames)-1):nrow(justNames),1]

## WIDE TO LONG ---------------

# Remove the dummy column and row

Alist[, which(names(Alist) == "dummyInst")] <- NULL
Alist <- Alist[-(which(is.na(Alist$author))),] 


library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
Alist2 <- melt(Alist, id.vars=c("author"))

Alist2 <- Alist2[Alist2$value == 1,]

## Add in first and last names
areNames <- c("firstName", "lastName", "FullName")
authors2 <- authors[,which(names(authors) %in% areNames)]

Alist2 <- merge(Alist2, authors2, by.x = "author", by.y = "FullName")

Alist2$value <- NULL

Alist2 <- Alist2[,c(1, 3, 4, 2)]
names(Alist2) <- c("FullName","firstName", "lastName",  "Address")
names(sworm) <- c("FullName","firstName", "lastName",  "Address")


allAuthors <- rbind(sworm, Alist2)
## 
## USING REMYS SCRIPT -------------------
# THE FUNCTION NEED A DATAFRAME WITH THE FOLLOWING COLUMNS:
#   "Address", "FullName", "lastName", "firstName"
#
# When an author have several institutions, they should be store in several rows
#
# first.author and last.author arguments should be a vector ordered as needed with the FullName of the authors as described in the dataset
list.author.fun <- function(dataframe,first.author = character(), last.author = character(), file.name = "authors.txt") {
        
        ### Prepare two datasets to work on
        address.list = data.frame(address = dataframe$Address, position = rep(NA, length(dataframe$Address)), name =dataframe$FullName)
        authors.list = data.frame(name = dataframe$FullName, first = dataframe$firstName,authors = dataframe$lastName, institute = rep(NA, length(dataframe$lastName)),  address = dataframe$Address)
        
        ### Order of the authors
        options(warn=-1) #remove the warnings
        
        #No first and last authors (alphabetical order _ last name > first name > address)
        if(length(first.author)==0&length(last.author)==0){
                authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
                # No first author
        }else if(length(first.author)==0){
                # last
                ## just 1 
                if(length(last.author)==1){
                        last = authors.list[authors.list$name == last.author,]
                        authors.list <- authors.list[ - (which(authors.list$name%in%last.author)),]
                        
                        ## severals
                }else{
                        last = authors.list[authors.list$name == last.author[1],]
                        authors.list <- authors.list[ -(which(authors.list$name %in% last.author[1])),]
                        for(i in 2:length(last.author)){
                                last = rbind(last, authors.list[authors.list$name == last.author[i],])
                                authors.list <- authors.list[ - (which(authors.list$name%in%last.author[i])),]
                        }
                }
                # order last part (ABC)
                authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
                #paste all 
                authors.list <- rbind(authors.list,last)
                
                # No last author
        }else if(length(last.author)==0){
                # remove first part
                ## if only one
                if(length(first.author)==1){
                        first = authors.list[authors.list$name == first.author,]
                        authors.list <- authors.list[ -(which(authors.list$name%in%first.author)),]
                        # if severals
                }else{
                        first = authors.list[authors.list$name == first.author[1],]
                        for(i in 2:length(first.author)){
                                first = rbind(first, authors.list[authors.list$name == first.author[i],])
                                authors.list <- authors.list[ - (which(authors.list$name%in%first.author[i])),]
                        }
                }
                # order the last part (ABC)
                authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
                #paste all 
                authors.list <- rbind(first,authors.list)
                
                # First and last authors 
        }else{
                # remove first and last
                ## first
                ### just 1
                if(length(first.author)==1){
                        first = authors.list[authors.list$name == first.author,]
                        authors.list <- authors.list[ -(which(authors.list$name%in%first.author)),]
                        ### several
                }else{
                        first = authors.list[authors.list$name == first.author[1],]
                        for(i in 2:length(first.author)){
                                first = rbind(first, authors.list[authors.list$name == first.author[i],])
                                authors.list <- authors.list[ -(which(authors.list$name %in% first.author[i])),]
                        }
                }
                ## last
                ### just 1
                if(length(last.author)==1){
                        last = authors.list[authors.list$name == last.author,]
                        authors.list <- authors.list[ -(which(authors.list$name %in% last.author)),]
                        ### several
                }else{
                        last = authors.list[authors.list$name == last.author[1],]
                        authors.list <- authors.list[ -(which(authors.list$name %in% last.author[1])),]
                        for(i in 2:length(last.author)){
                                last = rbind(last, authors.list[authors.list$name == last.author[i],])
                                authors.list <- authors.list[ -(which(authors.list$name %in% last.author[i])),]
                        }
                }
                # order the middle
                authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
                # paste all 
                authors.list <- rbind(first,authors.list,last)
        }
        options(warn=0)
        #### Ranking the institutions ####
        i = 1
        j = 1
        while (sum(is.na(address.list$position))>0) {
                inst = authors.list$address[i]
                if(is.na(address.list$position[address.list$address==inst])){
                        address.list$position[address.list$address==inst] = j
                        j = j + 1
                        i = i + 1
                }else{
                        i = i + 1
                }
        }
        
        #### Give institutions to authors  ####
        # remove duplicated names
        authors.list = authors.list[!duplicated(authors.list$name),]
        
        # add institution number
        for( i in authors.list$name){
                inst = address.list$position[address.list$name == i]
                inst = inst[order(inst)]
                authors.list$institute[authors.list$name==i] = paste0("(",paste(inst, collapse = ','),")")
        }
        
        address.list <- address.list[order(address.list$position),]
        address.list <- address.list[!duplicated(address.list$position),]
        
        ### Writing text output ###
        addresses <- paste0("^",address.list$position," ",address.list$address)
        names <- paste(paste0(authors.list$name, " ^", authors.list$institute), collapse = ", ")
        
        fileConn<-file(file.name)
        writeLines(c(names," ",addresses), fileConn,useBytes=T)
        close(fileConn)
}

#### try ####
# df  = data.frame(Address = c('lieu1','lieu1','lieu3','lieu2','lieu2','lieu4','lieu5'), 
#                  FullName = c('aA','dB','eB','eB','tF','tD','xQ'), 
#                  lastName = c('A','B','B','B',"F",'D','Q'), 
#                  firstName =  c('a','d','e','e','t','t','x'))
# f = c("xQ","eB")
# l =  c('aA','tD')


df = allAuthors %>%
        apply(., 2, as.character) %>%
        as.data.frame()

df <- df[-(grep("Robin", df$firstName)),]

list.author.fun(dataframe = df,first.author = as.character(fAuthor),last.author = as.character(lAuthor) ,file.name = "authors_list.txt")
