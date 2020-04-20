
## WORKING DIRECTORY --------------
setwd("C:/Users/hp39wasi/WORK/sWorm/EarthwormAnalysis")

## VARIABLES -----------------------

library(dplyr)


## LOAD THE MANUALLY CORRECTED ALLAUTHORS CSV

# allAuthors <- read.csv("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\Authorship\\AllAuthors_manualEdit.csv", encoding = "latin1")
allAuthors <- read.csv("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\Authorship\\AllAuthors_manualEdit_correction.csv", encoding = "latin1")

sworm <- read.csv("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\Authorship\\sWormPeople.csv",encoding = "latin1")
justNames <- sworm[!(duplicated(sworm[,1])),]
fAuthor <- justNames[1:(nrow(justNames)-2),1]
lAuthor <- justNames[(nrow(justNames)-1):nrow(justNames),1]


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
      first = authors.list[grep(first.author[1], authors.list$name),]
      authors.list <- authors.list[ -(grep(first.author[1], authors.list$name)),]
      ### several
      for(i in 2:length(first.author)){
        first = rbind(first, authors.list[grep(first.author[i], authors.list$name),])
        print(i)
        toRemove <- grep(first.author[i], authors.list$name)
        print(toRemove)
        authors.list <- authors.list[ -(toRemove),]
        print(nrow(authors.list))
      }
    }
    ## last
    ### just 1
    if(length(last.author)==1){
      last = authors.list[authors.list$name == last.author,]
      authors.list <- authors.list[ -(which(authors.list$name %in% last.author)),]
      ### several
    }else{
      
      last = authors.list[(grep(paste(last.author,collapse="|"), 
                                authors.list$name)),]
      authors.list <- authors.list[ -(grep(paste(last.author,collapse="|"), 
                                           authors.list$name)),]

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


list.author.fun(dataframe = df,first.author = as.character(fAuthor),last.author = as.character(lAuthor) ,file.name = "authors_list_correction.txt")
