#### try ####
# df  = data.frame(Address = c('lieu1','lieu1','lieu3','lieu2','lieu2','lieu4','lieu5'), 
#                  FullName = c('aA','dB','eB','eB','tF','tD','xQ'), 
#                  lastName = c('A','B','B','B',"F",'D','Q'), 
#                  firstName =  c('a','d','e','e','t','t','x'))
# f = c("xQ","eB")
# l =  c('aA','tD')
# list.author.fun(df,first.author = f,last.author = l ,file.name = "try2.txt")


# THE FUNCTION NEED A DATAFRAME WITH THE FOLLOWING COLUMNS:
#   "Address", "FullName", "lastName", "firstName"
# When an athir have several institutions, they should be store in several rows
# first.author and last.author arguments should be a vector ordered as needed of the fullname of the authors as described in the dataset
list.author.fun <- function(dataframe,first.author = character(), last.author = character(), file.name = "authors.txt") {

  address.list = data.frame(address = dataframe$Address, position = rep(NA, length(dataframe$Address)), name =dataframe$FullName)
  authors.list = data.frame(name = dataframe$FullName, first = dataframe$firstName,authors = dataframe$lastName, institute = rep(NA, length(dataframe$lastName)),  address = dataframe$Address)
  
  ### NEED one institution per row, duplicate if two ###
  ## TO COMPLETE LATER
  if(length(first.author)==0&length(last.author)==0){
    authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
  }else if(length(first.author)==0){
    
  }else if(length(last.author)==0){
    
  }else{
    # move first and last\
    # first
    if(length(first.author)==1){
      first = authors.list[authors.list$name == first.author,]
      authors.list <- authors.list[ - which(authors.list$name==first.author),]
    }else{
      first = authors.list[authors.list$name == first.author[1],]
      for(i in 2:length(first.author)){
        first = rbind(first, authors.list[authors.list$name == first.author[i],])
        authors.list <- authors.list[ - which(authors.list$name==first.author[i]),]
      }
    }
    # last
    if(length(last.author)==1){
      last = authors.list[authors.list$name == last.author,]
      authors.list <- authors.list[ - which(authors.list$name==last.author),]
    }else{
      last = authors.list[authors.list$name == last.author[1],]
      for(i in 2:length(last.author)){
        last = rbind(last, authors.list[authors.list$name == last.author[i],])
        authors.list <- authors.list[ - which(authors.list$name==last.author[i]),]
      }
    }
    # order the middle
    authors.list <- authors.list[order(authors.list$authors,authors.list$first,authors.list$address),]
    #paste all 
    authors.list <- rbind(first,authors.list,last)
  }
  
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
    print(address.list$position)
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
  
  address.list = address.list[order(address.list$position),]
  
  
  ### Writing text output ###
  addresses <- paste0("^",address.list$position," ",address.list$address)
  names <- paste(paste0(authors.list$name, " ^", authors.list$institute), collapse = ", ")
  
  fileConn<-file(file.name)
  writeLines(c(names," ",addresses), fileConn)
  close(fileConn)
}