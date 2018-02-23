ApplyQuantiles <- function(x) {
  	breaks = c(quantile(x, probs = seq(0, 1, by = 0.25)))
	cut(x, breaks = breaks + seq_along(breaks) * .Machine$double.eps, ## For urban which has a lot of non-unique values 
      labels=c(1,2,3,4),include.lowest=TRUE)
}


HellsDist <- function(site_data_col, raster_filename, inDir, land.use = TRUE, n.reps = 1000, UD = UD){
	ras <- readGDAL(paste(inDir,raster_filename,sep=""))
	ras <- ras$band1
	ras <- ras[!is.na(ras)]
	if(land.use == TRUE){
		if(max(ras) > 1){
			ras <- ras/1000}
	}
	qrt <- quantile(ras)
	ras <- ApplyQuantiles(ras)

	ras[is.na(ras)] <- 1 ## Lowest value is NA

	
	dat <- site_data_col
	if(land.use == TRUE){
		if(max(dat, na.rm = TRUE) > 1){
			dat <- dat/1000 }
	}

	quants <- c(1:4) ##Essentially the labels, but need to be continuous labels
	l <- n.reps
	HD.vals <- c()
	UD<- DiscreteDistribution(quants, c(0.25,0.25,0.25,0.25))


	for(i in 1:l){
		ras.s <- sample(ras, length(dat), replace = FALSE)
		ras.s <- table(ras.s)/sum(table(ras.s)) ## Proportion of data in each quantile
		DD<- DiscreteDistribution(quants, ras.s)	## Distribtion on quantiles
		HD.vals[i] <- HellingerDist(UD, DD) # Comparing actual distribution with uniform distribution
	}

	return(HD.vals)
}


HellsDistData <- function(site_data_col, raster_filename, inDir, land.use = TRUE, UD = UD){
	ras <- readGDAL(paste(inDir,raster_filename,sep=""))
	ras <- ras$band1
	ras <- ras[!is.na(ras)]
	if(land.use == TRUE){
		if(max(ras) > 1){
			ras <- ras/1000}
	}
	qrt <- quantile(ras)
	dat <- site_data_col
	if(land.use == TRUE){
		if(max(dat, na.rm = TRUE) > 1){
			dat <- dat/1000 }
	}

	h <- hist(dat, breaks = qrt, plot = FALSE)
	h$count <- h$count/sum(h$count)
	DD<- DiscreteDistribution(quants, h$count)
	HD<- HellingerDist(UD, DD)
	
	return(HD)
}
