## A function to remove single, separated outliers
##
## Input: trait: a vector with measured data
##
## Output: a vector with cleaned data
outlier <- function(trait) {
  
  # calculate number of expected values in 0.5*mad neighbourhood (arbitrary choice)
  mad <- mad(trait, na.rm=TRUE)
  subjIDs <- c()
  IDs <- 1:length(trait)
  nums <- max(c(0.0001*sum(!is.na(trait)), 20))

  count <- 0 # count variable: if for 100 observations in a row sufficient nearby observations are identified: stop
  
  sort.trait <- sort(trait)
  order.trait <- order(trait)
  
  # check lowest and highest 10000 observations for separated outliers
  for (i in 1:10000) { 
    bounds <- c(sort.trait[i]-0.5*mad, sort.trait[i]+0.5*mad)
    neigh <- sum((bounds[1] < trait) & (bounds[2] > trait), na.rm=TRUE)
    count <- count + 1
    if (neigh<=nums) {
      subjIDs <- c(subjIDs, IDs[order.trait[i]])
      count <- 0
    }
    if (count > 100) {
      break
    }
  }
  
  # check lowest and highest 10000 observations for separated outliers
  sort.traitD <- sort(trait, decreasing=TRUE)
  order.traitD <- order(trait, decreasing=TRUE)
  count <- 0
  for (i in 1:10000) { ##?? modify
    bounds <- c(sort.traitD[i]-0.5*mad, sort.traitD[i]+0.5*mad)
    neigh <- sum((bounds[1] < trait) & (bounds[2] > trait), na.rm=TRUE)
    count <- count + 1
    if (neigh<=nums) {
      subjIDs <- c(subjIDs, IDs[order.traitD[i]])
      count <- 0
    }
    if (count > 100) {
      break
    }
  }
  # set outliers to NA
  if (length(trait) > 0) {
    trait[which(IDs %in% subjIDs)] <- NA
  }
  
  return(trait)
  
}

## A function to remove the lowest and highest 0.1% of the samples
##
## Input: trait: a vector with measured data
##
## Output: a vector with cleaned data
perc_delete <- function(trait) {
  low_quan <- quantile(trait, na.rm=T, probs=0.001)
  up_quan <- quantile(trait, na.rm=T, probs=0.999)
  trait[which(trait < low_quan | trait > up_quan)] <- NA
  return(trait)
}


## A function to remove the delete observations with value 0 (does not work for log transformation)
##
## Input: trait: a vector with measured data
##
## Output: a vector with cleaned data
null_delete <- function(trait) { 
  trait[which(trait==0)] <- NA
  return(trait)
}


