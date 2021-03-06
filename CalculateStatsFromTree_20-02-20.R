# Created on 20-02-20
# Script to calculate a few statistics and numbers from the phylogenetic tree
# This script is NOT independent, you are REQUIRED to run DoPhylogeneticTrees_18-02-20.R before running this one

# Number of herds
length(table(herdNames))

# Isolate VNTR types
table(vntrNames)

# Herds with more than 1 isolate
sum(table(herdNames)>1)

# Find how many herds have multiple VNTRs
multivntr <- findHerdsWithMultiVNTR(herdNames, vntrNames)
table(multivntr)

# Max SNP distance
max(allDist, na.rm = T)

# Min/Max distances for INMV 1, 2 & 3
range(allDist[which(vntrNames == 1),which(vntrNames == 1)], na.rm=T)
range(allDist[which(vntrNames == 2),which(vntrNames == 2)], na.rm=T)
range(allDist[which(vntrNames == 3),which(vntrNames == 3)], na.rm=T)
range(allDist[which(vntrNames == 13),which(vntrNames == 13)], na.rm=T)

# Median SNP distances for INMV 1 & 2
median(allDist[which(vntrNames == 1),which(vntrNames == 1)], na.rm=T)
median(allDist[which(vntrNames == 2),which(vntrNames == 2)], na.rm=T)
median(allDist[which(vntrNames == 13),which(vntrNames == 13)], na.rm=T)

# INMV 116 proximity to other isolates
range(allDist[which(vntrNames == 116),], na.rm=T)

# INMV 3 isolates closest neighbour
min(allDist[which(vntrNames == 3)[1],], na.rm=T)
min(allDist[which(vntrNames == 3)[2],], na.rm=T)

# Look at Group A on Irish tree
irishA <- extract.clade(onlytree, node = 317)
length(irishA$tip.label)
distA <- cophenetic.phylo(irishA)
countiesA <- findSimpleCounties(rownames(distA))
length(table(countiesA))

# Look at Group H on Irish tree
irishH <- extract.clade(onlytree, node = 285)
length(irishH$tip.label)
distH <- cophenetic.phylo(irishH)
countiesH <- findSimpleCounties(rownames(distH))
length(table(countiesH))

# Amount of Cork isolates & SNP distances
length(grep("Cork", herdNames))
range(allDist[grep("Cork", herdNames), grep("Cork", herdNames)], na.rm=T)
median(allDist[grep("Cork", herdNames), grep("Cork", herdNames)], na.rm=T)

# Same for Clare 
length(grep("Clare", herdNames))
range(allDist[grep("Clare", herdNames), grep("Clare", herdNames)], na.rm=T)
median(allDist[grep("Clare", herdNames), grep("Clare", herdNames)], na.rm=T)

# Westmeath 1 stats
range(allDist[grep("Westmeath _ 1", herdNames)[4],grep("Westmeath _ 1", herdNames)],na.rm=T)
range(allDist[grep("Westmeath _ 1", herdNames)[5],grep("Westmeath _ 1", herdNames)],na.rm=T)

# Max dist within group A subgroup of 21
max(allDist[115:146, 115:146], na.rm = T)
table(herdNames[115:146])
table(countyNames[115:146])

#### European groups ####

# Group A europe
euroA <- extract.clade(reRoot, node = 472)
euroAmat <- cophenetic(euroA)

for(index in 1:nrow(euroAmat)){
  
  euroAmat[index, index] <- NA
}

euroAmat[upper.tri(euroAmat)] <- NA

# Group B europe
euroB <- extract.clade(reRoot, node = 425)
euroBmat <- cophenetic(euroB)

for(index in 1:nrow(euroBmat)){
  
  euroBmat[index, index] <- NA
}

euroBmat[upper.tri(euroBmat)] <- NA

# Group C
euroC <- extract.clade(reRoot, node = 403)
euroCmat <- cophenetic(euroC)

for(index in 1:nrow(euroCmat)){
  
  euroCmat[index, index] <- NA
}

euroCmat[upper.tri(euroCmat)] <- NA

# Group D
euroD <- extract.clade(reRoot, node = 391)
euroDmat <- cophenetic(euroD)

for(index in 1:nrow(euroDmat)){
  
  euroDmat[index, index] <- NA
}

euroDmat[upper.tri(euroDmat)] <- NA

# Group E
euroE <- extract.clade(reRoot, node = 379)
euroEmat <- cophenetic(euroE)

for(index in 1:nrow(euroEmat)){
  
  euroEmat[index, index] <- NA
}

euroEmat[upper.tri(euroEmat)] <- NA

# Group F
euroF <- extract.clade(reRoot, node = 334)
euroFmat <- cophenetic(euroF)

for(index in 1:nrow(euroFmat)){
  
  euroFmat[index, index] <- NA
}

euroFmat[upper.tri(euroFmat)] <- NA

# Group G
euroG <- extract.clade(reRoot, node = 303)
euroGmat <- cophenetic(euroG)

for(index in 1:nrow(euroGmat)){
  
  euroGmat[index, index] <- NA
}

euroGmat[upper.tri(euroGmat)] <- NA

# Group H
euroH <- extract.clade(reRoot, node = 296)
euroHmat <- cophenetic(euroH)

for(index in 1:nrow(euroHmat)){
  
  euroHmat[index, index] <- NA
}

euroHmat[upper.tri(euroHmat)] <- NA


# Run function to find similar isolates for group A
similarsA <- findSimilars(euroAmat, 30, counties)

# Group B
similarsB <- findSimilars(euroBmat, 30, counties)

# Group C
similarsC <- findSimilars(euroCmat, 30, counties)

# Group D
similarsD <- findSimilars(euroDmat, 30, counties)

# Group E
similarsE <- findSimilars(euroEmat, 30, counties)

# Run function to find similars for group F
similarsF <- findSimilars(euroFmat, 30, counties)

# Group G
similarsG <- findSimilars(euroGmat, 30, counties)

# Group H
similarsH <- findSimilars(euroHmat, 30, counties)

#### Functions ####

# Function to loop thru the matrix and find locations within the SNP threshold
findSimilars <- function(mat, threshold, counties){
  
  # Create result vector
  vector <- c()
  
  # Fetch the names for the matrix
  matNames <- rownames(mat)
  
  # Loop thru each row 
  for(row in 1:nrow(mat)){
    
    # Loop thru each col
    for(col in 1:ncol(mat)){
      
      # Skip yokey if it's NA
      if(is.na(mat[row,col])){
        
        next
      } else if(mat[row,col] > threshold){
        
        next
      } else{
        
        # Split to check location
        loc1 <- strsplit(matNames[row], split = "_")[[1]][2]
        loc2 <- strsplit(matNames[col], split = "_")[[1]][2]
        
        loc3 <- strsplit(matNames[row], split = "_")[[1]][1]
        loc4 <- strsplit(matNames[col], split = "_")[[1]][1]
        
        if(loc1 %in% counties ==  TRUE && loc2 %in% counties == TRUE){
          
          next
        } else if(loc3 == loc4){
          
          
          next
        }else if(loc1 %in% counties ==  TRUE && grepl("Ireland", loc4) == TRUE){
          
          
          next
        }else if(loc2 %in% counties ==  TRUE && grepl("Ireland", loc3) == TRUE){
          
          
          next
        } else { 
          # Paste together the row and col names and SNP count
          current <- paste(matNames[row], "~", matNames[col], "~", mat[row,col])
          
          vector <- append(vector, current)
        }
      }
    }
  }
  return(vector)
}

# Function to find how many herds have different VNTR types present
findHerdsWithMultiVNTR <- function(herds, VNTRs){
  
  herdstatus <- c()
  
  herdstocheck <- which(table(herds)>1)
  
  # Loop thru herds with more than 1 isolate
  for(index in names(herdstocheck)){
    
    places <- which(herds == index)
    
    currentvntrs <- rep(NA, length(places))
    
    for(place in 1:length(places)){
      
      currentvntrs[place] <- VNTRs[places[place]]
    }
    
    herdstatus <- append(herdstatus, length(unique(currentvntrs)))
    
  }
  
  return(herdstatus)
}

# Function to get names of counties from simple lbaels
findSimpleCounties <- function(inputvec){
  
  result <- rep(NA, length(inputvec))
  
  for(index in 1:length(inputvec)){
    
    result[index] <- strsplit(inputvec[index], " ")[[1]][2]
  }
  
  return(result)
}
