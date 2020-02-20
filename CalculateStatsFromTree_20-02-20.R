# Created on 20-02-20
# Script to calculate a few statistics and numbers from the phylogenetic tree
# This script is NOT independent, you are REQUIRED to run DoPhylogeneticTrees_18-02-20.R before running this one

# Number of herds
length(table(herdNames))

# Isolate VNTR types
table(vntrNames)

# Herds with more than 1 isolate
sum(table(herdNames)>1)

# Max SNP distance
max(allDist, na.rm = T)

# Min/Max distances for INMV 1, 2 & 3
range(allDist[which(vntrNames == 1),which(vntrNames == 1)], na.rm=T)
range(allDist[which(vntrNames == 2),which(vntrNames == 2)], na.rm=T)
range(allDist[which(vntrNames == 3),which(vntrNames == 3)], na.rm=T)

# INMV 116 proximity to other isolates
range(allDist[which(vntrNames == 116),], na.rm=T)

# INMV 3 isolates closest neighbour
min(allDist[which(vntrNames == 3)[1],], na.rm=T)
min(allDist[which(vntrNames == 3)[2],], na.rm=T)

# Amount of Cork isolates & SNP distances
length(grep("Cork", herdNames))
range(allDist[grep("Cork", herdNames), grep("Cork", herdNames)], na.rm=T)
median(allDist[grep("Cork", herdNames), grep("Cork", herdNames)], na.rm=T)

# Same for Clare 
length(grep("Clare", herdNames))
range(allDist[grep("Clare", herdNames), grep("Clare", herdNames)], na.rm=T)
median(allDist[grep("Clare", herdNames), grep("Clare", herdNames)], na.rm=T)

# Westmeath 1 stats
range(allDist[grep("Westmeath", herdNames)[4],grep("Westmeath", herdNames)],na.rm=T)
range(allDist[grep("Westmeath", herdNames)[5],grep("Westmeath", herdNames)],na.rm=T)

# Group A stats
groupA <- extract.clade(onlytree, node = 198)
groupAmat <- cophenetic(groupA)
max(groupAmat)

# Group A europe
euroA <- extract.clade(euOnlyTree, node = 365)
euroAmat <- cophenetic(euroA)

for(index in 1:nrow(euroAmat)){
  
  euroAmat[index, index] <- NA
}

euroAmat[upper.tri(euroAmat)] <- NA

# Groups B
euroB <- extract.clade(euOnlyTree, node = 334)
euroBmat <- cophenetic(euroB)

for(index in 1:nrow(euroBmat)){
  
  euroBmat[index, index] <- NA
}

euroBmat[upper.tri(euroBmat)] <- NA

# Group C
euroC <- extract.clade(euOnlyTree, node = 316)
euroCmat <- cophenetic(euroC)

for(index in 1:nrow(euroCmat)){
  
  euroCmat[index, index] <- NA
}

euroCmat[upper.tri(euroCmat)] <- NA

# Group D
euroD <- extract.clade(euOnlyTree, node = 278)
euroDmat <- cophenetic(euroD)

for(index in 1:nrow(euroDmat)){
  
  euroDmat[index, index] <- NA
}

euroDmat[upper.tri(euroDmat)] <- NA

# Group F
euroF <- extract.clade(euOnlyTree, node = 240)
euroFmat <- cophenetic(euroF)

for(index in 1:nrow(euroFmat)){
  
  euroFmat[index, index] <- NA
}

euroFmat[upper.tri(euroFmat)] <- NA

# Group G
euroG <- extract.clade(euOnlyTree, node = 229)
euroGmat <- cophenetic(euroG)

for(index in 1:nrow(euroGmat)){
  
  euroGmat[index, index] <- NA
}

euroGmat[upper.tri(euroGmat)] <- NA


# Run function to find similar isolates for group A
similarsA <- findSimilars(euroAmat, 30, counties)

# Group B
similarsB <- findSimilars(euroBmat, 30, counties)

# Group C
similarsC <- findSimilars(euroCmat, 30, counties)

# Group D
similarsD <- findSimilars(euroDmat, 30, counties)

# Run function to find similars for group F
similarsF <- findSimilars(euroFmat, 30, counties)

# Group G
similarsG <- findSimilars(euroGmat, 30, counties)


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