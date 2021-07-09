# Created 18-02-20
# Script to process and visualise phylogenetic trees
# Please run this FIRST, as it is a PREREQUISITE for the other scripts to work

# Load packages
library(ape)
library(phytools)
library(scales)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaMay2021Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_04-03-21"
pathNI <- "C:/Users/UCD/Documents/Lab/CVRL MAP/NIMetaOct2020.csv"

# Read in table of bryant isolates
isoBryantTable <- read.table(pathBryantIso,
                             header = TRUE,
                             sep = ",",
                             stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                             check.names=FALSE) # Names left as they are, no dots inserted

# Read in table of CVRL isolates
isoCVRLTable <- read.table(pathNewIso,
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)

# Read in table of NI isolates
isoNITable <- read.table(pathNI,
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

shortCounties <- c("AnNI", "ArNI", "CW", "CN", "CE", "C", "DL", "DoNI", "D", "FNI", "G", "KY", "KE", "KK", "LS", "LM",
                   "L", "DeNI", "LD", "LH", "MO", "MH", "MN", "OY", "RN", "SO", "T", "TNI", "W", "WH", "WX", "WW")

#### Tree file processing ####

# Read in tree
TheTree <- read.tree(pathTree)

# Re-root
reRoot <- root(TheTree, node = 297)
#rotatedTree <- rotateNodes(reRoot, c(409))

# Get the Bryant names
realNames <- getBryantLabels(isoBryantTable, reRoot)

# Update the names in the tree
reRoot$tip.label <- realNames

# Get metadata for CVRL isolates
realNames <- getCVRLLabels(isoCVRLTable, reRoot)

# Update the names in the tree
reRoot$tip.label <- realNames

# Get metadata for NI isolates
realNames <- getNILabels(isoNITable, reRoot)

# Update names
reRoot$tip.label <- realNames

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(reRoot$tip.label)
irishOnlytree <- drop.tip(reRoot, dropper)

# Get rid of non-relevant tips
dropem <- c(17,19,20,95,96,135,136,183,184)
onlytree <- drop.tip(irishOnlytree, dropem)

# Convert branch lengths to SNP values
onlytree$edge.length <- round(onlytree$edge.length * 7843)
reRoot$edge.length <- round(reRoot$edge.length * 7843)

# Find the distances between all isolates
allDist <- cophenetic(onlytree)
euDist <- cophenetic(reRoot)

# Round the distances
allDist <- round(allDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(allDist)){
  
  allDist[index, index] <- NA
}

# Get the herd names
herdNames <- getNames(allDist, "Herd")
euHerds <- getNames(euDist, "Herd")

# Get current county, birth county names and sameness
countyNames <- getNames(allDist, "CCounty")
birthcountyNames <- getNames(allDist, "BCounty")
sameness <- getNames(allDist, "Same")

# Get VNTR names
vntrNames <- getNames(allDist, "VNTR")

# Make VNTR colours
vntrTips <- makeVNTRCols(vntrNames)

# Make EU colours
euCols <- makeRegionColours(reRoot$tip.label)

# Simplify the labels
simpleMat <- deconstructLabels(onlytree$tip.label, counties, shortCounties)

# Assign simple labels
onlytree$tip.label <- simpleMat[1,]

#### Tree plotting (.png) ####

# Save plot as .png file (Ireland)
outputFile <- paste("VNTR_Tree_10-05-21.png", sep="")
png(outputFile, height=10000, width=6000)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 11, font = 1, label.offset = 0.2, 
           tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 5, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(317,285,267,259,249,222,213,199), frame = "n", cex=15, adj = c(1,1), col = "red")


# Add the SNP scale
add.scale.bar(x=10, y = 5, cex = 8, lwd = 15)
text(x=40, y =5, cex = 8, "SNPs")

# Add a legend
legend(x=9, y=160, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types", title.col = "black")

dev.off()

# Make EU plot
outputFile <- paste("EU-Tree_10-05-21.png", sep="")
png(outputFile, height=7000, width=5500)

# Plot EU tree
plot.phylo(reRoot, edge.width = 10, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(472,425,403,391,379,334,303,296), frame = "n", cex=15, adj = c(1,0), col = "red")


#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 9)

# Add the SNP scale
add.scale.bar(x=100, y=0, cex = 10, lwd = 15)
text(x=167, y=0, cex = 10, "SNPs")

# Add a legend
legend(x=130, y=170, legend = c("Ireland", "Great Britain", "Continental Europe"), 
       text.col = c("darkgreen", "firebrick4", "slateblue"), pch = 18,
       col = c("darkgreen", "firebrick4", "slateblue"),
       bty = "n", cex = 12, y.intersp = 0.8, title = "Region", title.col = "black")


dev.off()


#### Tree plotting (.pdf) ####

# Save plot as .pdf file (Ireland)
outputFile <- paste("VNTR_Tree_12-05-21.pdf", sep="")
pdf(outputFile, height=100, width=100)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 17, font = 2, label.offset = 2, tip.color = vntrTips,
           align.tip.label = TRUE, type="fan", cex = 8.5, no.margin = TRUE,x.lim = c(-175,175), y.lim = c(-175,175))
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(317,285,267,259,249,222,213,199), frame = "c", cex=10, col = "darkred")

tiplabels(pch = 18, frame = "n", col = simpleMat[2,], cex=15)
# Add the SNP scale
add.scale.bar(x=-175, y = -175, cex = 10, lwd = 17)
text(x=-100, y =-175, cex = 10, "SNPs")


# Add a legend
legend(x=-175, y=175, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types", title.col = "black")
legend(x=135, y=175, legend = c("2013", "2014", "2015", "2016", "2017", "2019"), 
       text.col = c("firebrick3", "goldenrod3", "gray40", "steelblue3", "palegreen3", "orchid"),
       pch = c(23,23,23,23,23,23), pt.bg = c("firebrick3", "goldenrod3", "gray40", "steelblue3", "palegreen3", "orchid"), pt.cex = 15,
       bty = "n", cex = 10, y.intersp = 0.8, title = "Isolation Year", title.col = "black")

dev.off()

# Make EU plot pdf
outputFile <- paste("EU-Tree_09-07-21.pdf", sep="")
pdf(outputFile, height=120, width=120)

# Plot EU tree
plot.phylo(reRoot, edge.width = 20, font = 2, label.offset = 2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="fan", cex = 30, no.margin = TRUE,x.lim = c(-175,175), y.lim = c(-175,175))
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(472,425,403,391,379,334,303,296), frame = "c", cex=13, col = "darkred")


#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 16)

# Add the SNP scale
add.scale.bar(x=-175, y=-175, cex = 15, lwd = 25)
text(x=-95, y=-175, cex = 15, "SNPs")

# Add a legend
legend(x=-175, y=175, legend = c("Ireland", "Great Britain", "Continental Europe"), 
       text.col = c("darkgreen", "firebrick4", "slateblue"), pch = 18,
       col = c("darkgreen", "firebrick4", "slateblue"),
       bty = "n", cex = 16, y.intersp = 0.8, title = "Region", title.col = "black")

dev.off()
#### Functions ####

# Function to get the labels names for bryant isolates
getBryantLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  chopVector <- TheTree$tip.label
  
  # Remove trailing _x.vcf.gz
  chopchopVector <- sapply(strsplit(chopVector, split = "_"), function(x) (x[1]))
  
  # Remove the > symbol
  nameVector <- sapply(strsplit(chopchopVector, split = ">"), function(x) (x[2]))
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"],
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
        nameVector[index] <- newname
      }else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"],
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
        nameVector[index] <- newname
      }else if(nameVector[index] == "Ref-1997") {
        
        nameVector[index] <- "MAP K10"
      }else{
        
        next
      }
      
    }
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to get the labels names fr CVRL isolates
getCVRLLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"AliquotFormat"] == nameVector[index]){
        
        herd <- strsplit(isoTable[row,"Herd Identifier"], split = " ")[[1]][2]
        
        newname <- paste(nameVector[index], "_", isoTable[row, "County"], "_", herd, "_",
                         isoTable[row,"INMV Group"], "_", isoTable[row, "Birth County"], "_", isoTable[row, "Number of moves to herd of sampling"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to get labels for NI isolates
getNILabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"SeqRef"] == nameVector[index]){
        
        county <- strsplit(isoTable[row,"Herd Ref"], split = " ")[[1]][1]
        herd <- strsplit(isoTable[row,"Herd Ref"], split = " ")[[1]][2]
        
        newname <- paste(isoTable[row,"AliquotFormat"], "_", county, "_", herd, "_",
                         isoTable[row,"INMV"], "_", isoTable[row, "Birth County"], "_", isoTable[row, "Herd of Birth"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
  
}

# Function to create vector with international tips to drop
toDropInternationalTips <- function(tiplabel){
  
  # Create vector to store index values of tips to be dropped
  dropVector <- c()
  
  # Loop thru tip labels and drop as required
  for(index in 1:length(tiplabel)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("UK", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Italy", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index) 
    } else if(grepl("Spain", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("France", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Scotland", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("England", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Wales", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Germany", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Netherlands", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Czech", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Greece", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Norway", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("NewZealand", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("USA", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Canada", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Venezuela", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("India", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Argentina", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("ERR0", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("SRR", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    }
  }
  return(dropVector)
}

# Function to pull out matrix names and simplify to chosen
getNames <- function(mat, chosen){
  
  # Create a vector for names
  rowcolNames <- rep(NA, length(colnames(mat)))
  
  if(chosen == "VNTR"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(rownames(mat)[index], split = "_")[[1]][4]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  } else if(chosen == "Herd"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      one <- strsplit(colnames(mat)[index], split = "_")[[1]][2]
      two <- strsplit(colnames(mat)[index], split = "_")[[1]][3]
      
      vRow <- paste(one,"_",two)
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
    
  } else if(chosen == "CCounty"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][2]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
    
  } else if(chosen == "BCounty"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][5]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  } else if(chosen == "Same"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][6]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  }
}

# Function to create tip labels colours based on VNTR
makeVNTRCols <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    if(is.na(colourVec[index]) == TRUE || colourVec[index] == "n/a"){
      
      colourVec[index] <- "grey30"
    } else if(colourVec[index] == "1"){
      
      colourVec[index] <- "red"
    } else if(colourVec[index] == "2"){
      colourVec[index] <- "deepskyblue3"
    } else if(colourVec[index] == "3"){
      colourVec[index] <- "darkorange3"
    } else if(colourVec[index] == "13"){
      colourVec[index] <- "black"
    } else if(colourVec[index] == "116"){
      colourVec[index] <- "darkgreen"
    }
  }
  return(colourVec)
}

# Function to generate colours based on region
makeRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkgreen"
    } else if(grepl("UK", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Italy", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue" 
    } else if(grepl("Spain", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("France", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Scotland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("England", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Wales", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Germany", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Netherlands", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Czech", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Greece", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Norway", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else {
      
      colourVec[index] <- "darkgreen"
    }
  }
  
  return(colourVec)
}

# Function to simplify the labels
deconstructLabels <- function(tiplabel, counties, shortCounties){
  
  # Create output frame
  outmat <- matrix(nrow = 2, ncol = length(tiplabel))
  
  
  # Loop thru the tips and cut them down
  for(index in 1:length(tiplabel)){
    
    # Split up the different parts of the tip label
    one <- strsplit(tiplabel[index], split = "_")[[1]][2]
    two <- strsplit(tiplabel[index], split = "_")[[1]][3]
    notdate <- strsplit(tiplabel[index], split = "_")[[1]][1]
    date <- strsplit(notdate, split = "-")[[1]][1]
    
    # Find which index in counties the tip county is and store shortened version
    short <- shortCounties[which(one == counties)]
    
    # Store herd
    herd <- paste(one,two)
    
    # Store birth location
    birth <- strsplit(tiplabel[index], split = "_")[[1]][5]

    
    # Store sameness
    same <- strsplit(tiplabel[index], split = "_")[[1]][6]
    
    # Check if it's the same county
    if(same == "n/a" || is.na(same) == TRUE || same == "Unknown" || same == "Not available"|| same == "Notavailable"){
      
      yoke <- paste(herd,"*", collapse = NULL)
    
      outmat[1,index] <- yoke
    
    } else if(same == "None" || same == "Same"){
      
      outmat[1,index] <- herd 
    }else {
      
      if(birth == "U.K. Import" || birth == "U.K.Import"){
        shortB <- "UK"
      } else{
      
        shortB <- shortCounties[which(birth == counties)]
      }
      birthstring <- paste("(",shortB,")", sep = "")
      
      thing <- paste(herd, birthstring)
      
      outmat[1,index] <- thing
    }
    
    if(date == "13"){
      
      co = "firebrick3"
    } else if(date == "14"){
      
      co = "goldenrod3"
    } else if(date == "15"){
      
      co = "gray40"
    } else if(date == "16"){
      
      co = "steelblue3"
    } else if(date == "17"){
      
      co = "palegreen3"
    } else if (date == "19"){
      
      co = "orchid"
    }
      
    outmat[2,index] <- co
  }
  
  return(outmat)
}

