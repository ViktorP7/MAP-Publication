# Created on 19-02-20
# Script to generate isolate map of Ireland and temporal (supplementary) plot
# This script is NOT independent, you are REQUIRED to run DoPhylogeneticTrees_18-02-20.R before running this one

# Load packages
library(sp)
library(raster)
library(gplots)
library(maptools)
library(shape)


# Create path strings
pathHerdStats <- "C:/Users/UCD/Documents/Lab/HerdTbStatistics_2010-2019.csv"
niData <- "C:/Users/UCD/Documents/Lab/NICattleData2018.csv"

# Read in herd statistics
statistics <- readTBStatisticsFile(pathHerdStats)

# Calculate the per county proportion test positive animals in most recent report - 2018Q3
summaryTables <- calculateSummaryStatisticsPerQuarter(statistics)

# Get the relevant year
relevantYear <- summaryTables$`2018Q2`

# Update with Northern Ireland info (https://www.daera-ni.gov.uk/articles/agricultural-census-northern-ireland)
northTable <- read.table(niData,
                         header = TRUE,
                         sep = ",",
                         stringsAsFactors=FALSE,
                         check.names=FALSE)

# Merge the two tables
merged2018Cattle <- rbind(relevantYear, northTable)

# Calculate the total of cattle and total of herds
cattleTotal <- sum(merged2018Cattle[1,3], merged2018Cattle[nrow(merged2018Cattle),3])
herdTotal <- sum(merged2018Cattle[1,2], merged2018Cattle[nrow(merged2018Cattle),2])


#### Processing of files for map building ####
# Set path variables
pathCoords <- "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Get the numbers of isolates
sampleNumbers <- numberSamples(herdNames)

# Get polygon coordinates for map plotting
polygonCoords <- getSpatialPolygons(counties, pathCoords)

# Calculate limits of plot
ranges <- mapSpatialLimits(polygonCoords, counties)

# Calculate the max number of samples
maxNSamples <- sampleMax(sampleNumbers)

#### Plot map as .png ####

# Save plot as .png file
outputFile <- paste("IrishMap_14-04-20.png", sep="")
png(outputFile, height=4500, width=4500)

par(bg=NA)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Add legend
legend("topleft", legend = c("1","2","3","4","5","6"), title="Isolates per Herd", bty="n", cex=15,
       pch = c(21,21,21,21,21,21), col = alpha("red", 0.9), pt.bg = alpha("red", 0.65), 
       pt.cex = c(5,10,15,20,25,30), ncol = 2)

# Add a second legend for the proportions of cattle/herds
legend("bottomright", legend = c(0,0.032,0.064,0.096,0.128,0.16), 
       title = "Proportion of Total Cattle Present", bty = "n", cex = 9, 
       pch = rep(22,8), col = "black",
       pt.bg = c(alpha("blue", 0),alpha("blue", 0.2),alpha("blue", 0.4),alpha("blue", 0.6),
                 alpha("blue", 0.8),alpha("blue", 1)),
       ncol = 2, pt.cex = 10)

# Plot county polygons and related sample data
polygonsSpatialData(polygonCoords, sampleNumbers, counties, shortCounties, herdNames, merged2018Cattle, 1140142, 3)

dev.off()

#### Plot map as .pdf ####

# Save plot as .pdf file
outputFile <- paste("IrishMap_14-04-20.pdf", sep="")
pdf(outputFile, height=75, width=75)

par(bg=NA)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Add legend
legend("topleft", legend = c("1","2","3","4","5","6"), title="Isolates per Herd", bty="n", cex=15,
       pch = c(21,21,21,21,21,21), col = alpha("red", 0.9), pt.bg = alpha("red", 0.65), 
       pt.cex = c(5,10,15,20,25,30), ncol = 2)

# Add a second legend for the proportions of cattle/herds
legend("bottomright", legend = c(0,0.032,0.064,0.096,0.128,0.16), 
       title = "Proportion of Total Cattle Present", bty = "n", cex = 10, 
       pch = rep(22,8), col = "black",
       pt.bg = c(alpha("blue", 0),alpha("blue", 0.2),alpha("blue", 0.4),alpha("blue", 0.6),
                 alpha("blue", 0.8),alpha("blue", 1)),
       ncol = 2, pt.cex = 10)

# Plot county polygons and related sample data
polygonsSpatialData(polygonCoords, sampleNumbers, counties, shortCounties, herdNames, merged2018Cattle, 1140142, 3)

dev.off()

#### Temporal Plot ####

# Get the count matrix
counterMatrix <- getYearCounts(herdNames, allDist)

# Save plot as .pdf file
outputFile <- paste("TemporalIsolates_14-04-20.pdf", sep="")
pdf(outputFile, height=20, width=15)

par(cex.main = 2.5, cex.lab = 1.5)
# Plot the figure
heatmap.2(counterMatrix, dendrogram = "none", trace = "none", Rowv = FALSE,
          Colv = FALSE, col = colorRampPalette(c("white", "blue")), denscol = "red",
          xlab = "", ylab = "", main = "Temporal Span of Isolates Across Herds", 
          offsetRow = -83, cexRow = 1.7, density.info = "none", keysize = 0.5, margins = c(7,7),
          key.xlab = "Number of Isolates", cexCol = 1.7, key.title = "Key")
mtext(text="YEAR", side=1, line=3.5, cex = 2)
mtext(text="HERD", side=4, line=0.5, cex = 2)


dev.off()


#### Movement panel plot ####

# Store the tip ranges for the groups
tipsA <- 110:115
tipsB <- 66:109
tipsC <- 60:64
tipsD <- 46:59
tipsE <- 30:45
tipsF <- 117:123
tipsG <- 10:28
tipsH <- 1:8

# Save plot as .pdf file
outputFile <- paste("IrishMoves_20-04-20.pdf", sep="")
pdf(outputFile, height=70, width=80)

par(bg=NA)
par(mfrow=c(2,4))
       
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsA], countyNames[tipsA], "Group A Map", vntrNames[tipsA], sameness[tipsA], herdNames[tipsA])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsB], countyNames[tipsB], "Group B Map", vntrNames[tipsB], sameness[tipsB], herdNames[tipsB])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsC], countyNames[tipsC], "Group C Map", vntrNames[tipsC], sameness[tipsC], herdNames[tipsC])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsD], countyNames[tipsD], "Group D Map", vntrNames[tipsD], sameness[tipsD], herdNames[tipsD])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsE], countyNames[tipsE], "Group E Map", vntrNames[tipsE], sameness[tipsE], herdNames[tipsE])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsF], countyNames[tipsF], "Group F Map", vntrNames[tipsF], sameness[tipsF], herdNames[tipsF])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsG], countyNames[tipsG], "Group G Map", vntrNames[tipsG], sameness[tipsG], herdNames[tipsG])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsH], countyNames[tipsH], "Group H Map", vntrNames[tipsH], sameness[tipsH], herdNames[tipsH])

dev.off()

# Save plot as .png file
outputFile <- paste("IrishMoves_20-04-20.png", sep="")
png(outputFile, height=4500, width=5500)

par(bg=NA)
par(mfrow=c(2,4))

plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsA], countyNames[tipsA], "Group A Map", vntrNames[tipsA], sameness[tipsA], herdNames[tipsA])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsB], countyNames[tipsB], "Group B Map", vntrNames[tipsB], sameness[tipsB], herdNames[tipsB])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsC], countyNames[tipsC], "Group C Map", vntrNames[tipsC], sameness[tipsC], herdNames[tipsC])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsD], countyNames[tipsD], "Group D Map", vntrNames[tipsD], sameness[tipsD], herdNames[tipsD])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsE], countyNames[tipsE], "Group E Map", vntrNames[tipsE], sameness[tipsE], herdNames[tipsE])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsF], countyNames[tipsF], "Group F Map", vntrNames[tipsF], sameness[tipsF], herdNames[tipsF])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsG], countyNames[tipsG], "Group G Map", vntrNames[tipsG], sameness[tipsG], herdNames[tipsG])
plotMoves(ranges, polygonCoords, counties, shortCounties, birthcountyNames[tipsH], countyNames[tipsH], "Group H Map", vntrNames[tipsH], sameness[tipsH], herdNames[tipsH])

dev.off()
  

#### Example herds plot ####

# Get WH1 indices
whgroup = which(herdNames %in% "Westmeath _ 1")
whg1 = whgroup[1:3]
whg2 = whgroup[4]
whg3 = whgroup[5]

# Create example tree
wh1tree <- makeExampleTree(whg1,whg2,whg3,allDist,onlytree)

# Assign colours 
wh1cols <- c(1:length(wh1tree$tip.label)) 
wh1cols[7] <- "gold"
wh1cols[6] <- "darkgreen"
wh1cols[c(1:5)] <- "red"

outputFile <- paste("WH1-0_06-05-20.pdf", sep="")
pdf(outputFile, height=85, width=75)

# Set margins to nothing
currentMar <- par()$mar
par(mar=c(0,0,0,0), fig=c(0,1,0,1))
par(bg=NA)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Plot county polygons and related sample data
makeExampleHerd(whg1,whg2,whg3,counties,polygonCoords,shortCounties,allDist,countyNames, birthcountyNames)


# Add a legend
legend("bottomright", legend = c("A", "B", "G", "A, B, G"), 
       text.col = c("gold", "darkgreen", "red", "blue"), 
       bty = "n", cex = 12, title = "Westmeath 1 Strains", title.col = "black")


# Set figure parameters to top right corner 
par(fig=c(0,0.38,0.75,1), new=T)


# Plot VNTR tree
plot.phylo(wh1tree, edge.width = 12, font = 1, label.offset = 0.2, 
           tip.color = wh1cols,
           align.tip.label = FALSE, type="phylogram", cex = 8)

# Add the SNP scale
add.scale.bar(x=10, y = 1.5, cex = 6, lwd = 15)
text(x=55, y =1.5, cex = 6, "SNPs")


dev.off()


# Get CE6 indices
cegroup = which(herdNames %in% "Clare _ 6")
ceg1 = cegroup[1]
ceg2 = cegroup[2:6]
ceg3 = NA


# Get LH1 indices
lhgroup = which(herdNames %in% "Louth _ 1")
lhg1 = lhgroup[2]
lhg2 = lhgroup[1]
lhg3 = lhgroup[3:5]


# Get C2 indices
cgroup = which(herdNames %in% "Cork _ 2")
cg1 = cgroup[2]
cg2 = cgroup[1]
cg3 = NA


#### Functions ####

# Function to count amount of isolates and herds per each county
numberSamples <- function(isolates){ 
  
  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()
  countyHerds <- list()
  
  # Split off the county name away from herd number
  for(index in 1:length(isolates)){ 
    
    # Split away the herd number
    countyName <- strsplit(isolates[index], split = " ")[[1]][1]
    herdSampled <- strsplit(isolates[index], split = " ")[[1]][3]
    
    # Check if we have encountered the current county before
    if(is.null(countySamples[[countyName]]) == TRUE){
      
      # add it into the list and count once
      countySamples[[countyName]] <- c(1, 1) 
      countyHerds[[countyName]] <- c(herdSampled)
      
      # if already there add an extra 1 to the count for each time seen
    } else { 
      
      countySamples[[countyName]][1] <- countySamples[[countyName]][1] + 1
      
      if(herdSampled %in% countyHerds[[countyName]] == FALSE){
        countyHerds[[countyName]] <- c(countyHerds[[countyName]], herdSampled)
        countySamples[[countyName]][2] <- length(countyHerds[[countyName]])
      }
    }
  }
  return(countySamples)
}

# Function to get spatial polygons
getSpatialPolygons <- function(counties, path) { # input is county list
  
  polygonCoords <- list() # empty list to store gps data for each county
  
  
  for(index in 1:length(counties)) { # for each county
    
    fileName <- paste(path, "PolygonCoords_", counties[index], ".txt", sep="" )
    # paste general path to counties, and the .txt ending to get path for each
    
    # read in the values for each county into the appropriate list segment
    polygonCoords[[counties[index]]] <- SpatialPolygons(list(Polygons(list(Polygon(read.table(fileName, 
                                                                                              header = TRUE, sep = "\t"))), "x")))
  }
  return(polygonCoords)
}

# Function to calculate map ranges
mapSpatialLimits <- function(polygonCoords, counties) {
  
  # Initialise vectors to store the mins and maxes of each county
  minX <- rep(NA, length(counties))
  maxX <- rep(NA, length(counties))
  minY <- rep(NA, length(counties))
  maxY <- rep(NA, length(counties))
  
  # For each county, calculate the mins and maxes, and populate the empty vectors
  for(index in 1:length(counties)) {
    
    minX[index] <- extent(polygonCoords[[counties[index]]])[1]
    maxX[index] <- extent(polygonCoords[[counties[index]]])[2]
    minY[index] <- extent(polygonCoords[[counties[index]]])[3]
    maxY[index] <- extent(polygonCoords[[counties[index]]])[4]
  }
  
  # Get the overall mins and maxes of X and Y to be able to set plot limits
  ranges <- c(min(minX), max(maxX), min(minY), max(maxY))
  
  return(ranges)
}

# Calculate max samples
sampleMax <- function(sampleNumbers) {
  currentmaxValue <- 1 # set current max = 1
  
  keys <- names(sampleNumbers) # assign names as keys to access list
  
  for(index in 1:length(keys)) { # for every element in the list
    if(sampleNumbers[[keys[index]]][1] > currentmaxValue) { # if they are greater than 1
      currentmaxValue = sampleNumbers[[keys[index]]][1] # current max is changed to that
    }
  }
  
  return(currentmaxValue)
}

# Fill polygons on map plot
polygonsSpatialData <- function(polygonCoords, sampleNumbers, counties, shortCounties, herdnames, yearData, someTotal, cattleherd) {
  
  for(index in 1:length(counties)) {
    
    # Calculate alpha = proportion of animals in current county
    # Cattleherd is either a 2 or a 3 depending on herd or cattle
    proportion <- yearData[grep(counties[index], yearData[,1], ignore.case = TRUE), cattleherd]/someTotal
    
    # Add county name
    if(is.null(sampleNumbers[[counties[index]]]) == TRUE ){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("blue", proportion))
      
    }else{ # counties that are present on the sample list
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("blue", proportion))
      
      # Create a vector of herds in the current county
      currentHerds <- unique(herdnames[grep(counties[index], herdNames)])
      
      # For each herd in the current county herds, count how many isolates and plot point
      for(herd in currentHerds){
        
        pointSize <- length(grep(herd, herdnames))
        
        points(spsample(polygonCoords[[counties[index]]], n=1, type = "random"),
               pch = 21, col = alpha("red", 0.9), bg = alpha("red", 0.65), cex = 5*pointSize)
      }
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
    }
  }
}

# Function to get a table of counts of isolates for each year
getYearCounts <- function(herdnames, matrix){
  
  # Create a new vector of length herdnames
  yearVector <- rep(NA, length(herdnames))
  
  # Loop thru columns and extract years
  for(col in 1:ncol(matrix)){
    
    yearVector[col] <- strsplit(colnames(matrix)[col], split = "-")[[1]][1]
  }
  
  # Paste the date into the herd vector
  comboVector <- paste(yearVector,"-",herdnames)
  
  # Tabulate the counts of each
  herdYears <- table(comboVector)
  
  # Make a matrix of length herds and width years
  countMatrix <- matrix(data = 0, nrow = length(unique(herdnames)), ncol = 5, 
                        dimnames = list(sort(unique(herdnames)), c("2012","2013","2014","2016","2017")))
  
  # Get the name of the table
  herdos <- names(herdYears)
  
  # Loop thru the names
  for(index in 1:length(herdos)){
    
    # Split off the year and herd
    year <- strsplit(herdos[index], split = "-")[[1]][1]
    place <- strsplit(herdos[index], split = "-")[[1]][2]
    
    # Match column
    currentCol <- grep(trimws(year), colnames(countMatrix))
    
    # Match row
    currentRow <- which(trimws(place) == rownames(countMatrix))
    
    # Add a count to the relevant matrix slot
    countMatrix[currentRow, currentCol] <- herdYears[index]
    
  }
  
  return(countMatrix)
}

# Function to read in stats file (author JC)
readTBStatisticsFile <- function(fileName){
  
  # Note that file was downloaded from: https://www.cso.ie/px/pxeirestat/Database/eirestat/Animal%20Disease%20Statistics/Animal%20Disease%20Statistics_statbank.asp?SP=Animal%20Disease%20Statistics&Planguage=0&ProductID=DB_DA
  # I selected all years and counties, downlaoded as csv and then removed quotes
  
  # Open a connection to a file to read (open="r")
  connection <- file(fileName, open="r")
  
  # Get all lines from file and store in vector
  fileLines <- readLines(connection)
  
  # Close file connection
  close(connection)
  
  # Intialise a dataframe to store the TB statistics
  statistics <- NULL
  county <- "NA"
  row <- 0
  
  # Loop through each of the lines in file
  for(i in 4:length(fileLines)){
    
    # Remove quotes
    fileLines[i] <- gsub(pattern="\"", replacement="", x=fileLines[i])
    
    # Split the current line into its parts
    parts <- strsplit(fileLines[i], split=",")[[1]]
    
    # If 24th line get the years initialise a dataframe to store the statistics
    if(i == 4){
      statistics <- data.frame(matrix(nrow=1, ncol=length(parts)), stringsAsFactors=FALSE)
      colnames(statistics) <- c("County", "Statistic", parts[-c(1,2)])
      next
    }
    
    # Check if found new county - name will present alone on new line
    if(length(parts) == 1){
      county <- parts[1]
      next
    }
    
    # Store the statistics from the current line
    row <- row + 1
    statistics[row, ] <- c(county, parts[-1])
  }
  
  # Convert the county names to upper case
  statistics$County <- toupper(statistics$County)
  
  return(statistics)
}

# Function to calculate summary stats (author JC)
calculateSummaryStatisticsPerQuarter <- function(statistics){
  
  # Keep only the statistics of interest
  info <- statistics[grepl(statistics$Statistic, pattern="Herds in County|Animals in County"), ]
  
  # Examine each of the year quarters and store a summary
  output <- list()
  for(column in colnames(info)[-c(1,2)]){
    
    # Initialise a data frame to store a data summary
    dataSummary <- data.frame("County"=NA, "HerdsInCounty"=-1, "AnimalsInCounty"=-1,
                              stringsAsFactors=FALSE)
    
    # Examine that TB statistics of interest
    row <- 0
    for(i in seq(from=2, to=nrow(info), by=2)){
      
      # Increment the row
      row <- row + 1
      
      # Store the information for the current county
      dataSummary[row, "County"] <- info[i, "County"]
      dataSummary[row, "HerdsInCounty"] <- as.numeric(info[i-1, column])
      dataSummary[row, "AnimalsInCounty"] <- as.numeric(info[i, column])
    }
    
    # Combine Tipperary (north & south), Cork (north & south), and Wicklow (east and west)
    corkRow <- which(dataSummary$County == "CORK NORTH")
    dataSummary[corkRow, "County"] <- "CORK"
    dataSummary[corkRow, 2:3] <- dataSummary[corkRow, 2:3] + dataSummary[corkRow + 1, 2:3]
    dataSummary <- dataSummary[-(corkRow + 1), ]
    
    tipperaryRow <- which(dataSummary$County == "TIPPERARY NORTH")
    dataSummary[tipperaryRow, "County"] <- "TIPPERARY"
    dataSummary[tipperaryRow, 2:3] <- dataSummary[tipperaryRow, 2:3] + dataSummary[tipperaryRow + 1, 2:3]
    dataSummary <- dataSummary[-(tipperaryRow + 1), ]
    
    wicklowRow <- which(dataSummary$County == "WICKLOW E")
    dataSummary[wicklowRow, "County"] <- "WICKLOW"
    dataSummary[wicklowRow, 2:3] <- dataSummary[wicklowRow, 2:3] + dataSummary[wicklowRow + 1, 2:3]
    dataSummary <- dataSummary[-(wicklowRow + 1), ]
    
    # Store the data summary for the current quarter
    output[[column]] <- dataSummary
  }
  
  return(output)
}

# Function to plot a movement map
plotMoves <- function(ranges, polygonCoords, counties, shortCounties, birthcounties, currentcounties, title, vntrs, sames, herds){

  # Create empty plot, with input of ranges
  plot(x=NA, y=NA,
       xlim = c(ranges[1], ranges[2]), 
       ylim = c(ranges[3], ranges[4]),
       main = title, xlab = "", ylab = "",
       bty = "n", axes = FALSE, cex.main = 7)

  legend("topleft", legend = c("INMV Count:", paste(names(table(vntrs)),"-",table(vntrs))), cex = 8, bty = "n")
  
  
  legend("bottomright", legend = c(1,2,3,4,5,6), 
         title = "Num. Moves to Herd", bty="n", cex=8,
         pch = c(17,17,17,17,17,17), col = alpha("gold", 0.9), pt.bg = alpha("gold", 0.65), 
         pt.cex = c(13,16,19,21,24,27), ncol = 2)
  
  polygonsMoveData(polygonCoords, counties, shortCounties, birthcounties, currentcounties, sames, herds)
  

}

# Fill polygons on map plot for movement arrow plot
polygonsMoveData <- function(polygonCoords, counties, shortCounties, birthcounties, currentcounties, sames, herds) {
  
  for(index in 1:length(counties)) {
    
    # Add county name
    if(counties[index] %in% unique(currentcounties) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("blue", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
      
    }else if(counties[index] %in% unique(birthcounties) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("red", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
    
    } else{

      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE)      
    }
  }
  
  # Initialise vector to store birth/current combos
  combo <- c()
  
  # Initalise vector to store intra county movements
  intra <- c()
  
  # Loop thru each entry in the names
  for(index in 1:length(currentcounties)){
    
    # Check if there is an NA in the current position in either birth or current
    if(is.na(currentcounties[index]) || is.na(birthcounties[index]) || birthcounties[index] == "n/a"){
      
      next
    } else if(currentcounties[index] == birthcounties[index]){
      
      if(sames[index] == "Same"){
        next
      } else{
        
        intra <- append(intra, paste(herds[index], "-", birthcounties[index])) 
      }
    } else{
      
      # Append the combo into the combo vector
      combo <- append(combo, paste(birthcounties[index], "_", currentcounties[index]))
    }
  }  
  
  # Loop thru the combo vector
  for(item in unique(combo)){
      
    # Weight arrow thickness based on how many occurences of a thing
    weight <- sum(combo == item)
      
    # Split apart the names
    birth <- strsplit(item, split = " _ ")[[1]][1]
    current <- strsplit(item, split = " _ ")[[1]][2]
      
    # Get the rough centre of each polygon
    meanB <- coordinates(polygonCoords[[birth]])
    meanC <- coordinates(polygonCoords[[current]])
      
    arrows(meanB[1], meanB[2], meanC[1], meanC[2], col = alpha("black", 0.4), length = 0.5, lwd = weight*3+15)
  }
  
  # Loop thru the intra vector
  for(thingy in unique(intra)){
    
    # Weigh point based on how many moves to that herd
    intraweight <- sum(intra == thingy)
    
    # Split apart the birth county
    intrabirth <- strsplit(thingy, split = "-")[[1]][2]
    
    # Plot random points to indicate movement to that herd
    points(spsample(polygonCoords[[trimws(intrabirth)]], n=1, type = "random"),
           pch = 17, col = alpha("gold", 0.9), bg = alpha("gold", 0.65), cex = 10+3*intraweight)
  }
}

# Function to find isolates within 10 SNPs of a given group of tips
proxFinder <- function(sometips, allDist){
  
  # Create vectors for the output
  out = c()
  
  # Check if any other isolates are within 10 SNPs of each group by looping thru matrix
  for(i in sometips){
    
    for(col in 1:ncol(allDist)){
      
      if(col %in% sometips){
        
        next
      }else if(col %in% out){
        
        next
        
      }else if(allDist[i,col] <= 10){
        
        out = append(out, col)
      }
    }  
  }
  return(out)
}

# Function to plot labelled arrows for groups
plotGroupArrows <- function(group, birthcountyNames, currentcounties, polygonCoords, allDist, colour){
  
  # Add movement arrows for each group
  for(loc in group){
    
    roll = runif(4, min = -0.15, max = 0.15)
    
    # Get current & birth locations
    birth = birthcountyNames[loc]
    current = currentcounties[loc]
    
    # Get the rough centre of each polygon
    meanB <- coordinates(polygonCoords[[birth]])
    meanC <- coordinates(polygonCoords[[current]])
    
    # Get the name of isolate and remove year prefix
    withprefix <- strsplit(rownames(allDist)[loc], split = "_")[[1]][1]
    noprefix <- strsplit(withprefix, split = "-")[[1]][2]
    
    Arrows(meanB[1]+roll[1], meanB[2]+roll[2], meanC[1], meanC[2], col = alpha(colour, 0.6), length = 5, lwd = 25, arr.type = "circle", arr.length = 2, arr.width = 2)
    midx = ((meanC[1]) - (((meanC[1]) - (meanB[1]+roll[1]))/2))
    midy = ((meanC[2]) - (((meanC[2]) - (meanB[2]+roll[2]))/2))
    m = atan(((meanC[2]) - (meanB[2]+roll[2]))/((meanC[1]) - (meanB[1]+roll[1])))*(180/pi)
    text(x=midx, y = midy, labels = noprefix, cex = 5,srt=m)
  }
}

# Function to plot example herd county and closely related locations
makeExampleHerd <- function(g1, g2, g3, counties, polygonCoords, shortCounties, allDist, currentcounties, birthcountyNames){
  
  # Find closely related isolates to groups
  p1 = proxFinder(g1, allDist)
  p2 = proxFinder(g2, allDist)
  
  if(is.na(g3) == FALSE){
    p3 = proxFinder(g3, allDist)
  }
  
  for(index in 1:length(counties)) {
    
    # Add county name
    if(counties[index] %in% unique(currentcounties[p1]) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("red", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
      
    }else if(counties[index] %in% unique(currentcounties[p2]) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("darkgreen", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
      
    } else if(is.na(g3) == FALSE && counties[index] %in% unique(currentcounties[p3]) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("gold", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
      
    } else {
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE)      
    }
    
    if(counties[index] %in% unique(currentcounties[g1]) == TRUE){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
           border = "black", add = TRUE, col = alpha("blue", 0.5))
      
      # Add county name
      pointLabel(coordinates(polygonCoords[[counties[index]]]), labels = shortCounties[index], cex = 7)
      
    }
  }
  
  # Add movement arrows for each group
  plotGroupArrows(g1, birthcountyNames, currentcounties, polygonCoords, allDist, "red")
  
  if(is.null(p1) == FALSE){
    
    # Add movement arrows for each group
    plotGroupArrows(p1, birthcountyNames, currentcounties, polygonCoords, allDist, "red")
  }
  
  # Add movement arrows for each group
  plotGroupArrows(g2, birthcountyNames, currentcounties, polygonCoords, allDist, "darkgreen")
  
  if(is.null(p2) == FALSE){
    
    plotGroupArrows(p2, birthcountyNames, currentcounties, polygonCoords, allDist, "darkgreen")
    
  }
  
  if(is.na(g3) == FALSE){
    
    plotGroupArrows(g3, birthcountyNames, currentcounties, polygonCoords, allDist, "gold")
    
    if(is.null(p3) == FALSE){
      plotGroupArrows(p3, birthcountyNames, currentcounties, polygonCoords, allDist, "gold")
    }  
  }
  
}

# Function to plot tree snippet for example herd
makeExampleTree <- function(g1, g2, g3, allDist, tree){
 
  # Find closely related isolates to groups
  p1 = proxFinder(g1, allDist)
  p2 = proxFinder(g2, allDist)
  
  if(is.na(g3) == FALSE){
    p3 = proxFinder(g3, allDist)
  } 
  
  # Combine groups to make vector of tips to keep
  if(is.na(g3) == FALSE){
    tipstokeep <- sort(c(g1,g2,g3,p1,p2,p3))
  } else{
    tipstokeep <- sort(c(g1,g2,p1,p2))
  }
  
  
  # Make vector of tips to dump
  tipstodump <- c(1:length(tree$tip.label))
  tipstodump = tipstodump[-tipstokeep]
  
  # Create example tree
  exampletree <- drop.tip(tree, tipstodump)
  
  return(exampletree)
}
