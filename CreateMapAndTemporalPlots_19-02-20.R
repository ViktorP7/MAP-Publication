# Created on 19-02-20
# Script to generate isolate map of Ireland and temporal (supplementary) plot
# This script is NOT independent, you are REQUIRED to run DoPhylogeneticTrees_18-02-20.R before running this one

# Load packages
library(sp)
library(raster)
library(gplots)
library(maptools)

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
outputFile <- paste("IrishMap.png", sep="")
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
outputFile <- paste("IrishMap.pdf", sep="")
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
outputFile <- paste("TemporalIsolates.pdf", sep="")
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
      
      # Get number of samples for current county
      nSamples <- sampleNumbers[[counties[index]]][1]
      nHerdsSamples <- sampleNumbers[[counties[index]]][2]
      
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
