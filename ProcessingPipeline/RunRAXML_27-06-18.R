# Modded for personal use

# Load the ape library
library(ape)

# Note the FASTA file name and path
fastaFile <- "/home/workspace/vperets/SequencingFiles/sequences_Prox-10_29-01-2020.fasta"

# Set the working directory
setwd("/home/workspace/vperets/SequencingFiles/")

# Build analysis name
analysisName <- "RaxML-R_30-01-20"

# Set the input parameters for RAXML
model <- "GTRCAT" # No rate heterogenity
seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
nBootstraps <- 100
nThreads <- 15

# Build the command for RAXML
command <- paste("raxmlHPC-PTHREADS", # Note on WINDOWS replace with: /path/to/raxmlHPC-PTHREADS.exe
                 " -f a", # Algorithm: Rapid boostrap inference
                 " -N ", nBootstraps,
                 " -T ", nThreads,
                 " -m ", model, " -V", # -V means no rate heterogenity
                 " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                 " -n ", analysisName,
                 " -s ", fastaFile, sep="")

# Run RAXML
system(command, intern=TRUE)

