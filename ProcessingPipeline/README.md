# WGS Data Processing Pipeline
This folder contains an assortment of scripts used to process the WGS data by calling the various tools described within the publication.
The processing pipeline was originally written by **Dr Joseph Crispell**, and is featured here with permission. 

## All original pipeline scripts can be found here: https://github.com/JosephCrispell/GeneralTools/tree/master/ProcessingPipeline

The pipeline was used in a pair of M. bovis WGS studies prior to MAP work, they are:
Crispell, J., Benton, C.H., Balaz, D., De Maio, N., Ahkmetova, A., Allen, A., Biek, R., Presho, E.L., Dale, J., Hewinson, G., Lycett, S.J., Nunez-Garcia, J., Skuce, R.A., 
Trewby, H., Wilson, D.J., Zadoks, R.N., Delahay, R.J., Kao, R.R., 2019. Combining genomics and epidemiology to analyse bi-directional transmission of Mycobacterium bovis in a multi-host system. Elife 8.
Crispell, J., Zadoks, R.N., Harris, S.R., Paterson, B., Collins, D.M., de-Lisle, G.W., Livingstone, P., Neill, M.A., Biek, R., Lycett, S.J., Kao, R.R., 
Price-Carter, M., 2017. Using whole genome sequencing to investigate transmission in a multi-host system: bovine tuberculosis in New Zealand. BMC Genomics 18, 180.
 
# Requirements
- Ubuntu 18
- fastqc
- cutadapt
- bwa
- samtools
- blastn
- bcftools
- java jdk
- perl

# Instructions

## Setup
Before running the pipeline, it is necessary to download and index the MAP K10 reference genome. You can find the reference genome here: https://www.ncbi.nlm.nih.gov/nuccore/NC_002944.2
Index the genome with bwa by typing the following:

``` bash
bwa index reference.fasta
```

You also need to download any paired sequence files you wish to run. Normally, these will be paired fastq.gz files (forward and reverse), but when downloading from Illumina BaseSpace,
the forward and reverse may sometimes be split up into two files each (four files per isolate in total). If this is the case, run the script below in the folder containing your files 
to combine them into a standard one forward and one reverse file:

``` bash
bash CombineFASTQsAcrossLanes_09-07-18.sh
```  

 ## Step 1 - Examine FASTQ file quality and define trimming parameters

Move to directory containing forward and reverse FASTQ files
Run `fastqc` using:
``` bash
fastqc --threads 4 *.fastq.gz
```

Examine each of the html files corresponding to each of the FASTQ files
Decide upon appropriate trimming parameters and enter them into a text file 
called PrinseqSettings_DATE.txt

## Step 2 - Align the FASTQ files against reference and create VCF files
Move to directory containing the forward and reverse FASTQ files and PrinseqSettings_DATE.txt file
Run `ProcessRawReads_29-07-19.sh` using:
``` bash
bash ProcessRawReads_29-07-19.sh fastq.gz cutadapt reference.fasta PickRandomReadsFromSAM_13-07-17.pl ExamineBlastOutput_13-07-17.pl
```
*NOTE:* Set the trimming parameters for cutadapt by creating a file called `Cutadapt_[DATE].txt` with the following structure:
```
UNIVERSALADAPTER1="AGATCGGAAGAG" # Illumina universal adapter sequence (source: https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
UNIVERSALADAPTER2="CTGTCTCTTATA" # Nextera Transposase sequence (source:https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
LENGTH=75                        # The minimum length of read to be accepted
QUAL=25                          # Trim low-quality bases from 5' and 3' ends
PROPN=0.5                        # Discard reads with more than this proportion of Ns
TRIML=15                         # Trim sequence at the 5' end by x positions
TRIMR=5                          # Trim sequence at the 3' end by x positions
```

## Step 3 - Merge the VCF files together
Move to directory containing the VCF files resulting from the previous step - "vcfFiles"
Run `MergeVCFFiles_10-01-19.jar` using:
```
java -Xmx22g -jar MergeVCFFiles_10-01-19.jar . reference.gff3
```

*NOTE:* The -Xmx22g parameter sets the Java heap size to 22gb - if you are running a large amount of VCF files with the merging tool, you will need to assign a large heap.
Please ensure your PC has sufficient RAM to assign to heap, as 22gb is beyond most ordinary computers.  

## Step 4 - Examine genome coverage of isolates
Run `ExamineGenomeCoverage_03-07-17.pl` using:
```
perl ExamineGenomeCoverage_03-07-17.pl ReadDepthThreshold genomeCoverage_DATE.txt
```

## Step 5 - Filter the merged VCF file
Run `FilterVariants_15-09-17.pl` using:
```
perl FilterVariants_DATE.pl 1 ReadDepth HighQualityBaseDepth MappingQuality ProportionReadsSupportingAllele SiteCoverageAcrossIsolates GeneralQuality FQ merged_DATE.txt
```

## Step 6 - Rescuing site information for poorer quality isolates
Run `RescueVariantPositionInfo_18-09-17.pl` using:
```
perl RescueVariantPositionInfo_18-09-17.pl nIsolatesAlleleMustBePresentIn HighQualityBaseDepth ProportionReadsSupportingAllele filtered_DATE.txt
```

## Step 7 - Calculate isolate coverage in filtered OR filtered-rescued file
Run `CalculateIsolateCoverageFromFiltered_20-09-17.pl` using:
```
perl CalculateIsolateCoverageFromFiltered_20-09-17.pl filtered_DATE.txt
```

## Step 8 - Create FASTA file from filtered
Run `CreateFastaFromFiltered_28-06-17.pl` using:
```
perl CreateFastaFromFiltered_28-06-17.pl 1 MinNumberSitesBetweenVariantPositions reference.fasta filtered-rescued_DATE.txt
```

## Tree-building
The resulting .fasta file can now be used as an input for RAxML. RAxML can be run via command line, or by running the R script `RunRAXML_27-06-18` 