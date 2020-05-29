#########################################
#                                       #
#   Calculates deltaP from mafs files   #
#                                       #
#   Author: Joshua Penalba              #
#   Date: 22 Oct 2016                   #
#   Downloaded this N. Lou modified ver #
#		03/24/20                        #
#   modified by Matt Hare 03/30/20      #
#########################################


# NOTES
# * Prior to calculating deltaP the following steps are recommended:
#   1. Run ANGSD with all populations with a -SNP_pval and -skipTriallelic flags.
#   2. Rerun ANGSD per population 
#       Use the -sites flag with a file corresponding to the recovered SNPs.
#       This will guarantee that sites with an allele fixed in one population is still included.
#       Remove the -SNP_pval flag.
#       IMPORTANT: Include an outgroup reference to polarize alleles.
#   3. Gunzip the resulting mafs files.
# 
# * Make sure the totLen only includes the chromosomes being analyzed.
# * minInd flag not added, assuming already considered in the ANGSD run.
# * Test for matching major and minor alleles not included as it would filter out sequencing errors. 
#   This has been accounted for in the allele frequency calculations.
#   This filter may give an underestimate of dxy.

# ADDITIONAL NOTES BY NICOLAS LOU
# * This script was slightly modified to allow for:
#   1. gzipped mafs file
#   2. specifying input and output directory
#   3. specifying output name

# example usage: Rscript get_deltaP.R -d /workdir/mph75/angsd -p HH_chr2forDxy_minMQ20ind13dp13_maxdp30_Mar29.mafs.gz -q SV_chr2forDxy_minMQ20ind13dp13_maxdp30_Mar29.mafs.gz -t 255296 -o HHSV_chr2_deltaP_minMQ20ind13dp13_maxdp30

### Creating an argument parser
library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-d","--dir"), type="character",default=NULL,help="path to input and output files",metavar="character"),
  make_option(c("-p","--popA"), type="character",default=NULL,help="path to unzipped mafs file for pop 1",metavar="character"),
  make_option(c("-q","--popB"), type="character",default=NULL,help="path to unzipped mafs file for pop 2",metavar="character"),
  make_option(c("-t","--totLen"), type="numeric",default=NULL,help="total sequence length for global per site Dxy estimate [optional]",metavar="numeric"),
  make_option(c("-o","--out_name"), type="character",default=NULL,help="output file name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

### Set working directory
setwd(opt$dir)

### Troubleshooting input
if(is.null(opt$popA) | is.null(opt$popB)){
  print_help(opt_parser)
  stop("One or more of the mafs paths are missing", call.=FALSE)
}

if(is.null(opt$totLen)){
  print("Total length not supplied. The output will not be a per site estimate.")
}

### Reading data in
allfreqA <- read_tsv(opt$popA)
allfreqB <- read_tsv(opt$popB)

### Manipulating the table and print deltaP table
allfreq <- merge(allfreqA, allfreqB, by=c("chromo","position"))
allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
# -> Actual deltaP calculation
allfreq <- transform(allfreq, deltaP=abs(knownEM.x-knownEM.y))
write.table(allfreq[,c("chromo","position","deltaP")], file=opt$out_name, quote=FALSE, row.names=FALSE, sep='\t')
print('Created deltaP_persite.txt')

### Print global deltaP
print(paste0('Global deltaP is: ',sum(allfreq$deltaP)))
if(!is.null(opt$totLen)){
  print(paste0('Global average deltaP per site is: ',sum(allfreq$deltaP)/opt$totLen))
}