# Genome-wide Association Study of Differential Allelic Expressed loci
# Authors: Joana Xavier <joana.xavier1@gmail.com>, Ramiro Magno <ramiro.magno@gmail.com> and Ana Teresa Maia <maia.anateresa@gmail.com>
# Cancer Functional Genomics Lab <http://www.maialab.org/>
# Description: This script blah blah...(TODO)
# 2017 May 09 | Faro, University of Algarve

library(data.table)
#install.packages("profvis")
library(profvis)
#install.packages("tidyverse")
library(tidyverse)
library(magrittr)
library(compiler)
#install.packages("microbenchmark")
library(microbenchmark)
#install.packages("foreach")
library(foreach)
#install.packages("doParallel")
library(doParallel)
#install.packages("progress")
library(progress)
library(purrr)
#install.packages ("purrrlyr")
library(purrrlyr)
#install.packages ("coin")
library(coin)

# TODO: tidy up package usage... some of these packages are no longer needed.

# Script root path
rootpath <- "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/"
# Input folder: contains all data necessary as input for this script.
inputpath <- file.path(rootpath, "daeQTL_analysis/input/")
# Script folder: R source code
path_scripts <- file.path(rootpath, "daeQTL_analysis/R")

# snp annotation: table of 3 cols: snp, chr and pos.
path_annotation <- file.path(inputpath, "annotation_data")
path_annotation <- file.path(inputpath, "annotation_data_chr5")

#annotation_data_chr5 <- read.csv(file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/annotation_data_chr5", header = TRUE)
#annotation_data_chr5 <- annotation_data_chr5 [,-1]
#write.csv(annotation_data_chr5, file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/annotation_data_chr5", row.names = FALSE)

# list of SNPs (ref Ids) that showed significant differential allelic expression
path_dae_snps <- file.path(inputpath, "dae_snps")
path_dae_snps <- file.path(inputpath, "dae_snps_chr5")

#dae_snps_chr5 <- read.csv(file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/dae_snps_chr5", header = TRUE)
#dae_snps_chr5 <- dae_snps_chr5 [,-1]
#write.csv(dae_snps_chr5, file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/dae_snps_chr5", row.names = FALSE)

# table of genotypes (AA, AB/BA, BB): row is snp, col is sample.
path_genotypes <- file.path (inputpath, "genotype_data")
path_genotypes <- file.path (inputpath, "genotype_data_chr5")

#genotype_data_chr5 <- read.csv(file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/genotype_data_chr5", header = TRUE)
#genotype_data_chr5 <- genotype_data_chr5 [,-1]
#write.csv(genotype_data_chr5, file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/genotype_data_chr5", row.names = FALSE)


# list of SNPs (ref Ids) that are putative daeQTLs
path_reg_snps <- file.path (inputpath, "tested_snps")
path_reg_snps <- file.path (inputpath, "tested_snps_chr5")

#tested_snps_chr5 <- read.csv(file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/tested_snps_chr5", header = TRUE)
#tested_snps_chr5 <- tested_snps_chr5 [,-1]
#write.csv(tested_snps_chr5, file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/tested_snps_chr5", row.names = FALSE)


# log2(RNA/DNA) data: row is snp and columns are samples.
path_snp_eratios <- file.path (inputpath, "dae_ratios")
path_snp_eratios <- file.path (inputpath, "dae_ratios_chr5")

#dae_ratios_chr5 <- read.csv(file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/dae_ratios_chr5", header = TRUE)
#dae_ratios_chr5 <- dae_ratios_chr5 [,-1]
#write.csv(dae_ratios_chr5, file = "/data/jxavier/Science/Analysis/DaeBreastMicroarrays/daeQTL_analysis/input/dae_ratios_chr5", row.names = FALSE)


## Output folder: contains data generated by this script.
outpath <- file.path(rootpath, "daeQTL_analysis/results")
# create output directory if does not exist yet
dir.create(outpath, showWarnings = FALSE)
# set outpath as working directory
setwd(outpath)


# list of SNPs that show differential allelic expression, i.e., daeSNPs that here are called tSNPs
tsnps <- readr::read_delim (file = path_dae_snps,
                            delim = " ",
                            col_names = "snp",
                            col_types = "c",
                            skip = 1)

# list of potential eQTL SNPs, i.e., rSNPs
rsnps <- readr::read_delim (file = path_reg_snps, 
                            delim = " ", 
                            col_names = "snp", 
                            col_types = "c",
                            skip = 1)

# list of genotyped/imputed SNPs and their genomic position (chromosome and position), 
# and added info on which SNPs are tSNPs and which are rSNPs.
snp_chr_pos <- fread(input = path_annotation, sep = ",",
                     header = FALSE, skip = 1, 
                     col.names = c("snp", "chr", "pos"),
                     colClasses = c("character", "integer", "integer"))
setkey(snp_chr_pos, snp)
# add tsnp and rsnp columns (Boolean)
snp_chr_pos[, tsnp := snp %in% tsnps$snp]
snp_chr_pos[, rsnp := snp %in% rsnps$snp]

# save to csv the snp genomic positions
snp_chr_pos_filename = "snp_chr_pos.csv"
fwrite(snp_chr_pos, file = file.path(outpath, snp_chr_pos_filename), 
       col.names = TRUE)


#
# Function: find_nearby_snps: fetches all snps around snp in snp_chr_pos_row.
#
# Arg 1: snp_chr_pos_row: a one-row tibble (a row from snp_chr_pos), for example:
#        # A tibble: 1 × 5
#                snp   chr    pos  tsnp  rsnp
#              <chr> <int>  <int> <lgl> <lgl>
#        1 rs3131972     1 817341 FALSE  TRUE
#
# Arg 2: snp_chr_pos, should look something like this:RPS23
#
# Arg 3: window_size: distance to look on either side of snp in snp_chr_pos_row
#                     effectively, a window of 2*window_size bps centred around
#                     snp_chr_pos_row$pos.
#
# Note: Automatically takes into account in which chromosome the SNPs lie.
#
# TODO: Improve this function by using a data.table approach instead of the
#       functional approach of dplyr.
find_nearby_snps <- function(snp_chr_pos_row, snp_chr_pos, window_size = 500000) {
  q_snp <- snp_chr_pos_row$snp
  q_pos <- snp_chr_pos_row$pos
  q_chr <- snp_chr_pos_row$chr
  dplyr::filter(snp_chr_pos, chr == q_chr & abs(q_pos - pos) <= window_size) %>%
    dplyr::select(chr, snp, pos) %>%
    dplyr::mutate(., tsnp = rep(q_snp, nrow(.))) %>%
    dplyr::mutate(., tsnp.pos = rep(q_pos, nrow(.))) %>%
    dplyr::select(chr, tsnp, tsnp.pos, snp, pos) %>%
    dplyr::rename(rsnp = snp, rsnp.pos = pos)
}

# autossomal chromosomes only, integers.
#This step will be edited in order to proceed only with chr5 that is in the data input
#chromosomes <- as.integer(1:22)
chromosomes <- as.integer(5)

# file to save a table of tSNPs and their respective rSNP neighbours
tsnp_rsnp_neighbours_filename <- "tsnp_rsnp_neighbours.csv"

# find neighbouring rSNPs around tSNPs, per chromosome.
local(
  for(chromosome in chromosomes) {
    
    is_first_chr <- chromosome == chromosomes[1]
    
    message("Finding neighbouring rSNPs in chromosome ", chromosome, "...")
    
    # filter snps based on chromosome, both tSNPs and rSNPS
    tsnps_tmp <- dplyr::filter(snp_chr_pos, tsnp == TRUE & chr == chromosome)
    rsnps_tmp <- dplyr::filter(snp_chr_pos, rsnp == TRUE & chr == chromosome)
    
    # find nearby rSNPs around tSNPs
    purrrlyr::by_row(.d = tsnps_tmp, ..f = find_nearby_snps, 
                     snp_chr_pos = rsnps_tmp, 
                     .labels = FALSE, .collate = "rows") %>% 
      # remove spurious column .row added by by_row()
      dplyr::select(-.row) %>%
      # save to file
      readr::write_csv(x = ., 
                       file = file.path(outpath, tsnp_rsnp_neighbours_filename),
                       col_names = is_first_chr, append = !is_first_chr)
  }
)
# tidy up temporary objects 
#rm(tsnps_tmp)
#rm(rsnps_tmp)
#rm(is_first_chr)

# tsnp-rsnp table, i.e. tsnp-rsnp to be tested for statistical association
tsnp_rsnp_pairs <- data.table::fread(
  input = file.path(outpath, tsnp_rsnp_neighbours_filename),
  sep = ",",
  header = TRUE)


# read-in the genotypes
genotypes <- data.table::fread(input = path_genotypes, sep = ",", header = TRUE)
data.table::setnames(genotypes, "V1", "snp")
data.table::setkey(genotypes, snp)

# read-in expression ratios
snp_eratios <- data.table::fread(input = path_snp_eratios, sep = ",", header = TRUE)
data.table::setnames(snp_eratios, "V1", "snp")
data.table::setkey(snp_eratios, snp)

# convert expression ratios in absolute values to be used in the statistical test between hets and homos
snp_eratios_abs <- abs(snp_eratios[,-1])
snp_eratios_abs <- cbind (snp_eratios$snp, snp_eratios_abs)
data.table::setnames(snp_eratios_abs, "V1", "snp")
data.table::setkey(snp_eratios_abs, snp)

# read-in function definitions related to association testing
source(file.path(path_scripts, "tsnp_rsnp_test.R"))

# test some dummy data
#system.time(dummy <- tsnp_rsnp_test(tsnp_rsnp_pairs,
#                                   snp_eratios_abs, genotypes))

#tsnp_rsnp_pairs_5 = subset(tsnp_rsnp_pairs, tsnp.pos >= 81800000 & tsnp.pos <= 83000000 )
#system.time(dummy2 <- tsnp_rsnp_test_parallel(tsnp_rsnp_pairs, snp_eratios_abs, genotypes, ncores = 20))
system.time(my_dae_assoc_results <- tsnp_rsnp_test_parallel(tsnp_rsnp_pairs, snp_eratios_abs, genotypes, ncores = 7))

#system.time(dummy <- tsnp_rsnp_test(tsnp_rsnp_pairs_5, snp_eratios_abs, genotypes))

table(my_dae_assoc_results$group)
#1        2        3        4 
#2022963    21484   410676 18572521  #with 44K daeSNPs and 500Kb each side

#Keep only group 2 (mapping approach 1) and group 4 (mapping approach 2) results
groups2_4 = subset (my_dae_assoc_results, group == 2 | group == 4)

write.csv(groups2_4, file = "~/daeQTL_analysis_Groups2&4_44K_500KbEachSide")
