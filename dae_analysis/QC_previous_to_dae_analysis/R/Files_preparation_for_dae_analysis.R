#Instal dependencies
install.packages("remotes")
remotes::install_local("~/MergedSNPs.Hsapiens.dbSNP149.GRCh38p7")
remotes::install_local("~/Chr.Rpts.Hsapiens.dbSNP149.GRCh38p7")
library(MergedSNPs.Hsapiens.dbSNP149.GRCh38p7)
library(Chr.Rpts.Hsapiens.dbSNP149.GRCh38p7)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocStyle")
library(BiocStyle)
#install.packages("poisbinom")
library(poisbinom)
install.packages("SNPassoc")
library(SNPassoc)
library(stringr)

#Define paths
inpath <- "dae_analysis/QC_previous_to_dae_analysis/input"
outpath <- "dae_analysis/QC_previous_to_dae_analysis/output"

###### LOAD DATA previously prepared by Roslin Russell and Mark Dunning ######
#data includes rna and dna ratios for all SNPs at all samples (including homozygous samples), AE ratios (ratios at cDNA normalized by ratios at gDNA) and genotype data for all SNPs at all samples for 64 samples

#File with Illumina microarrays genotypes for 64 DNA samples from normal breast tissue
#Col1 is SNP_ID, col 2 to 65 are sample IDs
Genotypes_original_AB = read.table(file <- file.path (inpath, "HumanExon510S_SNPs&Genotype_info_original.txt"), header = TRUE, sep = "")

#File with genotypes at 64 DNA samples from normal breast tissue coded as AA, AB or BB
load (file <- file.path (inpath, "gtypes2.Rdata"))

#File with intensities of gDNA samples for the Red (dnaR) and Green channels (dnaG) for 64 samples
#These files were already processed as described previously in Liu et al. (2012) Bioinformatics;28(8):1102-8. 
load (file <- file.path (inpath, "dnaR.new.Rdata"))
load (file <- file.path (inpath, "dnaG.new.Rdata"))

#File with intensities of cDNA samples for the Red (rnaR) and Green channels (rnaG) for 64 samples
# These files were already processed as described previously in Liu et al. (2012) Bioinformatics;28(8):1102-8. 
load (file <- file.path (inpath, "rnaR.new.Rdata"))
load (file <- file.path (inpath, "rnaG.new.Rdata"))

####### Calculation of Alellic ratios ##########

#Allelic ratios at gDNA
dnaLR <- dnaR.new - dnaG.new

#Alellic ratios at cDNA
rnaLR <- rnaR.new - rnaG.new

#Alelic ratios at cDNA normalized by gDNA (AE ratios used bottom down in the analysis)
rnaLR.dnaLR.ratio2 <- rnaLR - dnaLR

#################################################### Get dbSNP149_info ########################################################
#rename rownames of the files with information on genotypes and allele A and B
row.names.Genotypes_original_AB <- Genotypes_original_AB$SNP_ID
rownames(Genotypes_original_AB) <- row.names.Genotypes_original_AB

### There are SNPs that do not start with rs (such as CNV) and should be deleted from the analysis
row.names.raw <- rownames(rnaLR.dnaLR.ratio2)

#Function to convert snp ids (char) into numeric
#grepl returns TRUE if a string contains the pattern, otherwise FALSE

convert_names_function <- function (original.row.names) {
  #return rs ID without rs and the original ones with no rs become NAs
  as.numeric(as.character(str_match (original.row.names, "^rs(\\d+)")[,2]))
}

snp.names <- convert_names_function (row.names.raw)

#Use function current_snp () from MergedSNPs.Hsapiens.dbSNP149.GRCh38p7 package to convert old rs ID into new ones (from times to times rsIDs are merged into new ones)
current.snps.ids <- current_snp(snp.names, keep.orig.snp.if.not.merged = TRUE)
unique.current.snps.ids <- unique(current.snps.ids)
#To put rs again before the numbers
#Ramiro created this functin because the way I was doing previously turned the rs numbers into scientific numbers when the number of digits was very high
snp_numeric2snp_character <- function (snps, add.rs = TRUE) {
  if(add.rs) {
    snps_char <- sprintf("rs%.0f", snps)
    # sanitize NAs
    snps_char[snps_char == "rsNA"] <- NA
  }
  else {
    snps_char <- sprintf("%.0f", snps)
    snps_char[snps_char == "NA"] <- NA
  }
  return(snps_char)
}

current.snp.names.with.rs <- snp_numeric2snp_character (current.snps.ids)

unique.current.snp.ids <- unique (current.snps.ids)

new.names.correspondence <- data.frame (original.name = row.names.raw, recent.name = current.snp.names.with.rs)
col_name_to_replace <- which(names(new.names.correspondence) == "recent.name")
names(new.names.correspondence)[col_name_to_replace] <- "db149.name"

#I will save the 9 original files with extra information from these two columns (#old name in microarray annotation and new name according to dbSNP149)
#In all files the SNPs have the same order except for Genotypes_original_AB

rnaLR.dnaLR.ratio2.dbSNP149 = cbind (new.names.correspondence, rnaLR.dnaLR.ratio2)
row.names(rnaLR.dnaLR.ratio2.dbSNP149) <- NULL 
rnaLR.dnaLR.ratio2.dbSNP149$db149.name = as.character(rnaLR.dnaLR.ratio2.dbSNP149$db149.name)
#write.csv (rnaLR.dnaLR.ratio2.dbSNP149, file <- file.path (outpath, "rnaLR.dnaLR.ratio2.dbSNP149"), row.names = FALSE)


rnaLR.dbSNP149 = cbind (new.names.correspondence, rnaLR)
row.names(rnaLR.dbSNP149) <- NULL 
rnaLR.dbSNP149$db149.name = as.character(rnaLR.dbSNP149$db149.name)

dnaLR.dbSNP149 = cbind (new.names.correspondence, dnaLR)
row.names(dnaLR.dbSNP149) <- NULL 
dnaLR.dbSNP149$db149.name = as.character(dnaLR.dbSNP149$db149.name)

gtypes.dbSNP149 = cbind (new.names.correspondence, gtypes2)
rownames(gtypes.dbSNP149) <- NULL 
gtypes.dbSNP149$db149.name = as.character(gtypes.dbSNP149$db149.name)

rnaR.new.dbSNP149 = cbind (new.names.correspondence, rnaR.new)
rownames(rnaR.new.dbSNP149) <- NULL 
rnaR.new.dbSNP149$db149.name = as.character(rnaR.new.dbSNP149$db149.name)

rnaG.new.dbSNP149 = cbind (new.names.correspondence, rnaG.new)
rownames(rnaG.new.dbSNP149) <- NULL 
rnaG.new.dbSNP149$db149.name = as.character(rnaG.new.dbSNP149$db149.name)

dnaR.new.dbSNP149 = cbind (new.names.correspondence, dnaR.new)
rownames(dnaR.new.dbSNP149) <- NULL 
dnaR.new.dbSNP149$db149.name = as.character(dnaR.new.dbSNP149$db149.name)

dnaG.new.dbSNP149 = cbind (new.names.correspondence, dnaG.new)
rownames(dnaG.new.dbSNP149) <- NULL 
dnaG.new.dbSNP149$db149.name = as.character(dnaG.new.dbSNP149$db149.name)

new_new.names.correspondence <- new.names.correspondence[order(new.names.correspondence$original.name),] 
Genotypes_original_AB <- Genotypes_original_AB[order(Genotypes_original_AB$SNP_ID),] 
Genotypes_original_AB.dbSNP149 = cbind (new_new.names.correspondence, Genotypes_original_AB)
rownames(Genotypes_original_AB.dbSNP149) <- NULL 

Genotypes_original_AB.dbSNP149$db149.name = as.character(Genotypes_original_AB.dbSNP149$db149.name)
column_to_remove <- which (names(Genotypes_original_AB.dbSNP149) == "SNP_ID")
Genotypes_original_AB.dbSNP149 = Genotypes_original_AB.dbSNP149 [,-column_to_remove]

############################## 2nd part #########################################
#####################Extraction of dbSNP149 chromosome reports ################

## Get autossomic chromosome and position information for db149 IDs (Did not include X chromosome SNPs...)
#Remove SNPs with no position info
# Remove SNPs flagged as suspected
# In the final we get info for 463760 SNPs

#use read_chr_rpt () function from package "Chr.Rpts.Hsapiens.dbSNP149.GRCh38p7"

#chr1
chr1_original_data = read_chr_rpt("1")
chr1_our.snps.data = subset (chr1_original_data, snp_id %in% as.character(unique.current.snps.ids))
chr1_our.snps.data_good = chr1_our.snps.data [!is.na (chr1_our.snps.data$chr_position) & chr1_our.snps.data$suspect == 0,]
rm(chr1_original_data)

#chr2
chr2_original_data = read_chr_rpt("2")
chr2_our.snps.data = chr2_original_data[chr2_original_data$snp_id %in% unique.current.snps.ids,]
chr2_our.snps.data_good = chr2_our.snps.data [!is.na (chr2_our.snps.data$chr_position) & chr2_our.snps.data$suspect == 0,]
rm(chr2_original_data)

#chr3
chr3_original_data = read_chr_rpt("3")
chr3_our.snps.data = chr3_original_data[chr3_original_data$snp_id %in% unique.current.snps.ids,]
chr3_our.snps.data_good = chr3_our.snps.data [!is.na (chr3_our.snps.data$chr_position) & chr3_our.snps.data$suspect == 0,]
rm(chr3_original_data)

#chr4
chr4_original_data = read_chr_rpt("4")
chr4_our.snps.data = chr4_original_data[chr4_original_data$snp_id %in% unique.current.snps.ids,]
chr4_our.snps.data_good = chr4_our.snps.data [!is.na (chr4_our.snps.data$chr_position) & chr4_our.snps.data$suspect == 0,]
rm(chr4_original_data)

#chr5
chr5_original_data = read_chr_rpt("5")
chr5_our.snps.data = chr5_original_data[chr5_original_data$snp_id %in% unique.current.snps.ids,]
chr5_our.snps.data_good = chr5_our.snps.data [!is.na (chr5_our.snps.data$chr_position) & chr5_our.snps.data$suspect == 0,]
rm(chr5_original_data)

#chr6
chr6_original_data = read_chr_rpt("6")
chr6_our.snps.data = chr6_original_data[chr6_original_data$snp_id %in% unique.current.snps.ids,]
chr6_our.snps.data_good = chr6_our.snps.data [!is.na (chr6_our.snps.data$chr_position) & chr6_our.snps.data$suspect == 0,]
rm(chr6_original_data)

#chr7
chr7_original_data = read_chr_rpt("7")
chr7_our.snps.data = chr7_original_data[chr7_original_data$snp_id %in% unique.current.snps.ids,]
chr7_our.snps.data_good = chr7_our.snps.data [!is.na (chr7_our.snps.data$chr_position) & chr7_our.snps.data$suspect == 0,]
rm(chr7_original_data)

#chr8
chr8_original_data = read_chr_rpt("8")
chr8_our.snps.data = chr8_original_data[chr8_original_data$snp_id %in% unique.current.snps.ids,]
chr8_our.snps.data_good = chr8_our.snps.data [!is.na (chr8_our.snps.data$chr_position) & chr8_our.snps.data$suspect == 0,]
rm(chr8_original_data)

#chr9
chr9_original_data = read_chr_rpt("9")
chr9_our.snps.data = chr9_original_data[chr9_original_data$snp_id %in% unique.current.snps.ids,]
chr9_our.snps.data_good = chr9_our.snps.data [!is.na (chr9_our.snps.data$chr_position) & chr9_our.snps.data$suspect == 0,]
rm(chr9_original_data)

#chr10
chr10_original_data = read_chr_rpt("10")
chr10_our.snps.data = chr10_original_data[chr10_original_data$snp_id %in% unique.current.snps.ids,]
chr10_our.snps.data_good = chr10_our.snps.data [!is.na (chr10_our.snps.data$chr_position) & chr10_our.snps.data$suspect == 0,]
rm(chr10_original_data)

#chr11
chr11_original_data = read_chr_rpt("11")
chr11_our.snps.data = chr11_original_data[chr11_original_data$snp_id %in% unique.current.snps.ids,]
chr11_our.snps.data_good = chr11_our.snps.data [!is.na (chr11_our.snps.data$chr_position) & chr11_our.snps.data$suspect == 0,]
rm(chr11_original_data)

#chr12
chr12_original_data = read_chr_rpt("12")
chr12_our.snps.data = chr12_original_data[chr12_original_data$snp_id %in% unique.current.snps.ids,]
chr12_our.snps.data_good = chr12_our.snps.data [!is.na (chr12_our.snps.data$chr_position) & chr12_our.snps.data$suspect == 0,]
rm(chr12_original_data)

#chr13
chr13_original_data = read_chr_rpt("13")
chr13_our.snps.data = chr13_original_data[chr13_original_data$snp_id %in% unique.current.snps.ids,]
chr13_our.snps.data_good = chr13_our.snps.data [!is.na (chr13_our.snps.data$chr_position) & chr13_our.snps.data$suspect == 0,]
rm(chr13_original_data)

#chr14
chr14_original_data = read_chr_rpt("14")
chr14_our.snps.data = chr14_original_data[chr14_original_data$snp_id %in% unique.current.snps.ids,]
chr14_our.snps.data_good = chr14_our.snps.data [!is.na (chr14_our.snps.data$chr_position) & chr14_our.snps.data$suspect == 0,]
rm(chr14_original_data)

#chr15
chr15_original_data = read_chr_rpt("15")
chr15_our.snps.data = chr15_original_data[chr15_original_data$snp_id %in% unique.current.snps.ids,]
chr15_our.snps.data_good = chr15_our.snps.data [!is.na (chr15_our.snps.data$chr_position) & chr15_our.snps.data$suspect == 0,]
rm(chr15_original_data)

#chr16
chr16_original_data = read_chr_rpt("16")
chr16_our.snps.data = chr16_original_data[chr16_original_data$snp_id %in% unique.current.snps.ids,]
chr16_our.snps.data_good = chr16_our.snps.data [!is.na (chr16_our.snps.data$chr_position) & chr16_our.snps.data$suspect == 0,]
rm(chr16_original_data)

#chr17
chr17_original_data = read_chr_rpt("17")
chr17_our.snps.data = chr17_original_data[chr17_original_data$snp_id %in% unique.current.snps.ids,]
chr17_our.snps.data_good = chr17_our.snps.data [!is.na (chr17_our.snps.data$chr_position) & chr17_our.snps.data$suspect == 0,]
rm(chr17_original_data)

#chr18
chr18_original_data = read_chr_rpt("18")
chr18_our.snps.data = chr18_original_data[chr18_original_data$snp_id %in% unique.current.snps.ids,]
chr18_our.snps.data_good = chr18_our.snps.data [!is.na (chr18_our.snps.data$chr_position) & chr18_our.snps.data$suspect == 0,]
rm(chr18_original_data)

#chr19
chr19_original_data = read_chr_rpt("19")
chr19_our.snps.data = chr19_original_data[chr19_original_data$snp_id %in% unique.current.snps.ids,]
chr19_our.snps.data_good = chr19_our.snps.data [!is.na (chr19_our.snps.data$chr_position) & chr19_our.snps.data$suspect == 0,]
rm(chr19_original_data)

#chr20
chr20_original_data = read_chr_rpt("20")
chr20_our.snps.data = chr20_original_data[chr20_original_data$snp_id %in% unique.current.snps.ids,]
chr20_our.snps.data_good = chr20_our.snps.data [!is.na (chr20_our.snps.data$chr_position) & chr20_our.snps.data$suspect == 0,]
rm(chr20_original_data)

#chr21
chr21_original_data = read_chr_rpt("21")
chr21_our.snps.data = chr21_original_data[chr21_original_data$snp_id %in% unique.current.snps.ids,]
chr21_our.snps.data_good = chr21_our.snps.data [!is.na (chr21_our.snps.data$chr_position) & chr21_our.snps.data$suspect == 0,]
rm(chr21_original_data)

#chr22
chr22_original_data = read_chr_rpt("22")
chr22_our.snps.data = chr22_original_data[chr22_original_data$snp_id %in% unique.current.snps.ids,]
chr22_our.snps.data_good = chr22_our.snps.data [!is.na (chr22_our.snps.data$chr_position) & chr22_our.snps.data$suspect == 0,]
rm(chr22_original_data)

#Juntar todas as dataframes numa única e guardar
new.gene.anno = rbind (chr1_our.snps.data_good, chr2_our.snps.data_good, chr3_our.snps.data_good, chr4_our.snps.data_good, chr5_our.snps.data_good, chr6_our.snps.data_good, chr7_our.snps.data_good, chr8_our.snps.data_good, chr9_our.snps.data_good, chr10_our.snps.data_good, chr11_our.snps.data_good, chr12_our.snps.data_good, chr13_our.snps.data_good, chr14_our.snps.data_good, chr15_our.snps.data_good, chr16_our.snps.data_good, chr17_our.snps.data_good, chr18_our.snps.data_good, chr19_our.snps.data_good, chr20_our.snps.data_good, chr21_our.snps.data_good, chr22_our.snps.data_good)


#Place "rs" again in the numbers
annotated.snp.names.with.rs = snp_numeric2snp_character (new.gene.anno$snp_id)
annotation.snp.names.with.rs = gsub("^", "rs", new.gene.anno$snp_id)
rownames(new.gene.anno) = annotation.snp.names.with.rs

new.gene.anno = new.gene.anno [, c("chr", "chr_position", "mapped_other_assembly", "local_loci", "avg_het", "allele_orig","gmaf", "mapweight", "snp_type", "chr_hits", "ctg_hits", "total_hits", "ctg_accession", "ctg_version", "ctg_id", "ctg_position", "s.e._het", "max_prob", "validated", "genotypes", "linkouts", "orig_build", "upd_build", "suspect", "clin_sign")]
new.gene.anno$SNP = annotation.snp.names.with.rs

new.gene.anno = new.gene.anno [,c(26, 1:25)]
row.names(new.gene.anno) <- NULL 

write.table(new.gene.anno, file <- file.path (outpath, "new.gene.anno.txt"), row.names = FALSE)

######################################### Rownames alteration according to new rs IDs ######################################################

snp.anno.dbSNP149 <- new.gene.anno
rm (new.gene.anno)

### Step 1#### Start by removing SNPs that are not in the snp.anno.dbSNP149
## most of the SNPs deleted were cnv #
# some were SNPs flagged as suspected in the chromosome reports that were not included in the snp.anno.dbSNP149
#I will place the SNPs ID has row.names
rnaLR.dnaLR.ratio2.dbSNP149.filt1 = subset(rnaLR.dnaLR.ratio2.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)

#Now remove SNPs that are duplicated
rnaLR.dnaLR.ratio2.dbSNP149.filt2 = rnaLR.dnaLR.ratio2.dbSNP149.filt1 [!duplicated(rnaLR.dnaLR.ratio2.dbSNP149.filt1$db149.name),]
rownames (rnaLR.dnaLR.ratio2.dbSNP149.filt2) = rnaLR.dnaLR.ratio2.dbSNP149.filt2$db149.name

#Exclude the two columns that have SNPs ID
cols_to_remove <- which (names (rnaLR.dnaLR.ratio2.dbSNP149.filt2) %in% c("original.name", "db149.name"))
rnaLR.dnaLR.ratio2.dbSNP149.filt2 = rnaLR.dnaLR.ratio2.dbSNP149.filt2 [,c(-cols_to_remove)]

###Now filter the remaining files as well###
# And atribute snp names to rownames

rnaLR.dbSNP149.filt1 = subset(rnaLR.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
rnaLR.dbSNP149.filt2 = rnaLR.dbSNP149.filt1 [!duplicated(rnaLR.dbSNP149.filt1$db149.name),]
rownames (rnaLR.dbSNP149.filt2) = rnaLR.dbSNP149.filt2$db149.name
rnaLR.dbSNP149.filt2 = rnaLR.dbSNP149.filt2 [,c(-cols_to_remove)]

dnaLR.dbSNP149.filt1 = subset(dnaLR.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
dnaLR.dbSNP149.filt2 = dnaLR.dbSNP149.filt1 [!duplicated(dnaLR.dbSNP149.filt1$db149.name),]
rownames (dnaLR.dbSNP149.filt2) = dnaLR.dbSNP149.filt2$db149.name
dnaLR.dbSNP149.filt2 = dnaLR.dbSNP149.filt2 [,c(-cols_to_remove)]

gtypes.dbSNP149.filt1 = subset(gtypes.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
gtypes.dbSNP149.filt2 = gtypes.dbSNP149.filt1 [!duplicated(gtypes.dbSNP149.filt1$db149.name),]
rownames (gtypes.dbSNP149.filt2) = gtypes.dbSNP149.filt2$db149.name
gtypes.dbSNP149.filt2 = gtypes.dbSNP149.filt2 [,c(-cols_to_remove)]

rnaG.new.dbSNP149.filt1 = subset(rnaG.new.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
rnaG.new.dbSNP149.filt2 = rnaG.new.dbSNP149.filt1 [!duplicated(rnaG.new.dbSNP149.filt1$db149.name),]
rownames (rnaG.new.dbSNP149.filt2) = rnaG.new.dbSNP149.filt2$db149.name
rnaG.new.dbSNP149.filt2 = rnaG.new.dbSNP149.filt2 [,c(-cols_to_remove)]

rnaR.new.dbSNP149.filt1 = subset(rnaR.new.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
rnaR.new.dbSNP149.filt2 = rnaR.new.dbSNP149.filt1 [!duplicated(rnaR.new.dbSNP149.filt1$db149.name),]
rownames (rnaR.new.dbSNP149.filt2) = rnaR.new.dbSNP149.filt2$db149.name
rnaR.new.dbSNP149.filt2 = rnaR.new.dbSNP149.filt2 [,c(-cols_to_remove)]

dnaR.new.dbSNP149.filt1 = subset(dnaR.new.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
dnaR.new.dbSNP149.filt2 = dnaR.new.dbSNP149.filt1 [!duplicated(dnaR.new.dbSNP149.filt1$db149.name),]
rownames (dnaR.new.dbSNP149.filt2) = dnaR.new.dbSNP149.filt2$db149.name
dnaR.new.dbSNP149.filt2 = dnaR.new.dbSNP149.filt2 [,c(-cols_to_remove)]

dnaG.new.dbSNP149.filt1 = subset(dnaG.new.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
dnaG.new.dbSNP149.filt2 = dnaG.new.dbSNP149.filt1 [!duplicated(dnaG.new.dbSNP149.filt1$db149.name),]
rownames (dnaG.new.dbSNP149.filt2) = dnaG.new.dbSNP149.filt2$db149.name
dnaG.new.dbSNP149.filt2 = dnaG.new.dbSNP149.filt2 [,c(-cols_to_remove)]

Genotypes_original_AB.dbSNP149.filt1 = subset(Genotypes_original_AB.dbSNP149, db149.name %in% snp.anno.dbSNP149$SNP)
Genotypes_original_AB.dbSNP149.filt2 = Genotypes_original_AB.dbSNP149.filt1 [!duplicated(Genotypes_original_AB.dbSNP149.filt1$db149.name),]
rownames (Genotypes_original_AB.dbSNP149.filt2) = Genotypes_original_AB.dbSNP149.filt2$db149.name
Genotypes_original_AB.dbSNP149.filt2 = Genotypes_original_AB.dbSNP149.filt2 [,c(-cols_to_remove)]


######################################## Quality control #######################

rownames(snp.anno.dbSNP149) = snp.anno.dbSNP149$SNP

###################### 1st QC filter - Intensity greater than 9.5 to guarantee a minimum of expression ####################

rnaA = (rnaR.new.dbSNP149.filt2 + rnaG.new.dbSNP149.filt2) / 2 #The average intensity A = (logR + logG) / 2
rnaA.mean = apply(rnaA, 1, mean)
#filter on mean log2 intensity for each snp is greater or equal to 9.5
rnaA.mean.9.5 = which(rnaA.mean >= 9.5)
snps.qc1 = rownames(rnaA[rnaA.mean.9.5,])

rnaLR.dnaLR.ratio2.qc1 = rnaLR.dnaLR.ratio2.dbSNP149.filt2[snps.qc1,]

################################ 2nd QC filter - Which SNPs have more than 4 (at least 5) het values? #########################################

cat("\nWhich SNPs have more than 4 hetz?\n")
flag = 1; #change flag to = 1 only if you actually want to run this code
if(flag == 1){
  het.count = NULL  
  het.count2 = NULL
  for(i in 1:length(rownames(rnaLR.dnaLR.ratio2.qc1))){
    #print (i)
    snp.rs.id = rownames(rnaLR.dnaLR.ratio2.qc1)[i]
    #get all genotypes and rnaLR-dnaLR normalised values for the snp
    gtypes.snp = gtypes.dbSNP149.filt2[snp.rs.id,]
    #check if any of the genotypes are 'NC' and remove them
    nc.cols = which(gtypes.snp == "NC")
    if(length(nc.cols) > 0){
      nc.cols = which(gtypes.snp == "NC")
      gtypes.snp = gtypes.snp[-nc.cols]
    }   
    # how many het values and make a note if it's greater than 4
    if(length(gtypes.snp) > 1){
      het.cols = which(gtypes.snp == "AB")   
      het.count2[i] = length(het.cols)
      if(length(het.cols) > 4){
        het.count[i] = TRUE
      }else{ het.count[i] = FALSE }      
    }else{ 
      het.count[i] = FALSE
      het.count2[i] = 0
    }  
  }
  save (het.count, file = file.path (outpath, "het.count.RData"))
  save (het.count2, file = file.path (outpath, "het.count2.RData"))
}

#load (file <- file.path (outpath, "het.count.RData"))
#load (file <- file.path (outpath, "het.count2.RData"))


###############

rnaLR.dnaLR.ratio2.qc2 = rnaLR.dnaLR.ratio2.qc1[het.count,]


#What's the distribution of Het counts????
pdf(file = file.path (outpath, "hist-het.count.pdf"))
hist(het.count2)
dev.off()

qc2.snp.list = row.names (rnaLR.dnaLR.ratio2.qc2)

############################# 3rd filter - T-tests ###########################
##To check for allelic descrimination in cDNA #
#We want the heterozygous groups to descriminate from the homozygous

rnaLR.qc2 = rnaLR.dbSNP149.filt2[qc2.snp.list,]
gtypes.qc2 = gtypes.dbSNP149.filt2[qc2.snp.list,]


# AAs vs ABs AND AB vs BB

test1.new = NULL #use this to store all the p-values so qvalue can be performed i.e. an fdr multiple testing correction method
test2.new = NULL

for(i in 1:nrow(rnaLR.qc2)){
  print(i)
  aas = which(as.character(gtypes.qc2[i,]) == 1)
  bbs = which(as.character(gtypes.qc2[i,]) == 3)
  abs = which(as.character(gtypes.qc2[i,]) == 2)
  
  if (length(aas) > 2 & length(bbs) > 2 & length(abs) > 2){
    test1.new[i] = t.test(rnaLR.qc2[i,aas], rnaLR.qc2[i,abs], alternative = "greater")$p.value
    test2.new[i] = t.test(rnaLR.qc2[i,abs], rnaLR.qc2[i,bbs], alternative = "greater")$p.value
  } else {
    if(length(aas) >2 & length(bbs) <=2 & length (abs) > 2){
      test1.new[i] = t.test(rnaLR.qc2[i,aas], rnaLR.qc2[i,abs], alternative = "greater")$p.value
    } else {   
      if(length(aas) <=2 & length(bbs) >2 & length(abs) >2){
        test1.new[i] = t.test(rnaLR.qc2[i,bbs], rnaLR.qc2[i,abs], alternative = "less")$p.value
      }  
    }
  }
}


test1.0.05.rows = which(test1.new <= 0.05)
test2.0.05.rows = which(test2.new <= 0.05)


table = data.frame (qc2.snp.list, test1.new, test2.new)
write.table (table, file <- file.path (output, "qc2.snp.list&ttestResults.txt"))

bad.ttest.rows = which(test1.new > 0.05 | test2.new > 0.05)
rnaLR.dnaLR.ratio2.qc3 = rnaLR.dnaLR.ratio2.qc2[-bad.ttest.rows,]
qc3.snp.list = rownames (rnaLR.dnaLR.ratio2.qc3)

########################## 4th filter: Hardy Weinberg Equilibrium, MAF and No Calls ######################

gtypes.qc3 = gtypes.qc2[qc3.snp.list,]

#To use the SNP assoc package genotypes cannot be AB or BB
#Since for this analysis the true genotypes/alleles don't matter, only if they are homozygous for each of the alleles or heterozygous, I will replace the B by G
#Files transformed will have an extra column with SNP ID
gtypes.snp.transformed <- as.data.frame(sapply(gtypes.qc3,gsub,pattern="AB",replacement="AG"))
gtypes.snp.transformed1 <- as.data.frame(sapply(gtypes.snp.transformed,gsub,pattern="BA",replacement="AG"))
gtypes.snp.transformed2 <- as.data.frame(sapply(gtypes.snp.transformed1,gsub,pattern="NC",replacement="NA"))
gtypes.snp.transformed3 <- as.data.frame(sapply(gtypes.snp.transformed2,gsub,pattern="BB",replacement="GG"))
gtypes.snp.transformed3$SNP = qc3.snp.list

#Place the column with the SNP ID at the first position
gtypes.snp.transformed3 = gtypes.snp.transformed3[, c(length(colnames(gtypes.snp.transformed3)), 1:length(colnames(gtypes.snp.transformed3))-1)]

#get a transposed dataframe so that it gets SNPassoc input format
row.names(gtypes.snp.transformed3) = gtypes.snp.transformed3$SNP
genotypes.transposed = as.data.frame(t(gtypes.snp.transformed3))

#Exclude the first row that has the SNP IDs (consequence from transposing the dataframe)
genotypes.transposed = genotypes.transposed[-1,]

#We will do Quality Control for Hardy-Weinberg Equilibrium (HWE) deviations using SNPassoc
#Here I will use the original SNP chromosomic and position information althought is not updated anymore since it doesn't matter for the type of quality control I am aplying here
#SNP assoc requires a snp info file. In this case in not important the true chromosome and position.
snp.info.92K = snp.anno.dbSNP149[qc3.snp.list,1:3] #Criar esta lista para todas as listas
colnames (snp.info.92K) = c("SNP", "Chromosome", "Position")

genotypes.transposed$SampleID = row.names(genotypes.transposed)

genotypes.transposed = genotypes.transposed[,c(length(colnames(genotypes.transposed)),1:length(colnames(genotypes.transposed))-1)]

mysnpdata<-setupSNP(data, colSNPs=2:92636, sort=FALSE, info=snp.info.92K, sep="")
hweanalysis<-tableHWE(mysnpdata)
a=summary(mysnpdata)
dim(a)
#[1] 92635     4
write.table (a, file <- file.path (outpath, "QCinfo_92K.txt"),  sep = "\t", quote = FALSE)

###Now I have a dataframe with information on alleles (rsID) ,major.allele.freq, HWE (p-value), and % of missing genotypes###
#I want at least 90% sample calls for each SNP and HWE pvalue higher than 1E-05
d<-a[(a$HWE >0.00001 & a$missing <= 10),] 
dim(d)
#[1] 91546     4

## Filter the 92K file based on SNPs that passed quality control ##
rows.d<-row.names (a[(a$HWE >0.00001 & a$missing <= 10),])
rnaLR.dnaLR.ratio2.qc4 = rnaLR.dnaLR.ratio2.qc3[rows.d,]

qc4.snp.list = row.names (rnaLR.dnaLR.ratio2.qc4)

############## To delete SNPs that are insertions/deletions (indels) #############

bad.vector = c("I", "D")

bad.rows = Genotypes_original_AB.dbSNP149.filt2$Allele1 %in%  bad.vector
snps_DI = rownames(Genotypes_original_AB.dbSNP149.filt2) [bad.rows]

rnaLR.dnaLR.ratio2.qc5 = rnaLR.dnaLR.ratio2.qc4 [!rownames(rnaLR.dnaLR.ratio2.qc4) %in% snps_DI,]
dim(rnaLR.dnaLR.ratio2.qc5)
#[1] 91468    64

rnaLR.dnaLR.ratio2.91K = rnaLR.dnaLR.ratio2.qc5
write.csv (rnaLR.dnaLR.ratio2.91K, file <- file.path (outpath, "rnaLR.dnaLR.ratio2.91K"),  row.names = TRUE)

qc5.snp.list = rownames(rnaLR.dnaLR.ratio2.91K)

gtypes.91K = gtypes.dbSNP149.filt2[qc5.snp.list,]
write.csv (gtypes.91K, file <- file.path (outpath, "gtypes.91K"),  row.names = TRUE)


#we end up with 91K autossomic informative SNPs that will be used to assess DAE.
#let's get gene info for these SNPs according to the chr reports from ncbi

snp.info.91K = snp.anno.dbSNP149[qc5.snp.list,]
write.csv (snp.info.91K, file <- file.path (inpath1, "snp.info.91K"))

length(unique(snp.info.91K$local_loci))
#[1] 19233 unique genes

sessionInfo()

#R Under development (unstable) (2017-02-23 r72249)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 16.04.2 LTS

#Matrix products: default
#BLAS: /opt/R/3.4.0/lib/R/lib/libRblas.so
#LAPACK: /opt/R/3.4.0/lib/R/lib/libRlapack.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#loaded via a namespace (and not attached):
#  [1] compiler_3.4.0 magrittr_1.5   tools_3.4.0    stringi_1.1.2 


sessionInfo()

#R Under development (unstable) (2017-02-23 r72249)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 16.04.2 LTS

#Matrix products: default
#BLAS: /opt/R/3.4.0/lib/R/lib/libRblas.so
#LAPACK: /opt/R/3.4.0/lib/R/lib/libRlapack.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] Chr.Rpts.Hsapiens.dbSNP149.GRCh38p7_0.1.0   MergedSNPs.Hsapiens.dbSNP149.GRCh38p7_0.1.0 stringr_1.2.0                              

#loaded via a namespace (and not attached):
#  [1] compiler_3.4.0 magrittr_1.5   tools_3.4.0    stringi_1.1.2 


sessionInfo()

#R Under development (unstable) (2017-02-23 r72249)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 16.04.2 LTS

#Matrix products: default
#BLAS: /opt/R/3.4.0/lib/R/lib/libRblas.so
#LAPACK: /opt/R/3.4.0/lib/R/lib/libRlapack.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] Chr.Rpts.Hsapiens.dbSNP149.GRCh38p7_0.1.0   MergedSNPs.Hsapiens.dbSNP149.GRCh38p7_0.1.0 stringr_1.2.0                              

#loaded via a namespace (and not attached):
#  [1] compiler_3.4.0 magrittr_1.5   tools_3.4.0    stringi_1.1.2 

sessionInfo()
#R version 3.3.2 (2016-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 16.04.2 LTS

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
#[8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] biomaRt_2.30.0

#loaded via a namespace (and not attached):
#  [1] IRanges_2.8.2        parallel_3.3.2       DBI_0.6-1            tools_3.3.2          RCurl_1.95-4.8       memoise_1.1.0        Rcpp_0.12.11         Biobase_2.34.0       AnnotationDbi_1.36.2
#[10] RSQLite_1.1-2        S4Vectors_0.12.2     BiocGenerics_0.20.0  digest_0.6.12        stats4_3.3.2         bitops_1.0-6         XML_3.98-1.7


