library(beeswarm)
library(data.table)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(patchwork)
library(ggpubr)
library(lubridate)
library(latex2exp)

#paths
inpath1 <- "dae_analysis/dae_plots/input"
inpath2 <- "dae_analysis/dae_analysis/input"
outpath <- "dae_analysis/dae_plots/plots"

#input files
##file that has correspondence of A and B to the real allele
file.alleles.info.91k <- read.csv(file  = file.path (inpath1, "file.alleles.info.91k.csv"), header = TRUE)
##file with the AE ratios for aeSNPs (including maeSNPs and daeSNPs)
rnaLR.dnaLR.ratio2.91K  <- read.csv(file  = file.path (inpath2, "rnaLR.dnaLR.ratio2.91K"), header = TRUE, row.names = "X")
##file with genotypes for aeSNPs (including maeSNPs and daeSNPs)
gtypes.91K <- read.csv (file <- file.path (inpath2, "gtypes.91K"), header = TRUE, row.names = "X")
#file with aeSNP gene annotation
aesnps_annotation <- read.csv (file <- file.path (inpath1, "dae91k.snplist.ENSEMBL_annotation_october2019.csv"), header = TRUE)
#file with daeSNP gene annotation
daesnps_annotation <- read.csv (file <- file.path (inpath1, "dae44k.snplist.ENSEMBL_annotation_october2019.csv"), header = TRUE)

#Prepare files with correct format
names (rnaLR.dnaLR.ratio2.91K) <- names (gtypes.91K)
rownames (file.alleles.info.91k) <- file.alleles.info.91k$original.name
##get a transpose of the allelic ratios dataframe
my.ratios <- as.data.frame(t(rnaLR.dnaLR.ratio2.91K))
##Eliminating the AE ratios (log2 A/B) that don't come from heterozygous samples
rnaLR.dnaLR.ratio2.91K_AB <- rnaLR.dnaLR.ratio2.91K
rnaLR.dnaLR.ratio2.91K_AB [gtypes.91K != "AB"] <- NA

my.ratios$het_samples <- "AB"
my.ratios$het_samples <- as.factor(my.ratios$het_samples)

#Get a vector with daeSNPs for genes of interest
#The user should change this vector to include the gene symbol of the gene of interest
genes_of_interest <- c("ATG10", "RPS23")

my_dae_snps <- unique(subset (daesnps_annotation, tsnp_hgnc_symbol %in% genes_of_interest)$tsnp)
my_ae_snps <- unique(subset (aesnps_annotation, tsnp_hgnc_symbol %in% genes_of_interest)$tsnp)

#Beeswarm function for aeSNPs AE ratios plots
my_dae_plots_combined_function <- function(my_ratios_data, dae_snps_vector, alleles_data)
  
{
  pdf(paste0(file.path (outpath, "DAE_plot_by_gene_beeswarm_teste.pdf")))
  par(mfrow=c(3,4))
  for (i in 1:length(dae_snps_vector))
    
  {
    
    my_combined_plot <- beeswarm (data = my_ratios_data, my_ratios_data[,dae_snps_vector[i]]~my_ratios_data$het_samples, 
                                  method = c("swarm"), main = NULL, 
                                  ylim = c(-5,5), pch=16, xlab =paste0(dae_snps_vector[i]), 
                                  ylab = "AE ratio",  labels = paste0(alleles_data[dae_snps_vector[i],"Allele1"], alleles_data[dae_snps_vector[i],"Allele2"], sep = ""), cex.axis=0.8)
    abline(h=c(-0.58, 0.58),lty=2,col="black")
    
  }
  return (my_combined_plot)
  dev.off()
}

#Run function to generate AE ratios plots
#decide wether to plot either my_dae_snps (only daeSNPs in the genes of interest) or my_ae_snps (all aeSNPs in the genes of interest)
#the user should modify the second argument in the function accordingly
my_dae_plots_combined_function (my.ratios, my_dae_snps, file.alleles.info.91k)
dev.off()
