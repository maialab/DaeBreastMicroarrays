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
                                   
#create a vector with SNPs to plot (can include aeSNPs without DAE (non-dae SNPs))
#the user should edit this vector by adding rsIDs of interest to the vector
my_snps <- c("rs11240762", "rs12444974")

###############

#Beeswarm function for aeSNPs AE ratios plots
my_dae_plots_combined_function <- function(my_ratios_data, dae_snps_vector, alleles_data)
  
{
  pdf(paste0(file.path (outpath, "daeSNPs_plot_COQ6_ENTPD5.pdf")))
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


#Run function to generate AE ratios plots for SNPs of interest
my_dae_plots_combined_function (my.ratios, my_snps, file.alleles.info.91k)
dev.off()
