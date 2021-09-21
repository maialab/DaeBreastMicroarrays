if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

BiocManager::install("Gviz")
library(Gviz)
library(GenomicRanges)

#paths
inpath1 <- "dae_analysis/dae_plots/input"
inpath2 <- "dae_analysis/dae_analysis/input"

#Read data
##file with the AE ratios for aeSNPs (including maeSNPs and daeSNPs)
rnaLR.dnaLR.ratio2.91K  <- read.csv(file  = file.path (inpath2, "rnaLR.dnaLR.ratio2.91K"), header = TRUE, row.names = "X")
##file with genotypes for aeSNPs (including maeSNPs and daeSNPs)
gtypes.91K <- read.csv (file <- file.path (inpath2, "gtypes.91K"), header = TRUE, row.names = "X")
#file with aeSNP gene annotation
aesnps_annotation <- read.csv (file <- file.path (inpath1, "dae91k.snplist.ENSEMBL_annotation_october2019.csv"), header = TRUE)
#file with daeSNP gene annotation
daesnps_annotation <- read.csv (file <- file.path (inpath1, "dae44k.snplist.ENSEMBL_annotation_october2019.csv"), header = TRUE)
#Prepare files in the correct format

#PrnaLR.dnaLR.ratio2.91K has in the end of the sample name .GType and rnaLR.dnaLR.ratio2.91K does not. The names have to be uniformised
names (rnaLR.dnaLR.ratio2.91K) <- names (gtypes.91K)

##get a transpose of the allelic ratios dataframe to get the AE ratios for each SNP by columns
my.ratios <- as.data.frame(t(rnaLR.dnaLR.ratio2.91K))

##Eliminating the AE ratios (log2 A/B) that don't come from heterozygous samples
rnaLR.dnaLR.ratio2.91K_AB <- rnaLR.dnaLR.ratio2.91K
rnaLR.dnaLR.ratio2.91K_AB [gtypes.91K != "AB"] <- NA

#Calculating the abs AE ratios to plot
rnaLR.dnaLR.ratio2.91K_AB_abs <- abs(rnaLR.dnaLR.ratio2.91K_AB)

#calculating the mean values of the abs AE ratios
rnaLR.dnaLR.ratio2.91K_AB_abs$abs_AE_means <- rowMeans(rnaLR.dnaLR.ratio2.91K_AB_abs[,1:ncol(rnaLR.dnaLR.ratio2.91K_AB_abs)], na.rm=TRUE )
rnaLR.dnaLR.ratio2.91K_AB_abs$aeSNP <- rownames (rnaLR.dnaLR.ratio2.91K_AB_abs)

#keep in the same dataframe info about aeSNP (tsnp in the dataframe) namely chr_name and position and AE ratios abs values previously calculated
rnaLR.dnaLR.ratio2.91K_AB_abs_position <- merge (rnaLR.dnaLR.ratio2.91K_AB_abs, aesnps_annotation[,c("tsnp", "chr_name", "tsnp_position_hg38")], by.x = "aeSNP", by.y ="tsnp")

#Preparation of data to perform a plot with the absolute value of the allelic ratios in a specific region
#the user should change chr, start.pos and end.pos according to the chromossomic region in which he/she is interested to plot (hg38 annotation)
my.chr <- 5
start.pos <- 81900000
end.pos <- 82523802

df_region <- subset (rnaLR.dnaLR.ratio2.91K_AB_abs_position, chr_name == my.chr & tsnp_position_hg38 %in% start.pos:end.pos)

## To perform a plot thar combines the image of the transcripts in a region together with the local mean and abs AE ratios  #####################

#use biomart to get genes and transcripts chr, positions and IDs

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

my.locus.biomart.data1 = getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype", "hgnc_symbol", "strand", "ensembl_transcript_id", "ensembl_exon_id", "transcript_biotype"),
                               filters = c("chromosome_name","start", "end"),
                               values = list(my.chr, start.pos, end.pos), 
                               mart = ensembl)

my.locus.biomart.data2 = getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
                                  filters = c("chromosome_name","start", "end"), 
                                  values = list(my.chr, start.pos, end.pos), 
                                  mart = ensembl)

my.locus.biomart.data <-cbind (my.locus.biomart.data1[,c("ensembl_gene_id", "chromosome_name","gene_biotype", "hgnc_symbol", "strand", "transcript_biotype")], my.locus.biomart.data2 [,c("ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end")])

#in case we just want to plot specific type of transcripts such as "protein_coding", "processed_pseudogene" or "miRNA"
#my.locus.data = my.locus.biomart.data [(my.locus.biomart.data$gene_biotype == "protein_coding" | my.locus.biomart.data$gene_biotype == "miRNA"), ]


#Prepare the data with gene and transcripts in the correct format for gviz
row.length = length(my.locus.biomart.data$ensembl_gene_id)

my.gviz.plot.data = data.frame (chromosome= rep(paste0("chr",my.chr), times = row.length), start = my.locus.biomart.data$exon_chrom_start, end = my.locus.biomart.data$exon_chrom_end, width = abs(my.locus.biomart.data$exon_chrom_end - my.locus.biomart.data$exon_chrom_start), strand = my.locus.biomart.data$strand, feature = my.locus.biomart.data$gene_biotype, gene= my.locus.biomart.data$ensembl_gene_id,  exon= my.locus.biomart.data$ensembl_exon_id, transcript = my.locus.biomart.data$ensembl_transcript_id, symbol = my.locus.biomart.data$hgnc_symbol)

#change the way strands are defined (+ and - instead of 1 and -1)
my.gviz.plot.data [my.gviz.plot.data == "-1" ] <- "-"
my.gviz.plot.data [my.gviz.plot.data == "1" ] <- "+"
my.gviz.plot.data$strand = as.factor (my.gviz.plot.data$strand)


#To prepare the file with the DAE abs mean for gviz
dae.abs.mean = df_region[,c("aeSNP", "abs_AE_means", "tsnp_position_hg38")]
dae.abs.mean$chr = paste0("chr",my.chr) 
dae.abs.mean$start = df_region$tsnp_position_hg38
dae.abs.mean$strand = "*"
dae.abs.mean$width = "1"
dae.abs.mean$end = dae.abs.mean$start
dae.abs.mean = dae.abs.mean [,c("chr", "start", "end", "width", "strand", "abs_AE_means")]


#Prepare tracks for gviz
gtrack <- GenomeAxisTrack()

itrack <- IdeogramTrack(genome = "hg38" , chromosome = paste0("chr",my.chr))

grtrack <- GeneRegionTrack(my.gviz.plot.data, name="Gene Model",  
                           showId =TRUE, showExonId = FALSE, collapseTranscripts = FALSE)

dae_abs_track = DataTrack (dae.abs.mean, name="mean abs AE", genome = "hg38" , chromosome = paste0("chr",my.chr), start = start.pos, end = end.pos, ylim  =c(0, 5),  type = "h")

plotTracks (list(gtrack, dae_abs_track, grtrack), from=start.pos, to=end.pos, showBandId =TRUE, sizes = NULL, background.title ="#40464C")
dev.off()
