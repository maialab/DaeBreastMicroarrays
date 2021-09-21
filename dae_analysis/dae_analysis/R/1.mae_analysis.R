#Define paths

inpath <- "dae_analysis/dae_analysis/input"
outpath <- "dae_analysis/dae_analysis/output"

#Load the data
snp.anno.dbSNP149 <- read.table (file <- file.path (inpath, "new.gene.anno.txt"),  sep ="", header = TRUE)
rownames(snp.anno.dbSNP149) <- snp.anno.dbSNP149$SNP

rnaLR.dnaLR.ratio2.91K <- read.csv (file <- file.path (inpath, "rnaLR.dnaLR.ratio2.91K"), row.names = "X")

gtypes.91K <- read.csv (file <- file.path (inpath, "gtypes.91K"), row.names = "X")

imprintGs <- read.csv(file <- file.path (inpath, "ImprintedGenes_2016.csv"),  header = FALSE, col.names = "ID")
imprintGenes <- unique(imprintGs$ID)

###############################
#Mono-allelic expression analysis

# We will apply a ratio (between alleles) treshold to define DAE
# We have define it as samples showing rnaLR.dnaLR.ratios >= 0.58 or <= -0.58
# Now we will count how many samples display DAE for each SNP

###################################How many hets with DAE there is in each SNP??###################
cat("\nHow many hets?\n")
flag = 1; #change flag to = 1 only if you actually want to run this code
if(flag == 1){
  het.count.91K = NULL  
  DAEnumber.91K = NULL
  DAEpercent.91K = NULL
  for(i in 1:length(rownames(rnaLR.dnaLR.ratio2.91K))){
    print (i)
    snp.rs.id.91K = rownames(rnaLR.dnaLR.ratio2.91K)[i]
    #get all genotypes and rnaLR-dnaLR normalised values for the snp
    gtypes.snp.91K = gtypes.91K[snp.rs.id.91K,]
    het.cols.91K = which(gtypes.snp.91K == "AB")
    hets.DAE = which(rnaLR.dnaLR.ratio2.91K[i,het.cols.91K] >= log2(1.5) | rnaLR.dnaLR.ratio2.91K[i,het.cols.91K] <= -(log2(1.5)))
    het.count.91K [i] = length(het.cols.91K)
    DAEnumber.91K [i]= length (hets.DAE)
    DAEpercent.91K [i] = DAEnumber.91K [i] /het.count.91K [i]
  }  
}
save (het.count.91K, file = file.path (outpath, "het.count.91K.RData"))
save (DAEnumber.91K, file = file.path (outpath, "DAEnumber.91K.RData"))
save (DAEpercent.91K, file = file.path (outpath, "DAEpercent.91K.RData"))

pdf(file <- file.path (outpath, "hist-het.percent.91K.pdf"))
hist(DAEpercent.91K)
dev.off()

pdf(file <- file.path (outpath, "hist-het.count.91K.pdf"))
hist(DAEnumber.91K)
dev.off()

### Separate the list of SNPs with 100% of hets with DAE (either Scenario 1 from Laura Scott paper or mono-alellic expressed SNPs (maeSNPs) ###
DAEpercent.100.rows = which(DAEpercent.91K == 1)
DAEpercent.100.snplist = rownames (rnaLR.dnaLR.ratio2.91K[DAEpercent.100.rows,])
length(DAEpercent.100.snplist)
#[1] 1273

rnaLR.dnaLR.ratio2.DAEpercent.100 = rnaLR.dnaLR.ratio2.91K[DAEpercent.100.snplist, ]

## To separate monoallelic expressed SNPs from SNPs where all the samples display DAE but allways the same allele (DAE scenario 1 as defined by Xiao, Scott 2011 paper#
#We have define DAE as samples showing rnaLR.dnaLR.ratios >= 0.58 or <= -0.58
#In monoallelic expressed snps we expect that half the samples show more expression of one allele and the other hapf more expression of the other allele (from the tSNP)
#I considered as monoallelic expressed snps where there were at least 3 samples displaying DAE in both directions and none sample not showing DAE, since our minimum het samples is 5.
#If there is only one or two samples expressing more of one of the alleles it might be a outlier but we should check manually the plots 

################################## To identify monoallelic expressed genes ###########################################

cat("\nWhich are the imprinted SNPs?\n")
imprint.flag = NULL
for(i in 1:length (DAEpercent.100.snplist)){   
  print (i)
  snp.rs.id.100percent = DAEpercent.100.snplist[i]
  gtypes.snp.100percent = gtypes.91K [snp.rs.id.100percent,]
  het.cols.100percent = which (gtypes.snp.100percent == "AB")   
  greater = which (rnaLR.dnaLR.ratio2.DAEpercent.100 [snp.rs.id.100percent, het.cols.100percent] >= 0.58)
  less = which (rnaLR.dnaLR.ratio2.DAEpercent.100 [snp.rs.id.100percent, het.cols.100percent] <= -0.58)
  if (length (greater) > 2 & length (less) > 2){ imprint.flag[i] = TRUE }
  else{ imprint.flag[i] = FALSE}
}

imprinted.rows = which(imprint.flag == TRUE)

imprinted.snplist = rownames(rnaLR.dnaLR.ratio2.DAEpercent.100[imprinted.rows,])
length(imprinted.snplist)
#[1] 85

save(imprinted.snplist, file = file.path (outpath, "imprinted.snplist.RData"))

scenario1_snps <- setdiff(DAEpercent.100.snplist, imprinted.snplist)
length(scenario1_snps)
#1188
save(scenario1_snps, file = file.path (outpath, "scenario1.snplist.RData"))

#Gene Annotation
snp.info.imprinted = snp.anno.dbSNP149[imprinted.snplist,]
#looking individually at these 85 SNPs
write.csv(snp.info.imprinted, file = file.path (outpath, "snp.info.imprinted.csv"))

#We have to confirm if these SNPs are really located at these genes manually by going to uscs genome browser
mae.genes = as.character(unique(snp.info.imprinted$local_loci))

mae.genes.manual = c("COPG2", "MEG3", "DLK1", "PEG3","ZIM2", "PEG10","KCNQ1", "KCNQ1OT1", "PEG3", "PEG3-AS1", "ZIM2", "NAA60", "ZNF597", "IGF2", "IGF2-AS", "INS-IGF2",  "ZDBF2", "L3MBTL1",  "LOC107987208", "MEG8",  "PEG10", "SGCE", "IGF2", "INS-IGF2",  "LOC107986657", "PLAGL1", "SNORD116-5",  "ZNF331",  "MEG9",  "GNAS", "LOC101927932", "SNORD116-15", "SNORD116-16", "ZIM2", "ZIM2-AS1", "PLAGL1", "SNORD116-17", "SNORD116-18", "SNORD116-29", "MIR127", "MIR432", "RTL1",  "GNAS", "GNAS-AS1", "MIR656", "H19", "HOTS", "MRPL23",  "ARCN1",  "NAA60", "IPW", "MIR154", "MIR496", "ZNF597")

# My maeGenes and known imprinted genes ##
#which ones are novel?

novelImprintGenes = setdiff(mae.genes.manual, imprintGenes) #diz-nos quais estão na 1ª lista que não estão na 2ª
#[1] "PEG3-AS1"     "IGF2-AS"      "INS-IGF2"     "L3MBTL1"      "LOC107987208" "LOC107986657" "SNORD116-5"   "ZNF331"       "MEG9" 
#[10] "LOC101927932" "SNORD116-15"  "SNORD116-16"  "ZIM2-AS1"     "SNORD116-17"  "SNORD116-18"  "SNORD116-29"  "MIR127"       "MIR432"    
#[19] "GNAS-AS1"     "MIR656"       "HOTS"         "MRPL23"       "ARCN1"        "IPW"          "MIR154"       "MIR496"    

############################ To remove imprinted SNPs from the 91K list #################################

rnaLR.dnaLR.ratio2.91K.nonImp= rnaLR.dnaLR.ratio2.91K[ !(rownames(rnaLR.dnaLR.ratio2.91K) %in% imprinted.snplist), ] 
dim(rnaLR.dnaLR.ratio2.91K.nonImp)
#[1] 91383    64

write.csv (rnaLR.dnaLR.ratio2.91K.nonImp, file = file.path (outpath, "rnaLR.dnaLR.ratio2.91K.nonImp.csv"))

snplist.91K.nonImp = rownames (rnaLR.dnaLR.ratio2.91K.nonImp)

#Annotation
snp.info.91K.nonImp = snp.anno.dbSNP149[snplist.91K.nonImp,]
length(unique(snp.info.91K.nonImp$local_loci))
#[1] 19200
write.csv(snp.info.91K.nonImp, file  = file.path (outpath, "snp.info.91K.nonImp.csv"))

gtypes.91K.non.Imp <- gtypes.91K [snplist.91K.nonImp,]
write.csv(gtypes.91K.non.Imp, file  = file.path (outpath, "gtypes.91K.non.Imp.csv"))

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