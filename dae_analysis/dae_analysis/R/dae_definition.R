library(data.table)

common_snps <- function(genotypes, eratios) {

  # find common snps to both genotype and eratios data frames
  snps <- intersect(eratios$snp, genotypes$snp)
  if(length(snps) == 0) 
    stop("No common snps to genotypes and eratios!")
  
  snps
}


# Convert data frames: genotypes and eratios to matrices
dt2matrices <- function(genotypes, eratios) {
  
  # find common snps to both genotype and eratios data frames
  tsnps <- intersect(eratios$snp, genotypes$snp)
  if(length(tsnps) == 0) 
    stop("No common snps to genotypes and eratios!")
  
  # convert data frames to matrices and remove the 'snp' column
  genotypes2 <- as.matrix(genotypes[.(tsnps)][, snp := NULL])
  eratios2   <- as.matrix(eratios[.(tsnps)][, snp := NULL])
  
  list(tsnps = tsnps, genotypes = genotypes2, eratios = eratios2)
}


# colnames in genotypes and eratios must be the same
# first column must be 'snp'
summary_of_het_samples_showing_dae <- function(dae.threshold, 
                                               genotypes,
                                               eratios, 
                                               comparison = "abs-greater") {
  
  dae_comparison_opts <- c("abs-greater", "abs-less", "greater", "less", "none")
  if(!(comparison %in% dae_comparison_opts)) 
    stop("Argument 'comparison' is not one of: ", dae_comparison_opts)
  
  m <- dt2matrices(genotypes = genotypes, eratios = eratios)
  
  # bonafide samples: samples for which there is simultaneously genotype and
  # expression data, i.e., not missing (NA)
  bonafide_logical <- !(is.na(m[["genotypes"]]) | is.na(m[["eratios"]]))
  
  # heterozygous' samples
  het_logical <- m[["genotypes"]] == "AB" | m[["genotypes"]] == "BA"

  # differential allelic expressing samples
  dae_logical <- switch(comparison,
                        "abs-greater" = abs(m[["eratios"]]) >= dae.threshold,
                        "abs-less" = abs(m[["eratios"]]) <= dae.threshold,
                        "greater" = m[["eratios"]] >= dae.threshold,
                        "less" = m[["eratios"]] <= dae.threshold,
                        "none" = !is.na(m[["eratios"]]) # every expression is 
                        # converted to TRUE, except NA values that remain NA.
  )
  
  # heterozygous samples showing differential allelic expression
  het_dae_logical <- het_logical & dae_logical
  
  # number of bonafide samples
  no_bonafide_samples <- apply(bonafide_logical, 1, sum)
  
  # number of samples
  no_samples <- ncol(m[["genotypes"]]) # equivalent to ncol(m[["eratios"]])
  
  # number of heterozygous samples
  no_het_samples <- apply(het_logical, 1, sum, na.rm = TRUE)
  
  # number of heterozygous samples showing differential allelic expression
  no_het_dae_samples <- apply(het_dae_logical, 1, sum, na.rm = TRUE)
  
  # proportion of heterozygous samples amongst samples for which there is
  # simultaneously genotype and expression data (bonafide samples)
  prop_het_samples <- no_het_samples / no_bonafide_samples
  
  # proportion of heterozygous samples showing differential allelic expression
  # amongst samples for which there are simultaneously genotype and
  # expression data (bonafide samples)
  prop_het_dae_samples <- no_het_dae_samples / no_bonafide_samples
  
  # proportion of heterozygous samples showing differential allelic expression
  # amongst heterozygous samples
  prop_het_dae_het_samples <- no_het_dae_samples / no_het_samples
  
  results <- data.table(tsnp = m[["tsnps"]],
                        number_samples = no_samples,
                        number_bonafide_samples = no_bonafide_samples,
                        number_het_samples = no_het_samples,
                        number_het_dae_samples = no_het_dae_samples,
                        proportion_het_samples = prop_het_samples,
                        proportion_het_dae_samples = prop_het_dae_samples,
                        proportion_het_dae_het_samples = prop_het_dae_het_samples
                        )
  
  results
  
}

summary_filter <- function(datatable,
                           number_het_dae_threshold,
                           proportion_het_dae_het_threshold) {
  subset(datatable, 
         number_het_dae_samples > number_het_dae_threshold &
           proportion_het_dae_het_samples > proportion_het_dae_het_threshold)
}

number_of_tsnps_with_dae <- function(datatable,
                                     number_het_dae_threshold,
                                     proportion_het_dae_het_threshold) {

  sum(datatable$number_het_dae_samples >= number_het_dae_threshold &
    datatable$proportion_het_dae_het_samples >= proportion_het_dae_het_threshold)
  
}
