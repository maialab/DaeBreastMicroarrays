library(coin)
library(progress)

# Depends:
#   -- ggplot2, because of function cut_number

tsnp_rsnp_test <- function(tsnp_rsnp_pairs,
                           eratios, genotypes,
                           scenario1_threshold = log2(1.5), 
                           min_sample_homo = 2, 
                           min_sample_het = 2) {
  
  # add two new columns to tsnp_rsnp_pairs by reference: pval and group
  tsnp_rsnp_pairs[, c("pval", "group") := list(NA_real_, NA_integer_)]
  
  f <- function(i) {
    tSNP <- tsnp_rsnp_pairs[[i, "tsnp"]]
    SNP <- tsnp_rsnp_pairs[[i, "rsnp"]]
    
    # if tSNP does not exist in genotypes then return a list of two empty
    # vectors
    if(genotypes[.(tSNP), mult = "first", nomatch = 0L, .N] == 0) {
      # warning(paste("tSNP", tSNP, "is missing from the genotype matrix."))
      # list(het = vector(mode = "numeric", length = 0), homo = vector(mode = "numeric", length = 0))
      # c(pval = NA, group = 1)
      set(tsnp_rsnp_pairs, i = i, j = 6L, value = NA)
      set(tsnp_rsnp_pairs, i = i, j = 7L, value = 1L)
      return()
    }
    
    if(eratios[.(tSNP), mult = "first", nomatch = 0L, .N] == 0) {
      # warning(paste("tSNP", tSNP, "is missing from the DAE matrix."))
      # list(het = vector(mode = "numeric", length = 0), homo = vector(mode = "numeric", length = 0))
      # c(pval = NA, group = 1)
      set(tsnp_rsnp_pairs, i = i, j = 6L, value = NA)
      set(tsnp_rsnp_pairs, i = i, j = 7L, value = 1L)
      return()
    }
    
    SNP_indices <- genotypes[.(SNP), mult = "first", nomatch = 0L][, snp := NULL]
    
    # define indices for homo- and het-samples for SNP
    SNP_homo_indices <- as.vector(SNP_indices == "AA" | SNP_indices == "BB")
    SNP_het_indices <- as.vector(SNP_indices == "AB" | SNP_indices == "BA")
    
    # replace NAs with FALSEs
    SNP_homo_indices[is.na(SNP_homo_indices)] <- FALSE
    SNP_het_indices[is.na(SNP_het_indices)] <- FALSE
    
    tSNP_indices <- genotypes[.(tSNP), mult = "first", nomatch = 0L][, snp := NULL]
    # find which samples are heterozygous for the tSNP
    tSNP_het_indices <- as.logical(tSNP_indices == "AB" | tSNP_indices  == "BA" )
    # replace NAs with FALSEs
    tSNP_het_indices[is.na(tSNP_het_indices)] <- FALSE
    
    
    # check if all values in SNP_homo_indices are FALSE
    if(all(!(tSNP_het_indices & SNP_homo_indices))) {
      tSNP_DAE_SNP_homo = vector(mode = "numeric", length = 0)
    } else {
      #tSNP_DAE_SNP_homo <- eratios[.(tSNP), tSNP_het_indices & SNP_homo_indices, with=FALSE][, snp := NULL]
      tSNP_DAE_SNP_homo <- eratios[.(tSNP)][, snp := NULL][, tSNP_het_indices & SNP_homo_indices, with=FALSE]
    }
    
    # DAE values for het-samples tSNP
    het_indices <- tSNP_het_indices & SNP_het_indices
    if(all(!het_indices)) {
      tSNP_DAE_SNP_het = vector(mode = "numeric", length = 0)
    } else {
      tSNP_DAE_SNP_het <- eratios[.(tSNP)][, snp := NULL][, het_indices, with=FALSE]
    }
    
    # list(het = tSNP_DAE_SNP_het, homo = tSNP_DAE_SNP_homo)
    
    data <- c(as.numeric(tSNP_DAE_SNP_het), as.numeric(tSNP_DAE_SNP_homo))
    n_het <- length(tSNP_DAE_SNP_het)
    n_homo <- length(tSNP_DAE_SNP_homo)
    stratification <- factor(rep(c("het", "homo"), c(n_het, n_homo)))
    
    # list(y = data, x = stratification)
    
    # Case 1: when there are no more than min_sample_het of het samples
    if(n_het < min_sample_het) {
      set(tsnp_rsnp_pairs, i = i, j = 6L, value = NA)
      set(tsnp_rsnp_pairs, i = i, j = 7L, value = 1L)
      return()
    }
    
    # when there are no more than min_sample_homo of homo samples
    if(n_homo < min_sample_homo) {
      # Case 2: Using "scenario 1" definition according to Xiao et al. 2009
      if(all(abs(tSNP_DAE_SNP_het) >= scenario1_threshold)) {
        set(tsnp_rsnp_pairs, i = i, j = 6L, value = NA)
        set(tsnp_rsnp_pairs, i = i, j = 7L, value = 2L)
        return()
      }
      
      # Case 3: scenario 2 or 3, according to Xiao et al. 2009 
      # rSNPs could not explain DAE because all samples are het for rSNPs
      # and are DAE scenario 2 and 3)
      else {
        set(tsnp_rsnp_pairs, i = i, j = 6L, value = NA)
        set(tsnp_rsnp_pairs, i = i, j = 7L, value = 3L)
        return()
      }
    }
    
    # Case 4: bona fide rSNP/tSNP pair to be tested for association
    if(
      length(tSNP_DAE_SNP_het) >= min_sample_het && 
      length(tSNP_DAE_SNP_homo) >= min_sample_homo)
    {
      
      stat_test_result <- oneway_test(data ~ stratification,
                                      alternative = "greater", 
                                      distribution = asymptotic(maxpts = 25000))
      pval <- pvalue(stat_test_result)[1]
      names(pval) <- NULL
      
      set(tsnp_rsnp_pairs, i = i, j = 6L, value = pval)
      set(tsnp_rsnp_pairs, i = i, j = 7L, value = 4L)
      return()
    }
  }
  
  dt_nrow <- nrow(tsnp_rsnp_pairs)
  pb <- progress_bar$new(
    format = "Performing oneway_tests [:bar] :percent | Time left: :eta | Test :current of :total.",
    total = dt_nrow, clear = FALSE, width= 120)
  #for(i in seq_along(pairs_indices)) {
  for(i in 1:dt_nrow) {
    # association_data[[i]] <- f(pairs_indices[i])
    f(i)
    # setTxtProgressBar(pb, i)
    pb$tick()
  }
  return(tsnp_rsnp_pairs)
  
}

tsnp_rsnp_test_parallel <- function(tsnp_rsnp_pairs,
                                    eratios, genotypes,
                                    ncores = NULL,
                                    scenario1_threshold = log2(1.5), 
                                    min_sample_homo = 2, 
                                    min_sample_het = 2) {
  
  #available_cores <- parallel::detectCores()
  
  if(is.null(ncores)) {
    stop("Please specify the number of threads: ncores is NULL!")
  } else {
    no_cores <- as.integer(ncores)
  }
  
  if(!is.integer(no_cores) || length(no_cores) != 1) {
    stop("ncores argument is not valid: ", no_cores)
  }
  
  if(no_cores < 1)
    stop("ncores must be greater than one: ", no_cores)
  
  message("Preparing parallelisation...")
  
  
  dt_nrow <- nrow(tsnp_rsnp_pairs)
  indices <- seq_along(1:dt_nrow)
  #partitioning_factor <- split(indices, sort((indices)%%no_cores))
  partitioning_factor <- sort((indices)%%no_cores)
  # tsnp_rsnp_pairs_chunks is a list of data.table data frames
  tsnp_rsnp_pairs_chunks <- split(tsnp_rsnp_pairs, partitioning_factor)
  
  message("Running...")
  results <- parallel::mclapply(tsnp_rsnp_pairs_chunks,
                                tsnp_rsnp_test,
                                eratios = eratios,
                                genotypes = genotypes,
                                scenario1_threshold = scenario1_threshold,
                                min_sample_homo = min_sample_homo,
                                min_sample_het = min_sample_het,
                                mc.cores = no_cores
  )
  
  bind_rows(results)
}


