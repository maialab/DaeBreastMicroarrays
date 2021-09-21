library(purrr)
library(ggplot2)

snps_with_dae_plot <- function(number_of_tsnps_dataframe, 
                               ind_variable1, ind_variable2, ind_variable3,
                               ind_variable2_range, ind_variable3_scalar, 
                               ylabel_left = "No. of tSNPs w/ DAE",
                               ylabel_right = "No. of tSNPs w/ DAE (%)",
                               xlabel = ind_variable1,
                               ind_variable2_label = ind_variable2,
                               xrange = c(NA_real_, NA_real_),
                               yrange = c(NA_real_, NA_real_),
                               ggplot2_extra = NULL,
                               tsnps_max_100 = NULL) {

  if(length(unique(c(ind_variable1, ind_variable2, ind_variable3))) != 3)
    stop("ind_variable{1,2,3} must be all different!")
  
  ind_variable_opts <- c("dae_thr", "number_het_dae_thr", "prop_het_dae_thr")
  
  if(!all(c(ind_variable1, ind_variable2, ind_variable3) %in% ind_variable_opts))
    stop("ind_variable{1,2,3} are not one of ", paste(ind_variable_opts, collapse = ", "), ", but are instead: ",
         paste(c(ind_variable1, ind_variable2, ind_variable3), collapse = ", "), ".")

  number_of_tsnps_dataframe2 <-
  number_of_tsnps_dataframe %>%
  filter_(paste(ind_variable3, "==", ind_variable3_scalar)) %>%
    filter_(paste(ind_variable2, "%in%", paste0("c(", paste(ind_variable2_range, collapse = ","),")")))
  
  if(is.null(tsnps_max_100)) {
    tsnps_max <- max(number_of_tsnps_dataframe2$number_of_tsnps)
  } else {
    tsnps_max <- tsnps_max_100
  }
  
  ggplot() +
    geom_line(data = number_of_tsnps_dataframe2,
              aes_string(x = ind_variable1, y = "number_of_tsnps", col = paste("factor(", ind_variable2, ")"))) +
    xlab(xlabel) +
    scale_y_continuous(ylabel_left,
                       limits = yrange,
                       breaks = ,
                       sec.axis = sec_axis(~ . * 100 / tsnps_max, name = ylabel_right)) +
    labs(colour = ind_variable2_label) +
    xlim(xrange) +
    theme(aspect.ratio = 1./1.65)

}
