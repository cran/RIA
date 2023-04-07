# #' @title Aggregates GLRLM-based statistics based-on supplied function
# #'
# #' @description Calculates aggregated statistics of GLRLM matrix statistics calculated on
# #' GLRLM matrices evaluated in all different directions.
# #' 
# #' @param RIA_data_in \emph{RIA_image}.
# #' 
# #' @param statistic string, defining the statistic to be calculated on the array of GLRLM statistics.
# #' By default, statistic is set to \emph{"mean"}, however any function may be provided. The proper
# #' syntax is: function(X, attributes). The supplied string must contain a "X", which will be replaced
# #' with the array of the GLRLM statistics value. Further attributes of the function may also be given.
# #' For example, if you wish to calculate the median of all GLRLMs calculated in different directions,
# #' then it must be supplied as: \emph{median(X, na.rm = TRUE)}.
# #'
# #' @param verbose_in logical, indicating whether to print detailed information.
# #' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
# #'
# #' @return \emph{RIA_image} containing the statistical information.
# #'
# #' @examples \dontrun{
# #' #Discretize loaded image and then calculate GLCM statistics for all matrices
# #' RIA_image <- discretize(RIA_image, bins_in = c(4, 8), equal_prob = TRUE,
# #' use_orig = TRUE, write_orig = FALSE)
# #' RIA_image <- glrlm_all(RIA_image, use_type = "discretized")
# #' RIA_image <- glrlm_stat(RIA_image)
# #' 
# #' #Calculate the average of the different GLCM matrices in the different directions
# #' RIA_image <- glrlm_stat_all(RIA_image)
# #' }
# #' 
# #' @references Márton KOLOSSVÁRY et al.
# #' Radiomic Features Are Superior to Conventional Quantitative Computed Tomographic
# #' Metrics to Identify Coronary Plaques With Napkin-Ring Sign
# #' Circulation: Cardiovascular Imaging (2017).
# #' DOI: 10.1161/circimaging.117.006843
# #' \url{https://pubmed.ncbi.nlm.nih.gov/29233836/}
# #' 
# #' Márton KOLOSSVÁRY et al.
# #' Cardiac Computed Tomography Radiomics: A Comprehensive Review on Radiomic Techniques.
# #' Journal of Thoracic Imaging (2018).
# #' DOI: 10.1097/RTI.0000000000000268
# #' \url{https://pubmed.ncbi.nlm.nih.gov/28346329/}
# #' @encoding UTF-8


glrlm_stat_all <- function(RIA_data_in, statistic = "mean(X, na.rm = TRUE)", verbose_in = TRUE)
{
  if(!any(class(RIA_data_in) == "RIA_image")) {message("PROCESSING OF RIA_image OBJECTS ARE SUPPORTED, OTHER CLASSES MIGHT CAUSE PROBLEMS! PLEASE LOAD DATA USING RIA load_ FUNCTIONS")}
  
  
  #create names to save to
  names_raw <- names(RIA_data_in$stat_glrlm)
  
  
  stat_abr <-  which(strsplit(statistic, "")[[1]]=="(")-1
  stat_abr <-  substr(statistic, start = 1, stop = stat_abr)
  
  names_dcr <- substr(names_raw, start = 1, stop = 2)
  names_dcr <- unique(names_dcr)
  
  bins <- NULL
  names_bin <- unlist(gregexpr('_', names_raw))
  for (i in seq(1, 2*length(names_raw), 2)) {
    bins <- c(bins, substring(names_raw[(i+1)/2], names_bin[i]+1, names_bin[i+1]-1))
  }
  names_bin <- unique(bins)
  
  names_out <- NULL
  for (i in 1: length(names_dcr))
  {
    for (j in 1: length(names_bin))
    {
      names_out <- append(names_out, paste0(names_dcr[i], "_b", names_bin[j], "_", stat_abr))
    }
  }
  
  #identify stat names
  m_number <- NULL
  
  D3s <-     matrix(c( 1, 0, 0,
                       0, 1, 0,
                       1, 1, 0,
                       1,-1, 0,
                       
                       1, 0, 1,
                       0, 1, 1,
                       1, 1, 1,
                       1,-1, 1,
                       
                       1, 0,-1,
                       0, 1,-1,
                       1, 1,-1,
                       1,-1,-1,
                       
                       0, 0, 1
  ), nrow = 13, ncol = 3, byrow = TRUE)
  
  
  D2s <- matrix(c(1, 0, 0,
                  0, 1, 0,
                  1, 1, 0,
                  1,-1, 0
  ), nrow = 4, ncol = 3, byrow = TRUE)
  
  
  
  names_in <- list()
  for (i in 1: length(names_dcr))
  {
    for (j in 1: length(names_bin))
    {
      names_in_each <- NULL
      if(length(dim(RIA_data_in$data$orig)) == 2) {offsets <- D2s
      } else if (length(dim(RIA_data_in$data$orig)) == 3) {offsets <- D3s
      } else {stop(paste0("DATA LOADED IS ", length(dim(RIA_data_in$data$orig)), " DIMENSIONAL. ONLY 2D AND 3D DATA ARE SUPPORTED!"))}
      
      for (l in 1: dim(offsets)[1])
      {
        names_in_each <- c(names_in_each, paste0(names_dcr[i], "_", names_bin[j], "_", offsets[l,1], offsets[l,2], offsets[l,3]))
      }
      names_in <- c(names_in, list(names_in_each))
      
    }
  }
  
  
  #statistical names
  stat_names <- names(RIA_data_in$stat_glrlm[[names_in[[1]][1]]])
  stat_abr_plus <- paste0("stat_glrlm_", stat_abr)
  stat_abr_str <- gsub("X", "stat_data",statistic)
  
  
  
  
  for (i in 1: length(names_out))
  {
    RIA_data_in[[stat_abr_plus]][[names_out[i]]] <- list()
    
    for (j in 1: length(stat_names))
    {
      
      stat_data <- NULL
      for (k in 1: length(names_in[[i]])){
        
        stat_data <- append(stat_data, eval(parse(text = paste0("RIA_data_in$stat_glrlm$`", names_in[[i]][k], "`$", stat_names[j]))))
      }
      
      func_stat_data <- eval(parse(text = stat_abr_str))
      
      RIA_data_in[[stat_abr_plus]][[names_out[i]]][[stat_names[j]]] <- func_stat_data
      
    }
    
    if (verbose_in == TRUE) {message(paste0("AGGREGATED STATISTICS WAS ADDED TO '", names_out[i], "' SLOT OF RIA_image$", stat_abr_plus, "\n"))}
    
  }
  
  return (RIA_data_in)
  
}
