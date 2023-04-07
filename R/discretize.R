# #' @title Discretizes RIA image to a given number of bins
# #'
# #' @description  Discretizes \emph{RIA_image} into \emph{bins_in} number of bins.
# #' The \emph{equal_prob} parameter is used to indicate whether to create bins containing
# #' the same number of values. If FALSE then equal sized bins will be created.
# #' Discretized images will be saved into the \emph{$data$modif} slot of \emph{RIA_image}
# #' as well as the \emph{discretized} slot of \emph{RIA_image}.
# #' The name will be automatically created based on the type of discretization
# #' (ep: equal probability; es: equal size) and the number of bins specified,
# #' for example: \emph{$dicotomized$es_8} will store the discretized image after
# #' equal sized discretization into 8 bins. This way many different discretized images using
# #' different bin numbers can be saved to the same object for further analysis.
# #' The \emph{RIA_log} will be updated with cut points.
# #'
# #' @param RIA_data_in \emph{RIA_image}.
# #'
# #' @param bins_in integer vector, number of bins specified.
# #'
# #' @param equal_prob logical, indicating to cut data into bins with equal relative frequencies.
# #' If FALSE, then equal interval bins will be used.
# #'
# #' @param use_orig logical, indicating to use image present in \emph{RIA_data$orig}.
# #' If FALSE, the modified image will be used stored in \emph{RIA_data$modif}.
# #'
# #' @param write_orig logical, indicating to write cropped image  to \emph{RIA_data$orig}.
# #' If FALSE, the modified image will be used stored in \emph{RIA_data$modif}.
# #'
# #' @param verbose_in logical, indicating whether to print detailed information.
# #' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
# #'
# #' @return \emph{RIA_image} with values discretized to bin values.
# #'
# #' @examples \dontrun{
# #' #Discretize into 8 bins, each containing equal number of elements
# #' RIA_image <- discretize(RIA_image, bins_in = 8, equal_prob = TRUE,
# #'  use_orig = TRUE, write_orig = FALSE)
# #'
# #' #Discretize into 6 bins, each with the same width
# #' RIA_image <- discretize(RIA_image, bins_in = 6, equal_prob = FALSE,
# #'  use_orig = TRUE, write_orig = FALSE)
# #'
# #' #Discretize into 2,4,8,16,32 bins, each containing equal number of elements
# #' RIA_image <- discretize(RIA_image, bins_in = 2^(1:5), equal_prob = FALSE,
# #'  use_orig = TRUE, write_orig = FALSE)
# #'  
# #' #D
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

discretize <- function(RIA_data_in, bins_in=8, equal_prob = FALSE, use_orig = TRUE, write_orig = FALSE, verbose_in = TRUE)
{
  data_in <- check_data_in(RIA_data_in, use_type = "single", use_orig = use_orig, verbose_in = verbose_in)
  
  if(equal_prob) {equal_txt <- "EQUAL PROBABILITY"
  } else {equal_txt <- "EQUALLY SIZED"}
  
  data_NA <- as.vector(data_in)
  data_NA <- data_NA[!is.na(data_NA)]
  
  if(length(data_NA) == 0) {stop("WARNING: SUPPLIED RIA_image DOES NOT CONTAIN ANY DATA!!!")}
  
  data_out <- NULL
  
  for(j in 1:length(bins_in))
  {
    data_in_mod <- data_in
    if(!equal_prob) dichot <- cut(data_NA, bins_in[j], include.lowest  = TRUE, dig.lab = 4)
    if(equal_prob)  dichot <- stats::quantile(data_NA, seq(0, 1, 1/bins_in[j]))
    
    log <- NULL
    shift <- -9999
    
    for (i in 1 : bins_in[j])
    {
      if(!equal_prob)
      {
        loc_comma <- regexpr(',', base::levels(dichot)[i])
        loc_end   <- regexpr(']', base::levels(dichot)[i])
        lower <- as.numeric(substr(base::levels(dichot)[i], 2, loc_comma-1))
        upper <- as.numeric(substr(base::levels(dichot)[i], loc_comma+1, loc_end-1))
        if(i==bins_in[j]) {data_in_mod[(data_in_mod >= lower & data_in_mod <= upper)] <- shift + i
        } else {data_in_mod[(data_in_mod >= lower & data_in_mod < upper)] <- shift + i}
      }
      
      if(equal_prob)
      {
        lower <- dichot[i]
        upper <- dichot[i+1]
        if(i==bins_in[j]) {data_in_mod[(data_in_mod >= lower & data_in_mod <= upper)] <- shift + i
        } else {data_in_mod[(data_in_mod >= lower & data_in_mod < upper)] <- shift + i}
      }
      
      if (i == 1) log <- as.numeric(upper)
      if(i < bins_in[j] & i != 1) log <- c(log, as.numeric(upper))
      
    }
    
    data_in_mod <- data_in_mod + abs(shift)
    
    if(any(class(RIA_data_in) == "RIA_image") )
    {
      if(write_orig) {RIA_data_in$data$orig <- data_in_mod
      } else {RIA_data_in$data$modif<- data_in_mod}
      
      
      if(!equal_prob) dichot_txt <- paste0("es_", bins_in[j])
      if(equal_prob)  dichot_txt <- paste0("ep_", bins_in[j])
      RIA_data_in$discretized[[dichot_txt]] <- data_in_mod
      
      
      if(!equal_prob) cuts_txt <- paste0("cuts_es_", bins_in[j])
      if(equal_prob)  cuts_txt <- paste0("cuts_ep_", bins_in[j])
      RIA_data_in$log[[cuts_txt]] <- log
      
      if(!equal_prob) RIA_data_in$log$events <- append(RIA_data_in$log$events, paste0("Discretized_equal_sized_", bins_in[j]))
      if(equal_prob)  RIA_data_in$log$events <- append(RIA_data_in$log$events, paste0("Discretized_equal_prob_", bins_in[j]))
      
    } else {
      if(!equal_prob) dichot_txt <- paste0("es_", bins_in[j])
      if(equal_prob)  dichot_txt <- paste0("ep_", bins_in[j])
      
      data_out[[dichot_txt]] <- data_in_mod
    }
    
  }
  
  if(verbose_in) {message(paste0("SUCCESSFULLY DISCRETIZED DATA INTO ", bins_in, " NUMBER OF ", equal_txt, " BINS\n"))}
  
  if(any(class(RIA_data_in) == "RIA_image")) return(RIA_data_in)
  else return(data_out)
  
  
}
