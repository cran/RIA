#' @title GLRLM-based statistics
#' @export
#'
#' @description Calculates GLRLM-based statistics for given GLRLM matrix.
#'
#' @param RIA_data_in \emph{RIA_image}.
#'
#' @param use_type string, can be \emph{"single"} which runs the function on a single image,
#' which is determined using \emph{"use_orig"} or \emph{"use_slot"}. \emph{"glrlm"}
#' takes all datasets in the \emph{RIA_image$glrlm} slot and runs the analysis on them.
#'
#' @param use_orig logical, indicating to use image present in \emph{RIA_data$orig}.
#' If FALSE, the modified image will be used stored in \emph{RIA_data$modif}. However, GLRLM matrices
#' are usually note present in either slots, therefore giving the slot name using \emph{use_slot} is advised.
#'
#' @param use_slot string, name of slot where data wished to be used is. Use if the desired image
#' is not in the \emph{data$orig} or \emph{data$modif} slot of the \emph{RIA_image}. For example,
#' ig the desired dataset is in \emph{RIA_image$glrlm$ep_4}, then \emph{use_slot} should be
#' \emph{glrlm$ep_4}. The results are automatically saved. If the results are not saved to
#' the desired slot, then please use \emph{save_name} parameter.
#'
#' @param save_name string, indicating the name of subslot of \emph{$glrlm} to save results to.
#' If left empty, then it will be automatically determined.
#'
#' @param verbose_in logical, indicating whether to print detailed information.
#' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
#'
#' @return \emph{RIA_image} containing the statistical information.
#'
#' @examples \dontrun{
#' #Discretize loaded image and then calculate GLRLM statistics
#' RIA_image <- discretize(RIA_image, bins_in = 8, equal_prob = TRUE)
#' RIA_image <- glrlm(RIA_image, use_orig = FALSE, use_slot = "discretized$ep_8",
#' right = TRUE, down = TRUE, forward = FALSE)
#' RIA_image <- glrlm_stat(RIA_image, use_orig = FALSE, use_slot = "glrlm$ep_8_110")
#' 
#' #Batch calculation of GLRLM-based statistics on all calculated GLRLMs
#' RIA_image <- glrlm_stat(RIA_image, use_type = "discretized")
#' }
#' 
#' @references Márton KOLOSSVÁRY et al.
#' Radiomic Features Are Superior to Conventional Quantitative Computed Tomographic
#' Metrics to Identify Coronary Plaques With Napkin-Ring Sign
#' Circulation: Cardiovascular Imaging (2017).
#' DOI: 10.1161/circimaging.117.006843
#' \url{https://pubmed.ncbi.nlm.nih.gov/29233836/}
#' 
#' Márton KOLOSSVÁRY et al.
#' Cardiac Computed Tomography Radiomics: A Comprehensive Review on Radiomic Techniques.
#' Journal of Thoracic Imaging (2018).
#' DOI: 10.1097/RTI.0000000000000268
#' \url{https://pubmed.ncbi.nlm.nih.gov/28346329/}
#' @encoding UTF-8



glrlm_stat <- function(RIA_data_in, use_type = "single", use_orig = FALSE, use_slot = "glrlm$es_8_111", save_name = NULL, verbose_in = TRUE)
{
  data_in_orig <- check_data_in(RIA_data_in, use_type = use_type, use_slot = use_slot, verbose_in = verbose_in)
  
  if(any(class(data_in_orig) != "list")) data_in_orig <- list(data_in_orig)
  
  list_names <- names(data_in_orig)
  if(!is.null(save_name) & (length(data_in_orig) != length(save_name))) {stop(paste0("PLEASE PROVIDE THE SAME NUMBER OF NAMES AS THERE ARE IMAGES!\n",
                                                                                     "NUMBER OF NAMES:  ", length(save_name), "\n",
                                                                                     "NUMBER OF IMAGES: ", length(data_in), "\n"))
  }
  
  for (i in 1: length(data_in_orig))
  {
    data_in <- data_in_orig[[i]]
    data_NA <- as.vector(data_in)
    data_NA <- data_NA[!is.na(data_NA)]
    if(length(data_NA) == 0) {message("WARNING: SUPPLIED RIA_image DOES NOT CONTAIN ANY DATA!!!")}
    
    SRE            <- sre(data_in)
    LRE            <- lre(data_in)
    GLN            <- gln(data_in)
    RLN            <- rln(data_in)
    RP             <- rp(data_in)
    LGLRE           <- lglre(data_in)
    HGLRE           <- hglre(data_in)
    SRLGLE          <- srlgle(data_in)
    LRHGLE          <- lrhgle(data_in)
    SRHGLE          <- srhgle(data_in)
    LRLGLE          <- lrlgle(data_in)
    
    metrics <- list(
      SRE            <- SRE,
      LRE            <- LRE,
      GLN            <- GLN,
      RLN            <- RLN,
      RP             <- RP,
      LGLRE          <- LGLRE,
      HGLRE          <- HGLRE,
      SRLGLE         <- SRLGLE,
      LRHGLE         <- LRHGLE,
      SRHGLE         <- SRHGLE,
      LRLGLE         <- LRLGLE
    )
    
    stat_names <- c(    "SRE",
                        "LRE",
                        "GLN",
                        "RLN",
                        "RP",
                        "LGLRE",
                        "HGLRE",
                        "SRLGLE",
                        "LRHGLE",
                        "SRHGLE",
                        "LRLGLE"
    )
    
    names(metrics) <- stat_names
    
    if(use_type == "single") {
      if(any(class(RIA_data_in) == "RIA_image") )
      {
        if(is.null(save_name)) {
          txt <- automatic_name(RIA_data_in, use_orig, use_slot)
          RIA_data_in$stat_glrlm[[txt]] <- metrics
          
        }
        if(!is.null(save_name)) {RIA_data_in$stat_glrlm[[save_name]] <- metrics
        }
      }
    }
    
    if(use_type == "glrlm") {
      if(any(class(RIA_data_in) == "RIA_image"))
      {
        if(is.null(save_name[i])) {
          txt <- list_names[i]
          RIA_data_in$stat_glrlm[[txt]] <- metrics
        }
        if(!is.null(save_name[i])) {RIA_data_in$stat_glrlm[[save_name[i]]] <- metrics
        }
      }
    }
    
    if(is.null(save_name)) {txt_name <- txt
    } else {txt_name <- save_name}
    if(verbose_in) {message(paste0("GLRLM STATISTICS WAS SUCCESSFULLY ADDED TO '", txt_name, "' SLOT OF RIA_image$stat_glrlm\n"))}
    
  }
  
  if(any(class(RIA_data_in) == "RIA_image") ) return(RIA_data_in)
  else return(metrics)
}

