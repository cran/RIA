#' @title Export radiomics calculations of RIA image to csv
#' @export
#'
#' @description  Exports given slots of statistics from RIA_image. Names of slots have to be defined
#' which the user wishes to export using the \emph{stats} parameter. Using the \emph{group_name}
#' parameter the user can lable the cases with a group ID, for example "Case", which can be
#' used as a grouping variable for further analysis.
#'
#' @param RIA_image \emph{RIA_image} with calculated statistics.
#' 
#' @param save_to string, path of folder to save results to.
#' 
#' @param save_name string, path of folder to save results to.
#' 
#' @param group_name string, a ID defining which group the case belongs to.
#' 
#' @param stats string vector, identifing which slots to export
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

save_RIA <- function(RIA_image, save_to = "C:/", save_name = "RIA_stat", group_name = "Case", stats = c("stat_fo", "stat_glcm_mean", "stat_glrlm_mean", "stat_geometry")) {
  
  df_out <- data.frame(); df_out[1,1] <- group_name; colnames(df_out) <- "Group_name"
  
  #Header_info
  df_out <- cbind(df_out, t(list_to_df(RIA_image$header))); space_loc <- regexpr(' ', RIA_image$header$PixelSpacing)[1];
  
  if(is.numeric(RIA_image$header$PixelSpacing)) {
    df_out[,"PixelSpacing"] <- RIA_image$header$PixelSpacing
  } else {
    df_out[,"PixelSpacing"] <- as.numeric(substr(RIA_image$header$PixelSpacing, 1, space_loc-1))
  }
  
  for(i in 1:length(stats)) {
    
    df_out <- cbind(df_out, save_RIA_base(RIA_image, stats[i]))
  }
  
  utils::write.csv(df_out, file = paste0(save_to, save_name, ".csv"))
  
}
