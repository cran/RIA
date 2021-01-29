#' @title Calculates all radiomic statistics on supplied RIA_image
#' @export
#'
#' @description Calculates specified radiomic statistics on \emph{RIA_image}. Parameters of
#' radiomic functions may be set. By default the the images are discretized to 2, 8, 32 and 128 bins
#' using equally sized binning. First-order statistics are calculated on the
#' original image and if asked then on all discretizations. Symmetric GLCMs are calculated for all directions
#' at a distance of 1, 2 and 3 for all discretizations. GLRLMs are also calculated for all
#' discretizations. Geometry-based statistics are calculated for the original image as well as all
#' discretizations is requested. 
#' 
#' @param RIA_data_in \emph{RIA_image}.
#' 
#' @param bins_in integer vector, number of bins specified.
#' 
#' @param equal_prob logical or string, indicating to cut data into bins with equal relative frequencies.
#' If FALSE, then equal interval bins will be used. If \emph{"both"} is supplied, the both equally probable
#' and equal interval bins will be created.
#' 
#' @param fo_discretized logical, indicating whether to calculate first-order statistics on 
#' discretized images.
#' 
#' @param distance integer, distance between the voxels being compared.
#' 
#' @param statistic string, defining the statistic to be calculated on the array of GLCM statistics.
#' By default, statistic is set to \emph{"mean"}, however any function may be provided. The proper
#' syntax is: function(X, attributes). The supplied string must contain a "X", which will be replaced
#' with the array of the GLCM statistics value. Further attributes of the function may also be given.
#' For example, if you wish to calculate the median of all GLCMs calculated in different directions,
#' then it must be supplied as: \emph{median(X, na.rm = TRUE)}.
#' 
#' @param geometry_discretized logical, indicating whether to calculate geometry-based statistics on 
#' discretized images.
#'
#' @param verbose_in logical, indicating whether to print detailed information.
#' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
#'
#' @return \emph{RIA_image} containing the statistical information.
#'
#' @examples \dontrun{
#' #Discretize loaded image and then calculate all radiomic statistics
#' DICOM <- radiomics_all(DICOM, equal_prob = "both", bins_in= c(32,64), distance = c(1:2))
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

radiomics_all <- function(RIA_data_in, bins_in=c(8, 16, 32), equal_prob = "both",
                          fo_discretized = FALSE,
                          distance = c(1), statistic = "mean(X, na.rm = TRUE)",
                          geometry_discretized = TRUE,
                          verbose_in = TRUE) {
  
  if(equal_prob == "both") {
    RIA_data_in <- discretize(RIA_data_in, bins_in=bins_in, equal_prob = TRUE, verbose_in = verbose_in)
    RIA_data_in <- discretize(RIA_data_in, bins_in=bins_in, equal_prob = FALSE, verbose_in = verbose_in)
  } else {
    RIA_data_in <- discretize(RIA_data_in, bins_in=bins_in, equal_prob = equal_prob, verbose_in = verbose_in)
  }
  
  RIA_data_in <- first_order(RIA_data_in, use_type = "single", use_orig = TRUE, verbose_in = verbose_in)
  if(fo_discretized) {RIA_data_in <- first_order(RIA_data_in, use_type = "discretized", verbose_in = verbose_in)}
  
  for (i in 1: length(distance)) {
    RIA_data_in <- glcm_all(RIA_data_in, use_type = "discretized", distance = distance[i], verbose_in = verbose_in)
  }
  RIA_data_in <- glcm_stat(RIA_data_in, use_type = "glcm", verbose_in = verbose_in)
  RIA_data_in <- glcm_stat_all(RIA_data_in, statistic = statistic, verbose_in = verbose_in)
  
  RIA_data_in <- glrlm_all(RIA_data_in, use_type = "discretized", verbose_in = verbose_in)
  RIA_data_in <- glrlm_stat(RIA_data_in, use_type = "glrlm", verbose_in = verbose_in)
  RIA_data_in <- glrlm_stat_all(RIA_data_in, statistic = statistic, verbose_in = verbose_in)
  
  RIA_data_in <- geometry(RIA_data_in, use_orig = TRUE, calc_sub = FALSE, verbose_in = verbose_in)
  if(geometry_discretized) {RIA_data_in <- geometry(RIA_data_in, use_type = "discretized", calc_sub = TRUE, verbose_in = verbose_in)}
  RIA_data_in
}