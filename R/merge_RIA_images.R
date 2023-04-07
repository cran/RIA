#' @title Merges multiple loaded images into one volume
#' @export
#'
#' @description  Merges multiple \emph{RIA_image} class objects loaded using any of the load functions.
#' All images need to have the same dimensions. Further, during loading the images should not be cropped
#' to assure that the orientation and position of the data is maintained. Data of the new combined image is updated sequentially,
#' using data from the \emph{data$orig} slot, that is only parts of the image that do not have data
#' (which are converted to NA during the load process) are updated in the order of provided \emph{RIA_images}. If multiple images
#' contain data in for the same element, the first value is used in the new image. Data in the \emph{data$log} slot
#' is updated based on the new combined image, while data in the \emph{data$header} slot is copied from the first provided image.
#'
#' @param RIA_data_in List of Multiple \emph{RIA_images}.
#' @param crop_in logical, indicating whether to crop the merged image to smallest bounding box.
#' @param verbose_in logical indicating whether to print detailed information.
#' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
#'
#' @return \emph{RIA_image} containing the merged volume with updated log and header data
#'
#' @examples \dontrun{
#' #Load multiple images and combine them
#' d1 <- load_nifti(ABC_p1.nii.gz, crop_in = FALSE)
#' d2 <- load_nifti(ABC_p2.nii.gz, crop_in = FALSE)
#' d  <- merge_RIA(list(d1, d2))
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

merge_RIA_images <- function(RIA_data_in, crop_in = TRUE, verbose_in = TRUE) {
  all_RIA_image <- lapply(RIA_data_in, class)
  if(!all(unlist(lapply(all_RIA_image, function(x) {"RIA_image" %in% x})))) {message("PROCESSING OF RIA_image OBJECTS ARE SUPPORTED, OTHER CLASSES MIGHT CAUSE PROBLEMS! PLEASE LOAD DATA USING RIA load_ FUNCTIONS")}
  all_dim_equal <- lapply(RIA_data_in, function(x) {dim(x$data$orig)})
  if(!all(unlist(lapply(1:length(all_dim_equal), function(x){all(all_dim_equal[[1]] == all_dim_equal[[x]])})))) {
    message("NOT ALL PROVIDED IMAGES HAVE THE SAME DIMENSIONS! IMAGES CAN ONLY BE MERGED IF ALL PROVIDED IMAGE SIZES ARE THE SAME!")
  }
  if(length(RIA_data_in) <2) {message("PLEASE PROVIDE AT LEAST 2 IMAGES TO MERGE")}
  
  #Create base image
  RIA_image <- list(data = NULL, header = list(), log = list())
  class(RIA_image) <- append(class(RIA_image), "RIA_image")
  RIA_image$data$orig  <- RIA_data_in[[1]]$data$orig
  RIA_image$data$modif <- NULL
  RIA_image$header <- RIA_data_in[[1]]$header
  RIA_image$log <- RIA_data_in[[1]]$log
  class(RIA_image$header) <- append(class(RIA_image$header), "RIA_header")
  class(RIA_image$data) <- append(class(RIA_image$data), "RIA_data")
  class(RIA_image$log) <- append(class(RIA_image$log), "RIA_log")
  RIA_image$log$orig_vol_mm <- NA
  RIA_image$log$orig_surf_mm <- NA
  RIA_image$log$surface_volume_r <- NA
  
  #Combine images
  for(i in 2:length(RIA_data_in)) {
    RIA_image$data$orig[is.na(RIA_image$data$orig)] <- RIA_data_in[[i]]$data$orig[is.na(RIA_image$data$orig)]
  }
  
  #Crop data
  if(crop_in)
  {
    RIA_image <- crop(RIA_image, zero_value = NA, write_orig = TRUE, verbose_in = verbose_in)
  }
  
  #Add volumetric data
  RIA_image$log$orig_vol_mm <- volume(RIA_image$data$orig, xy_dim = RIA_image$log$orig_xy_dim, z_dim = RIA_image$log$orig_z_dim)
  RIA_image$log$orig_surf_mm <- surface(RIA_image$data$orig, xy_dim = RIA_image$log$orig_xy_dim, z_dim = RIA_image$log$orig_z_dim)
  RIA_image$log$surface_volume_r <- ifelse(RIA_image$log$orig_vol_mm != 0, RIA_image$log$orig_surf_mm/RIA_image$log$orig_vol_mm, 0)
  
  if(verbose_in) {message("SUCCESSFULLY MERGED IMAGES INTO SINGLE RIA_image\n")}
  return(RIA_image)
}
