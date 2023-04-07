#' @title Loads NIfTI images to RIA image format
#' @export
#' @description  Loads NIfTI images to a \emph{RIA_image} object.
#' \emph{RIA_image} is a  list with three mandatory attributes.
#' \itemize{
#'  \item \bold{RIA_data} is a \emph{RIA_data} object, which has two potential slots.
#'  \emph{$orig} contains the original image after loading
#'  \emph{$modif} contains the image that has been modified using functions.
#'  \item \bold{RIA_header} is a \emph{RIA_header} object, which is list of header information.
#'  \item \bold{RIA_log} is a \emph{RIA_log} object, which is a list updated by RIA functions
#'  and acts as a log and possible input for some functions.
#'  }
#'  Further attributes may also be added by RIA functions.
#'
#' @param filename string, file path to directory containing \emph{NIfTI} file.
#' 
#' @param image_dim integer, dimensions of the image.
#' 
#' @param mask_filename string vector, file path to optional directory containing \emph{NIfTI} file
#' of mask image. If multiple are supplied, then those voxels are kept which have one of the values of \emph{keep_mask_values}
#' in any of the supplied masks. 
#' 
#' @param keep_mask_values integer vector or string, indicates which value or values of the mask image
#' to use as indicator to identify voxels wished to be processed. Usually 1-s indicate voxels
#' wished to be processed. However, one mask image might contain several segmentations, in which
#' case supplying several integers is allowed. Furthermore, if the same string is supplied to
#' \emph{filename} and \emph{mask_filename}, then the integers in \emph{keep_mask_values} are used
#' to specify which voxel values to analyze. This way the provided image can be segmented to specific
#' components. For example, if you wish to analyze only the low-density non-calcified component
#' of coronary plaques, then \emph{keep_mask_values} can specify this by setting it to: -100:30.
#' If  a single string is provided, then each element of the mask will be examined against the statement in the string.
#' For example, if \emph{'>0.5'} is provided i.e. the mask is probabilities after a DL algorithm, then all
#' voxels with values >0.5 in the mask image will be kept. This can be a complex logical expression.
#' The data on which the expression is executed is called \emph{data} or \emph{data_mask}, depending on whether
#' you wish to filter the original image, that is the original image is supplied as a mask, or if you have
#' unique mask files respectively. Therefore for complex logical expressions you can define for example:
#' \emph{'>-100 & data<30'} to consider data values between -100 and 30, or \emph{'>0.5 & data_mask<0.75'}
#' to select voxels based-on mask values between 0.5 and 0.75 for example if they represent a probability mask.
#' 
#' @param switch_z logical, indicating whether to change the orientation of the images in the Z axis. Some
#' software reverse the order of the manipulated image in the Z axis, and therefore the images of the mask
#' image need to be reversed.
#'
#' @param crop_in logical, indicating whether to crop \emph{RIA_image} to smallest bounding box.
#'
#' @param replace_in logical, whether to replace smallest values indicated by \emph{zero_value},
#' which are considered to indicate no signal, to NA.
#'
#' @param center_in logical, whether to shift data so smallest value
#' is equal to \emph{min_to} input parameter.
#'
#' @param zero_value integer, indicating voxels values which are considered
#' not to have any information. If left empty,
#' then the smallest HU value in the image will be used, if \emph{replace_in} is TRUE.
#'
#' @param min_to integer, value to which data is shifted to if \emph{center_in} is TRUE.
#'
#' @param verbose_in logical, indicating whether to print detailed information.
#' Most prints can also be suppresed using the \code{\link{suppressMessages}} function.
#'
#' @param reorient_in \emph{reorient} parameter input of \code{\link[oro.nifti]{readNIfTI}}.
#' 
#' @param ... additional arguments to  \code{\link[oro.nifti]{readNIfTI}},
#' \code{\link[oro.nifti]{nifti_header}}.
#'
#' @details \emph{load_nifti} is used to transform NIfTI datasets into the RIA environment.
#' \emph{RIA_image} object was developed to facilitate and simplify radiomics calculations by keeping
#' all necessary information in one place.
#' \cr
#' \cr
#' \emph{RIA_data} stores the image that is converted to numerical 3D arrays using
#' \code{\link[oro.nifti]{readNIfTI}}.
#' The function stores the original loaded image in  \emph{RIA_data$orig},
#' while all modified images are stored in \emph{RIA_data$modif}.
#' By default, the original image \emph{RIA_data$orig} is untouched by functions
#' other than those operating in \emph{load_nifti}. While other functions
#' operate on the \emph{RIA_data$modif} image by default.
#' \cr
#' Due to memory concerns, there can only be one \emph{RIA_data$orig} and \emph{RIA_data$modif}
#' image present at one time in a \emph{RIA_image}. Therefore, if image manipulations are performed,
#' then the \emph{RIA_data$modif} will be overwritten. However, functions can save images
#' into new slots of \emph{RIA_image}, for example discretized images can be saved to the \emph{discretized} slot of \emph{RIA_image}.
#' \cr
#' \emph{load_nifti} not only loads the image based on parameters that can be set for
#' \code{\link[oro.nifti]{readNIfTI}}, but also can perform
#' minimal manipulations on the image itself.
#' \cr
#' \emph{crop_in} logical variable is used to indicate, whether to crop the image to the
#' smallest bounding box still containing all the information. If TRUE, then all X, Y and potentially
#' Z slices containing no information will be removed. This allows significant reduction of necessary
#' memory to store image data.
#' \cr
#' \emph{zero_value} parameter is used to indicate HU values which contain no information. If left empty,
#' then the smallest value will be considered as indicating voxels without a signal.
#' \cr
#' \emph{replace_in} logical can be used to change values that are considered to have no signal to NA.
#' This is necessary to receive proper statistical values later on.
#' \cr
#' \emph{center_in} logical is used to indicate whether the values should be shifted.
#' Some vendors save HU values as positive integers to spare memory and minimalize file sizes.
#' Therefore, in some instances shift of the scale is needed. By default,
#' the values are shifted by -1024, but in other cases a different constant might be required,
#' which can be set using the \emph{min_to} input.
#' \cr
#' \cr
#' \emph{RIA_header} is a list containing the most basic patient and examination information
#' present in the NIfTI file.
#' \cr
#' \cr
#' \emph{RIA_log} is a list of variables, which give an overview of what has been done with the image.
#' If the whole \emph{RIA_image} is supplied to a function, the information regarding the manipulations
#' are written into the \emph{$events} array in chronological order. Furthermore, some additional
#' information is also saved in the log, which might be needed for further analysis.
#'
#' @return  Returns a \emph{RIA_image} object. \emph{RIA_image} is a list with three mandatory attributes.
#' \itemize{
#'  \item \bold{RIA_data} is a \emph{RIA_data} object containing the image in \emph{$orig} slot.
#'  \item \bold{RIA_header} is a \emph{RIA_header} object, which is s list of meta information.
#'  \item \bold{RIA_log} is a \emph{RIA_log} object, which is a list updated by RIA functions
#'  and acts as a log and possible input for some functions.
#'  }
#'
#' @examples \dontrun{
#'  #Image will be croped to smallest bounding box, and smallest values will be changed to NA,
#'  while 1024 will be substracted from all other data points.
#'  RIA_image <- load_nifti("/Users/Test/Documents/Radiomics/John_Smith/NIfTI_folder/sample.nii")
#'  }
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


load_nifti <- function(filename, image_dim = 3, mask_filename = NULL, keep_mask_values = 1, switch_z = FALSE, 
                       crop_in = TRUE, replace_in = TRUE, center_in = FALSE,  zero_value = NULL, min_to = -1024,
                       verbose_in = TRUE,
                       reorient_in = TRUE, ...
)
{
  if(verbose_in) {message(paste0("LOADING NIFTI FILES FROM: ", filename, "\n"))}
  
  dcmImages <- oro.nifti::readNIfTI(filename, verbose = FALSE, reorient = reorient_in)
  
  
  ###create 3D matrix - crop to smallest bounding box - change 0/-1024 values to NA - center around 0
  dim_string <- paste0( c(rep(",", image_dim-1), rep(",1", (dcmImages@dim_)[1]-image_dim)), collapse="")
  data  <- dcmImages@.Data
  data <- eval(parse(text= paste0("data[", dim_string, "]")))
  
  ###create RIA_image structure
  RIA_image <- list(data = NULL, header = list(), log = list())
  if(length(dim(data)) == 3 | length(dim(data)) == 2) {class(RIA_image) <- append(class(RIA_image), "RIA_image")
  } else {stop(paste0("NIFTI LOADED IS ", length(dim(data)), " DIMENSIONAL. ONLY 2D AND 3D DATA ARE SUPPORTED!"))}
  
  
  if(is.null(zero_value)) zero_value <- min(data, na.rm = TRUE)
  
  #mask image
  if(!is.null(mask_filename)) {
    if(identical(filename, mask_filename)) {
      if(verbose_in) {message(paste0("CANCELING OUT VALUES OTHER THAN THOSE SPECIFIED IN 'keep_mask_values' PARAMETER \n"))}
      
      if(suppressWarnings(any(is.na(as.numeric(keep_mask_values))))) { #not just numeric values
        data[!eval(parse(text = paste0("data", keep_mask_values)))] <- zero_value
      } else{
        data[!data %in% keep_mask_values] <- zero_value
      }
      
    } else {
      for(i in 1:length(mask_filename)) {
        mask_filename_i <- mask_filename[i]
        
        if(verbose_in) {message(paste0("LOADING NIFTI IMAGES OF MASK IMAGE FROM: ", mask_filename, "\n"))}
        dcmImages_mask <- oro.nifti::readNIfTI(mask_filename_i, verbose = FALSE, reorient = reorient_in)
        data_mask  <- dcmImages_mask@.Data
        data_mask <- eval(parse(text= paste0("data_mask[", dim_string, "]")))
        
        if(!all(dim(data) == dim(data_mask))) {
          stop(paste0("DIMENSIONS OF THE IMAGE AND THE MASK ARE NOT EQUAL!\n",
                      "DIMENSION OF IMAGE: ", dim(data)[1], " ",  dim(data)[2], " ", dim(data)[3], "\n",
                      "DIMENSION OF MASK:  ", dim(data_mask)[1], " ", dim(data_mask)[2], " ", dim(data_mask)[3], "\n"))
        } else {
          if(switch_z) {data_mask[,,dim(data_mask)[3]:1] <- data_mask
          message("MASK IMAGE WAS TRANSFORMED TO ACHIEVE PROPER ORIENTATION OF THE ORIGINAL AND THE MASK IMAGE.\n")
          }
          
          if(suppressWarnings(any(is.na(as.numeric(keep_mask_values))))) { #not just numeric values
            data[!eval(parse(text = paste0("data_mask", keep_mask_values)))] <- zero_value
          } else{
            data_mask[!(data_mask %in% keep_mask_values)] <- NA
            if(i == 1) {
              data_mask_all <- data_mask
            }else {
              data_mask_all[is.na(data_mask_all)] <- data_mask[is.na(data_mask_all)]
            }
          }
        }
      }
      if(suppressWarnings(!any(is.na(as.numeric(keep_mask_values))))) { #If only numerical values in masks
        data[!(data_mask_all %in% keep_mask_values)] <- zero_value
      }
    }
  }
  
  
  RIA_image$data$orig  <- data
  RIA_image$data$modif <- NULL
  class(RIA_image$header) <- append(class(RIA_image$header), "RIA_header")
  class(RIA_image$data) <- append(class(RIA_image$data), "RIA_data")
  class(RIA_image$log) <- append(class(RIA_image$log), "RIA_log")
  RIA_image$log$events  <- "Created"
  RIA_image$log$orig_dim  <- dim(data)
  RIA_image$log$directory  <- filename
  
  
  if(!is.null(mask_filename)) {
    if(identical(filename, mask_filename)) {
      RIA_image$log$events  <- paste0("Filtered_using_values_", paste0(keep_mask_values, collapse = "_"))
    } else {
      RIA_image$log$events  <- paste0("Filtered_using_mask_values_", paste0(keep_mask_values, collapse = "_"))
    }
  }
  
  
  
  ###Crop data
  if(crop_in)
  {
    if(verbose_in) {message(paste0("SMALLEST VALUES IS ", zero_value, ", AND WILL BE CONSIDERED AS REFERENCE POINT TO IDENTIFY VOXELS WITHOUT ANY SIGNAL\n"))}
    if(verbose_in & center_in == FALSE) message(paste0("MIGHT CONSIDER RESCALING, SINCE SMALLEST VALUE IS NOT -1024, AND THUS VOXEL VALUES MIGHT NOT BE CORRECT\n"))
    
    RIA_image <- crop(RIA_image, zero_value, write_orig = TRUE, verbose_in = verbose_in)
  }
  
  
  ###Replace values
  if(replace_in)
  {
    if(verbose_in) {message(paste0("SMALLEST VALUES IS ", zero_value, ", AND WILL CHANGE TO NA\n"))}
    
    RIA_image <- change_to(RIA_image, zero_value_in = zero_value, verbose_in = verbose_in)
  }
  
  ###Shift to
  if(center_in & (min(data, na.rm = T) != min_to))
  {
    if(verbose_in) {message(paste0("SMALLEST VALUES IS not ", min_to, " THEREFORE SHIFTING VALUES TO ACHIVE THIS\n"))}
    RIA_image <- shift_to(RIA_image, to = min_to, min_value_in = zero_value, verbose_in = verbose_in)
  }
  
  ###Create dataframe of standard or specific NIFTI information
  header <- create_header_nifti(filename)
  RIA_image$header <- header
  #Add original volume of abnormality
  xy_dim <- as.numeric(RIA_image$header$PixelSpacing)
  z_dim <-  as.numeric(RIA_image$header$SpacingBetweenSlices)
  RIA_image$log$orig_vol_mm <- volume(RIA_image$data$orig, xy_dim = xy_dim, z_dim = z_dim)
  RIA_image$log$orig_surf_mm <- surface(RIA_image$data$orig, xy_dim = xy_dim, z_dim = z_dim)
  RIA_image$log$surface_volume_r <- ifelse(RIA_image$log$orig_vol_mm != 0, RIA_image$log$orig_surf_mm/RIA_image$log$orig_vol_mm, 0)
  RIA_image$log$orig_xy_dim <- xy_dim
  RIA_image$log$orig_z_dim  <- z_dim
  
  
  
  if(verbose_in) {message(paste0("SUCCESSFULLY LOADED ", RIA_image$header$PatientsName, "'s NIFTI IMAGES TO RIA IMAGE CLASS\n"))}
  data_NA <- as.vector(RIA_image$data$orig)
  data_NA <- data_NA[!is.na(data_NA)]
  
  if(length(data_NA) == 0) {message("WARNING: RIA_image$data DOES NOT CONTAIN ANY DATA!!!\n")}
  
  return(RIA_image)
}
