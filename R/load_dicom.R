#' @title Loads DICOM images to RIA image format
#' @export
#' @description  Loads DICOM images to a \emph{RIA_image} object.
#' \emph{RIA_image} is a  list with three mandatory attributes.
#' \itemize{
#'  \item \bold{RIA_data} is a \emph{RIA_data} object, which has two potential slots.
#'  \emph{$orig} contains the original image after loading and is a 3D array of integers
#'  created with \code{\link[oro.dicom]{create3D}}.
#'  \emph{$modif} contains the image that has been modified using functions.
#'  \item \bold{RIA_header} is a \emph{RIA_header} object, which is list of DICOM header information.
#'  \item \bold{RIA_log} is a \emph{RIA_log} object, which is a list updated by RIA functions
#'  and acts as a log and possible input for some functions.
#'  }
#'  Further attributes may also be added by RIA functions.
#'
#' @param filename string, file path to directory containing \emph{dcm} files.
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
#' @param header_add dataframe, with three columns: Name, Group and Element containing the name,
#' the group and the element code of the DICOM fields wished to be added to the\emph{RIA_header}.
#'
#' @param header_exclude dataframe, with three columns: Name, Group and Element containing the name,
#' the group and the element code of the DICOM fields wished to be excluded
#' from the default header elements present in \emph{DICOM_codes} rda file.
#'
#' @param verbose_in logical, indicating whether to print detailed information.
#' Most prints can also be suppresed using the \code{\link{suppressMessages}} function.
#'
#' @param recursive_in \emph{recursive} parameter input of \code{\link[oro.dicom]{readDICOM}}.
#'
#' @param exclude_in \emph{exclude} parameter input of \code{\link[oro.dicom]{readDICOM}}.
#'
#' @param mode_in \emph{mode} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @param transpose_in \emph{transpose} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @param pixelData_in \emph{pixelData} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @param mosaic_in \emph{mosaic} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @param mosaicXY_in \emph{mosaicXY} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @param sequence_in \emph{sequence} parameter input of \code{\link[oro.dicom]{create3D}}.
#'
#' @details \emph{load_dicom} is used to transform DICOM datasets into the RIA environment.
#' \emph{RIA_image} object was developed to facilitate and simplify radiomics calculations by keeping
#' all necessary information in one place.
#' \cr
#' \cr
#' \emph{RIA_data} stores the DICOM image that is converted to numerical 3D arrays using
#' \code{\link[oro.dicom]{readDICOM}} and \code{\link[oro.dicom]{create3D}}.
#' The function stores the original loaded image in  \emph{RIA_data$orig},
#' while all modified images are stored in \emph{RIA_data$modif}.
#' By default, the original image \emph{RIA_data$orig} is untouched by functions
#' other than those operating in \emph{load_dicom}. While other functions
#' operate on the \emph{RIA_data$modif} image by default.
#' \cr
#' Due to memory concerns, there can only be one \emph{RIA_data$orig} and \emph{RIA_data$modif}
#' image present at one time in a \emph{RIA_image}. Therefore, if image manipulations are performed,
#' then the \emph{RIA_data$modif} will be overwritten. However, functions can save images
#' into new slots of \emph{RIA_image}, for example the \code{\link[RIA]{discretize}} function can save
#' discretized images to the \emph{discretized} slot of \emph{RIA_image}.
#' \cr
#' \emph{load_dicom} not only loads the DICOM image based on parameters that can be set for
#' \code{\link[oro.dicom]{readDICOM}} and \code{\link[oro.dicom]{create3D}}, but also can perform
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
#' needed for further analysis. The default DICOM set is present in \emph{DICOM_codes},
#' which can be edited to anyones needs. But if we wish only to add of remove specific
#' DICOM header rows, then the \emph{header_add} and \emph{header_exclude} can be used.
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
#'  \item \bold{RIA_header} is a \emph{RIA_header} object, which is s list of DICOM information.
#'  \item \bold{RIA_log} is a \emph{RIA_log} object, which is a list updated by RIA functions
#'  and acts as a log and possible input for some functions.
#'  }
#'
#' @examples \dontrun{
#'  #Image will be croped to smallest bounding box, and smallest values will be changed to NA,
#'  while 1024 will be substracted from all other data points.
#'  RIA_image <- load_dicom("C://Users//Test//Documents//Radiomics//John_Smith//DICOM_folder//")
#'  }


load_dicom <- function(filename, crop_in = TRUE, replace_in = TRUE, center_in = TRUE,  zero_value = NULL, min_to = -1024,
                       header_add = NULL, header_exclude = NULL, verbose_in = TRUE,
                       recursive_in = TRUE, exclude_in = "sql",
                       mode_in = "integer", transpose_in = TRUE, pixelData_in = TRUE,
                       mosaic_in = FALSE, mosaicXY_in = NULL, sequence_in = FALSE
                       )
{
  if(verbose_in) {message(paste0("LOADING DICOM FILES FROM: ", filename, "\n"))}

  dcmImages <- oro.dicom::readDICOM(filename, recursive = recursive_in, exclude = exclude_in, verbose = verbose_in)


  ###create 3D matrix - crop to smallest bounding box - change 0/-1024 values to NA - center around 0
  if(length(dcmImages$img)==1) {
  data  <- suppressWarnings(oro.dicom::create3D(dcmImages, mode = mode_in, transpose = transpose_in, pixelData = pixelData_in,
                        mosaic = mosaic_in, mosaicXY = mosaicXY_in, sequence = sequence_in))
  } else {
  data  <- oro.dicom::create3D(dcmImages, mode = mode_in, transpose = transpose_in, pixelData = pixelData_in,
                                 mosaic = mosaic_in, mosaicXY = mosaicXY_in, sequence = sequence_in)
  }

  ###create RIA_image structure
  RIA_image <- list(data = NULL, header = list(), log = list())
  if(length(dim(data)) == 3 | length(dim(data)) == 2) {class(RIA_image) <- append(class(RIA_image), "RIA_image")
  } else {stop(paste0("DICOM LOADED IS ", length(dim(data)), " DIMENSIONAL. ONLY 2D AND 3D DATA ARE SUPPORTED!"))}

  RIA_image$data$orig  <- data
  RIA_image$data$modif <- NULL
  class(RIA_image$header) <- append(class(RIA_image$header), "RIA_header")
  class(RIA_image$data) <- append(class(RIA_image$data), "RIA_data")
  class(RIA_image$log) <- append(class(RIA_image$log), "RIA_log")
  RIA_image$log$events  <- "Created"
  RIA_image$log$orig_dim  <- dim(data)
  RIA_image$log$directory  <- filename


  ###Crop data
  if(is.null(zero_value)) zero_value <- min(RIA_image$data$orig, na.rm = TRUE)

  if(crop_in)
  {
    if(verbose_in) {message(paste0("SMALLEST VALUES IS ", zero_value, ", AND WILL BE CONSIDERED AS REFERENCE POINT TO IDENTIFY VOXELS WITHOUT ANY SIGNAL\n"))}
    if(verbose_in & center_in == FALSE) message(paste0("MIGHT CONSIDER RESCALING, SINCE SMALLEST VALUE IS NOT -1024, AND THUS HU VALUES MIGHT NOT BE CORRECT\n"))

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

  ###Create dataframe of standard or specific DICOM information
  header <- create_header(filename, header_add, header_exclude)
  RIA_image$header <- header
  #Add original volume of abnormality
  space_loc <- regexpr(' ', RIA_image$header$PixelSpacing)[1]
  xy_dim <- as.numeric(substr( RIA_image$header$PixelSpacing, 1, space_loc-1))
  z_dim <-  as.numeric(RIA_image$header$SpacingBetweenSlices)
  RIA_image$log$orig_vol_mm <- volume(RIA_image$data$orig, xy_dim = xy_dim, z_dim = z_dim)
  RIA_image$log$orig_surf_mm <- surface(RIA_image$data$orig, xy_dim = xy_dim, z_dim = z_dim)
  RIA_image$log$surface_volume_r <- ifelse(RIA_image$log$orig_vol_mm != 0, RIA_image$log$orig_surf_mm/RIA_image$log$orig_vol_mm, 0)
  RIA_image$log$orig_xy_dim <- xy_dim
  RIA_image$log$orig_z_dim  <- z_dim



  if(verbose_in) {message(paste0("SUCCESSFULLY LOADED ", RIA_image$header$PatientsName, "'s DICOM IMAGES TO RIA IMAGE CLASS\n"))}
  data_NA <- as.vector(RIA_image$data$orig)
  data_NA <- data_NA[!is.na(data_NA)]

  if(length(data_NA) == 0) {message("WARNING: RIA_image$data DOES NOT CONTAIN ANY DATA!!!\n")}

  return(RIA_image)
}
