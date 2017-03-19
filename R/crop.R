# Crops input image to smallest boundung box containing values
# Crops RIA_image of slices containing no information. RIA_log will be updated with
# logical arrays which indicate which slices in the given direction were kept after cropping.
#
# crops the image to the smallest bounding box still containing all the information.
# All X, Y and potentially Z slices containing no information will be removed.
# This allows significant reduction of necessery memory to store image data.
# The \emph{zero_value} parameter is used to indicate HU values which contain no information.
#
# RIA_data_in RIA_image
#
# zero_value_in the value that is considered not to contain any data
#
# use_orig indicating to use image present in RIA_data$orig.
# If FALSE, the modified image will be used stored in RIA_data$modif.
#
# write_orig indicating to write cropped image  to RIA_data$orig.
# If FALSE, the modified image will be used stored in RIA_data$modif.
#
# verbose_in logical indicating whether to print detailed information.
# Most prints can also be suppresed using the suppressMessages function
#
# return RIA_image object cropped to smallest bouding box, whith RIA_log updated.
#
# examples
# RIA_image$data$orig will be replaced by the modified image containg smaller dimentions,
# since slices in any direction containing no information will be cropped
# RIA_image <- crop(RIA_image)
#
# RIA_image$data$orig will be kept and the modified image containg smaller dimentions
# will be stored in RIA_image$data$modif
# RIA_image <- crop(RIA_image, use_orig = TRUE, write_orig = FALSE)
#
# Crop 3D matrix
# matrix <- array(rnorm(1000), dim = c(10, 10, 10))
# base   <- array(0, dim = c(10+2, 10+2, 10+2))
# base[2:11, 2:11, 2:11] <- matrix
# crop(base)



crop <- function(RIA_data_in, zero_value = 0, use_orig = TRUE, write_orig = TRUE, verbose_in = TRUE)
{
  data_in <- check_data_in(RIA_data_in, use_type = "single", use_orig = use_orig, verbose_in = verbose_in)

  if(verbose_in) message("CALCULATING IN DIRECTION X")
  tpb <- utils::txtProgressBar(min = 0, max = dim(data_in)[1], style = 3)

  logic_x = rep(0, dim(data_in)[1])

  for (i in 1:dim(data_in)[1])
  {
    flat = data_in[i, 1:dim(data_in)[2], 1:dim(data_in)[3]]
    ref = matrix(zero_value, dim(data_in)[2], dim(data_in)[3])

    if(max(abs(flat- ref)) != 0) {logic_x[i] = 1}

    if(verbose_in) utils::setTxtProgressBar(tpb, i); if(i == dim(data_in)[1]) {close(tpb)}
  }


  if(verbose_in) message("CALCULATING IN DIRECTION Y")
  tpb <- utils::txtProgressBar(min = 0, max = dim(data_in)[2], style = 3)

  logic_y = rep(0, dim(data_in)[2])

  for (i in 1:dim(data_in)[2])
  {
    flat = data_in[1:dim(data_in)[1], i, 1:dim(data_in)[3]]
    ref = matrix(zero_value, dim(data_in)[1], dim(data_in)[3])

    if(max(abs(flat- ref)) != 0) {logic_y[i] = 1}

    if(verbose_in) utils::setTxtProgressBar(tpb, i); if(i == dim(data_in)[1]) {close(tpb)}
  }

  X_1 = Inf; X_2 = Inf; Y_1 = Inf; Y_2 = Inf; Z_1 = Inf; Z_2 = Inf

  if(dim(data_in)[3] != 1)
  {
  if(verbose_in) message("CALCULATING IN DIRECTION Z")
  tpb <- utils::txtProgressBar(min = 0, max = dim(data_in)[3], style = 3)

  logic_z = rep(0, dim(data_in)[3])

  for (i in 1:dim(data_in)[3])
  {
    flat = data_in[1:dim(data_in)[1], 1:dim(data_in)[2], i]
    ref = matrix(zero_value, dim(data_in)[1], dim(data_in)[2])

    if(max(abs(flat- ref)) != 0) {logic_z[i] = 1}

    if(verbose_in) utils::setTxtProgressBar(tpb, i); if(i == dim(data_in)[1]) {close(tpb)}
  }
  Z_1 <- suppressWarnings(min(which(logic_z %in% 1))); if(Z_1 == Inf | Z_1 == -Inf) {data_in_mod <- data_in[,,]}
  Z_2 <- suppressWarnings(max(which(logic_z %in% 1))); if(Z_2 == Inf | Z_2 == -Inf) {data_in_mod <- data_in[,,]}
  }

  X_1 <- suppressWarnings(min(which(logic_x %in% 1))); if(X_1 == Inf | X_1 == -Inf) {data_in_mod <- data_in[,,]}
  X_2 <- suppressWarnings(max(which(logic_x %in% 1))); if(X_2 == Inf | X_2 == -Inf) {data_in_mod <- data_in[,,]}
  Y_1 <- suppressWarnings(min(which(logic_y %in% 1))); if(Y_1 == Inf | Y_1 == -Inf) {data_in_mod <- data_in[,,]}
  Y_2 <- suppressWarnings(max(which(logic_y %in% 1))); if(Y_2 == Inf | Y_2 == -Inf) {data_in_mod <- data_in[,,]}

  if( (X_1 != Inf & X_1 != -Inf | X_2 != Inf & X_2 != -Inf | Y_1 != Inf & Y_1 != -Inf | Y_2 != Inf & Y_2 != -Inf) & dim(data_in)[3] == 1) {data_in_mod <- data_in[X_1:X_2, Y_1:Y_2, 1]
  }
  if( dim(data_in)[3] != 1 && (X_1 != Inf & X_1 != -Inf & X_2 != Inf & X_2 != -Inf & Y_1 != Inf & Y_1 != -Inf & Y_2 != Inf & Y_2 != -Inf & Z_1 != Inf & Z_1 != -Inf & Z_2 != Inf & Z_2 != -Inf) ) {data_in_mod <- data_in[X_1:X_2, Y_1:Y_2, Z_1:Z_2]
  }

  if(length(dim(data_in_mod)) == 2) {data_in_mod <- array(data_in_mod, dim = c(dim(data_in_mod),1))}

  if(verbose_in) { if(all(dim(data_in_mod) ==  dim(data_in))) {
    message("CROPPING WAS NOT DONE, SINCE NO DATA WAS OBSERVED IN THE DICOM DATASET, THEREFORE ORIGINAL DATA IS COPIED TO MODIF SLOT OF RIA_image");
  } else {message("CROPPING DOME IN ALL DIRECTIONS")}
  }

  if(any(class(RIA_data_in) == "RIA_image"))
  {
    if(write_orig){RIA_data_in$data$orig <- data_in_mod
    } else {RIA_data_in$data$modif<- data_in_mod}

    if(!any(class(RIA_data_in$data) == "RIA_data")) class(RIA_data_in$data) <- append(class(RIA_data_in$data), "RIA_data")

    RIA_data_in$log$logic_x <- logic_x; RIA_data_in$log$logic_y <- logic_y
    if((Z_1 != Inf & Z_1 != -Inf) & (Z_2 != Inf & Z_2 != -Inf)) {RIA_data_in$log$logic_z <- logic_z
    } else {RIA_data_in$log$logic_z <- NULL}
    RIA_data_in$log$events <- append(RIA_data_in$log$events, "Cropped")
    return(RIA_data_in)
  } else return(data_in_mod) #Non RIA_image currently not supported
}
