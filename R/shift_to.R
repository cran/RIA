# Shift image intensities by given amount
#
# Shifts image intensities by given amount. The amount of shift will be determined so that the
# smallest original value indicated by "min_value_in will" have the value "to" after transformation.
# Used by load_dicom to shift values if necessary.
#
# RIA_data_in RIA_image, matrix or array
#
# to if min_value_in provided then the amount will be determined so that the smallest value will equal
# to the to parameter. If not, then all values will be added to
#
# min_value_in integer provided indicating smallest value
#
# use_orig indicating to use image present in RIA_data$orig. If FALSE, the modified image
# will be used stored in RIA_data$modif.
#
# write_orig indicating to write cropped image  to RIA_data$orig. If FALSE, the modified image
# will be used stored in RIA_data$modif.
#
# verbose_in logical indicating whether to print detailed information.
# Most prints can also be suppresed using the suppressMessages function.
#
# Return: RIA_image object croped to smallest bouding box, whith RIA_log updated.



shift_to <- function(RIA_data_in, to = -1024, min_value_in = NULL, use_orig = TRUE, write_orig = TRUE, verbose_in = TRUE)
{
  data_in <- check_data_in(RIA_data_in, use_type = "single", use_orig = use_orig, verbose_in = verbose_in)


  if (!is.null(min_value_in)) {shift <- (-min_value_in + to)
  }  else {shift <- to}
  RIA_data_mod <- data_in + shift


  if(any(class(RIA_data_in) == "RIA_image") )
  {
    if(write_orig) {RIA_data_in$data$orig <- RIA_data_mod
    } else {RIA_data_in$data$modif<- RIA_data_mod}

    if(!any(class(RIA_data_in$data) == "RIA_data")) class(RIA_data_in$data) <- append(class(RIA_data_in$data), "RIA_data")

    RIA_data_in$log$shift <- shift
    RIA_data_in$log$events <- append(RIA_data_in$log$events, paste("Shifted_by_", shift, sep = ""))
    return(RIA_data_in)
  }
  else return(RIA_data_mod)
}
