# Change given value of image to a specified value
#
# Change given value of \emph{RIA_image} to a given value keeping the dimentions of the object. Used by \code{\link{load_dicom}} to change smallest values to NA, which are considered to indicate values without a signal.
#
# RIA_data_in RIA_image
#
# zero_value_in the value that will be changed
#
# value_to the value to change to
#
# use_orig indicating to use image present in RIA_data$orig.
# If FALSE, the modified image will be used stored in RIA_data$modif.
#
# write_orig indicating to write cropped image  to RIA_data$orig.
# If FALSE, the modified image will be used stored in RIA_data$modif.
#
# verbose_in logical indicating whether to print detailed information.
# Most prints can also be suppresed using the suppressMessages function.
#
# return RIA_image object where given values were changed to NA, whith RIA_log updated.

change_to<- function(RIA_data_in, zero_value_in = 0, value_to = NA, use_orig = TRUE, write_orig = TRUE, verbose_in = TRUE)
{
  data_in <- check_data_in(RIA_data_in, use_type = "single", use_orig = use_orig, verbose_in = verbose_in)

  data_in[data_in == zero_value_in] <- value_to


  if(any(class(RIA_data_in) == "RIA_image"))
  {
    if(write_orig) {RIA_data_in$data$orig <- data_in
    } else {RIA_data_in$data$modif<- data_in}

    if(!any(class(RIA_data_in$data) == "RIA_data")) class(RIA_data_in$data) <- append(class(RIA_data_in$data), "RIA_data")

    RIA_data_in$log$zero_value <- zero_value_in
    RIA_data_in$log$changed_to <- value_to
    RIA_data_in$log$events <- append(RIA_data_in$log$events, "Changed_to_NA")
    return(RIA_data_in)
  }
  else return(data_in)
}
