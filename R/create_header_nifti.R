# Retrieves NIFTI header information from images in directory
#
# Returns a RIA_header object containing the NIFTI header information as a list
# (c): Márton Kolossváry, 2018


create_header_nifti <- function(directory)
{
  
  dcmHeader <- oro.nifti::nifti_header(directory, verbose = FALSE)
  dcm_capture <- utils::capture.output(dcmHeader)
  
  header_list <- list()
  for(i in 2:length(dcm_capture)) {
    
    space     <- gregexpr('  ', dcm_capture[i])[[1]][2]
    if(is.na(space)) {
      space     <- gregexpr(':', dcm_capture[i])[[1]][1]-1
    }
    var  <- substr(dcm_capture[i], 3, space-1)
    colon     <- gregexpr(':', dcm_capture[i])[[1]][1]
    length_var <- nchar(dcm_capture[i])
    val <- substr(dcm_capture[i], colon+2, length_var)
    
    header_list[[var]] <- val
  }
  
  PixelSpacing <- as.numeric(substr( header_list$`Pixel Dimension`, 1, regexpr(' ', header_list$`Pixel Dimension`)[1]-1))
  if(length(gregexpr('x', header_list$`Pixel Dimension`)[[1]]) > 2) {
    SpacingBetweenSlices <- as.numeric(substr( header_list$`Pixel Dimension`, gregexpr('x', header_list$`Pixel Dimension`)[[1]][2]+2,
                                               gregexpr('x', header_list$`Pixel Dimension`)[[1]][3]-2))
  } else {
    SpacingBetweenSlices <- as.numeric(substr( header_list$`Pixel Dimension`, gregexpr('x', header_list$`Pixel Dimension`)[[1]][2]+2,
                                               nchar(header_list$`Pixel Dimension`)))
  }

  header_list[["PixelSpacing"]] <- PixelSpacing
  header_list[["SpacingBetweenSlices"]] <- SpacingBetweenSlices
  
  if(!any(class(header_list) == "RIA_header")) class(header_list) <- append(class(header_list), "RIA_header")
  
  return(header_list)
}
