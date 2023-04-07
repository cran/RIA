# Retrieves NIFTI header information from images in directory
#
# Returns a RIA_header object containing the NIFTI header information as a list
# (c): Márton Kolossváry, 2018


create_header_npy <- function(directory, PixelSpacing, SpacingBetweenSlices)
{
  header_list <- list()

  header_list[["PixelSpacing"]] <- PixelSpacing
  header_list[["SpacingBetweenSlices"]] <- SpacingBetweenSlices
  
  if(!any(class(header_list) == "RIA_header")) class(header_list) <- append(class(header_list), "RIA_header")
  
  return(header_list)
}
