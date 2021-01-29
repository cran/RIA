# Retrieves nrrd header information from images in directory
#
# reuturn Returns a RIA_header object contining the NRRD header information as a list
# (c): Márton Kolossváry, 2018


create_header_nrrd <- function(directory)
{
  
  dcmHeader <- nat::read.nrrd(directory, Verbose = FALSE, ReadData = FALSE, AttachFullHeader = TRUE)
  dcm_capture <- attributes(dcmHeader)
  dcm_capture <- dcm_capture$header
  dcm_capture <- attributes(dcm_capture)
  var_names <- dcm_capture$names
  
  
  header_list <- list()
  header_list[["case"]] <- dcm_capture$headertext[1]
  header_list[["file"]] <- dcm_capture$path
  for(i in 4:(length(var_names)+3)) {
    
    colon     <- gregexpr(':', dcm_capture$headertext[i])[[1]][1]
    length_var <- nchar(dcm_capture$headertext[i])
    val <- substr(dcm_capture$headertext[i], colon+2, length_var)
    
    header_list[[var_names[i-3]]] <- val
  }
  
  dimentions <- nat::nrrd.voxdims(directory)
  
  PixelSpacing <- mean(dimentions[1], dimentions[2], na.rm = T)
  SpacingBetweenSlices <- dimentions[3]
  header_list[["PixelSpacing"]] <- PixelSpacing
  header_list[["SpacingBetweenSlices"]] <- SpacingBetweenSlices
  
  if(!any(class(header_list) == "RIA_header")) class(header_list) <- append(class(header_list), "RIA_header")
  
  return(header_list)
}
