# Retrieves DICOM header information from DICOM images in directory
#
# addition dataframe with three columns: Name Group Element containing the name,
# the group and the element code of the DICOM fields wished to be added to the predefined set
# exclusion dataframe with three columns: Name Group Element containing the name,
# the group and the element code of the DICOM fields wished to be excluded from the predefined set.
# Returns a RIA_header object containing the DICOM header information as a list
# (c): Márton Kolossváry, 2018

create_header <- function(directory, addition = NULL, exclusion = NULL)
{
  
  if(substr(directory, nchar(directory), nchar(directory)) != "/") directory <- paste0(directory, "/")
  
  singe_dicom <- paste(directory, list.files(directory)[1], sep = "")
  dcmImages<- oro.dicom::readDICOMFile(singe_dicom)
  header_path <- paste("dcmImages$hdr", sep = "")
  
  header <- as.data.frame(eval(parse(text = header_path)))
  
  DICOM_codes_df <- as.data.frame(RIA::DICOM_codes)
  if (!is.null(addition)) DICOM_codes_df <- rbind(DICOM_codes_df, addition)
  if (!is.null(exclusion)) DICOM_codes_df <- DICOM_codes_df[((!DICOM_codes_df$Group %in% exclusion$Group) | (!DICOM_codes_df$Element %in% exclusion$Element)),]
  header_df <- data.frame(matrix(NA, nrow = length(DICOM_codes_df$Name), ncol = 2))
  
  options(warn=-1)
  for (i in 1: dim(DICOM_codes_df)[1])
  {
    one_row <- subset(header, as.numeric(DICOM_codes_df$Group[i]) == as.numeric(header$group) & as.numeric(DICOM_codes_df$Element[i]) == as.numeric(header$element))[1,]
    header_df[i,2] <- one_row$value
  }
  options(warn=0)
  
  header_l <- stats::setNames(split(header_df[,2], seq(nrow(header_df))), DICOM_codes_df$Name)
  if(!any(class(header_l) == "RIA_header")) class(header_l) <- append(class(header_l), "RIA_header")
  
  return(header_l)
}
