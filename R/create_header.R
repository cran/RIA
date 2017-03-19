# Retrieves DICOM header information from images in directory
#
# @description  Creates RIA_header object given the directory containing the DICOM files
# based on predefined DICOM code set specified in DICOM_codes.rda. By editing
# DICOM_codes.rda, one can modify the preset default values to be exported.
#
# directory filepath to DICOM directory
#
# addition dataframe with three columns: Name Group Element containing the name,
# the group and the element code of the DICOM fields wished to be added to the predefined set
#
# exclusion dataframe with three columns: Name Group Element containing the name,
# the group and the element code of the DICOM fields wished to be excluded from the predefined set
#
# reuturn Returns a RIA_header object contining the DICOM header information as a list
#
# examples
# Loads predefined variables specified by DICOM_codes.rda from the DICOM header of the image
# header <- create_header("C://Users//Test//Documents//Radiomics//John_Smith//DICOM_folder//")
#
# Loads predefined variables specified by DICOM_codes.rda except patient identifying information
# exclude <- as.data.frame(DICOM_codes[3:6,])
# header <- create_header("C://Users//Test//Documents//Radiomics//John_Smith//DICOM_folder//",
# exclusion = exclude)
#
# Loads predefined set of variables specified by DICOM_codes.rda except patient information
# but adds Manufacturer
# add <- as.data.frame(array(c("Manufacturer", "0008", "0070"), dim = c(1,3)))
# colnames(add) <- c("Name", "Group", "Element")
# exclude <- as.data.frame(DICOM_codes[3:6,])
# header <- create_header("C://Users//Test//Documents//Radiomics//John_Smith//DICOM_folder//",
# addition = add, exclusion = exclude)

create_header <- function(directory, addition = NULL, exclusion = NULL)
{

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
