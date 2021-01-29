# Changes list to dataframe
#
# Takes in list. List elements can have multiple numeric elements
#
# list_in list to be transformed to dataframe
# sci_not whether to present results with scientific notation
#
# return dataframe of values
# (c): Márton Kolossváry, 2018


list_to_df <- function(list_in, sci_not = FALSE)
{
  if(!sci_not) options(scipen=999)
  df <- do.call(rbind, lapply(list_in, data.frame, stringsAsFactors=FALSE))
  colnames(df) <- "Values"
  return(df)
}
