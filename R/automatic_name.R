# General function to automatically assign subslot name
# (c): Márton Kolossváry, 2018

automatic_name <- function(RIA_data_in, orig_in, use_slot)
{
  if(!is.null(use_slot)) {
    loc_dollar <- base::gregexpr(pattern = "$", text = use_slot, fixed = TRUE)
    txt <- base::substr(use_slot, max(loc_dollar[[1]])+1, nchar(use_slot))
    txt <- gsub("`", "", txt)
    return(txt)
  }
  
  if(orig_in) {txt <- "orig"; return(txt)}
  last_event <- RIA_data_in$log$events[length(RIA_data_in$log$events)]
  
  if(length(grep("equal_sized", last_event)) >0) {save_txt <- "es_"
  } else if(length(grep("equal_prob", last_event)) >0) {save_txt <- "ep_"
  } else {save_txt <- "orig"}
  
  num_ind <- unlist(gregexpr('[1-9]', last_event))
  num_txt <- substr(last_event, num_ind[1], num_ind[length(num_ind)])
  
  if(save_txt ==  "orig") txt <- save_txt
  else txt <- paste0(save_txt, num_txt)
  
  return(txt)
}
