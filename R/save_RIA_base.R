#Helpter function for save_RIA
# (c): Márton Kolossváry, 2018

save_RIA_base <- function(RIA_image, RIA_name) {
  
  names_stat <- names(eval(parse(text = paste0("RIA_image$", RIA_name))))
  df_out <- NULL
  
  for (i in 1: length(names_stat)) {
    if(RIA_name == "stat_geometry") {
      
      names_stat_geo <- names(eval(parse(text = paste0("RIA_image$", RIA_name, "$", names_stat[i]))))
      
      for (j in 1: length(names_stat_geo)) {
        name_data_geo <- t(list_to_df(eval(parse(text = paste0("RIA_image$",RIA_name, "$", names_stat[i], "$", names_stat_geo[j])))))
        colnames(name_data_geo) <- paste0(names_stat_geo[j], "_", colnames(name_data_geo), "__", names_stat[i])
        df_out <- cbind(df_out, name_data_geo)
      }
      
      
    } else {
      name_data <- t(list_to_df(eval(parse(text = paste0("RIA_image$",RIA_name, "$", names_stat[i])))))
      colnames(name_data) <- paste0(colnames(name_data), "__", names_stat[i])
      df_out <- cbind(df_out, name_data)
    }
    
  }
  return(df_out)
}
