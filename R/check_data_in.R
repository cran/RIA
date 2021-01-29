# General function to check validitiy of RIA_image object
# (c): Márton Kolossváry, 2018

check_data_in <- function(RIA_data_in, use_type = "single", use_orig = TRUE, use_slot = NULL, verbose_in = TRUE)
{
    if(!any(class(RIA_data_in) == "RIA_image")) {message("PROCESSING OF RIA_image OBJECTS ARE SUPPORTED, OTHER CLASSES MIGHT CAUSE PROBLEMS! PLEASE LOAD DATA USING load_dicom")}
    
    if(use_type == "single") {
        if(!is.null(use_slot)) {data_in <- eval(parse(text = paste0("RIA_data_in$", use_slot)))
        } else if(any(class(RIA_data_in) == "RIA_image") & use_orig) {data_in <- RIA_data_in$data$orig
        } else if(any(class(RIA_data_in) == "RIA_image") & !use_orig) {data_in <- RIA_data_in$data$modif
        } else { #non RIA_image class
            if(!(any(class(RIA_data_in) == "matrix") | any(class(RIA_data_in) == "array"))) {
                if(verbose_in) {message(" "); message(paste("INPUT IS NOT MATRIX OR ARRAY, WILL USE as.matrix FUNCTION. PLEASE CHECK IF THIS CAUSES ANY PROBLEMS!", sep = ''))
                    data_in <- as.matrix(RIA_data_in)}}
            else {data_in <- RIA_data_in}
        }
        
        if(length(dim(data_in)) < 2 | length(dim(data_in)) > 3) stop(paste0("DATA LOADED IS ", length(dim(data_in)), " DIMENSIONAL. ONLY 2D AND 3D DATA ARE SUPPORTED!"))
        
        data_NA <- as.vector(data_in)
        data_NA <- data_NA[!is.na(data_NA)]
        
        if(length(data_NA) == 0) {message(" "); message("WARNING: RIA_image$data DOES NOT CONTAIN ANY DATA!!!"); message(" ")}
        
        
        return(data_in)
        
    }
    
    if(use_type == "discretized") {
        return(RIA_data_in$discretized)
    }
    
    if(use_type == "glcm") {
        return(RIA_data_in$glcm)
    }
    
    if(use_type == "glrlm") {
        return(RIA_data_in$glrlm)
    }
    
}