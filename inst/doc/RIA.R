## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- results = "hide"--------------------------------------------------------
library(RIA)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM = load_dicom(filename = "C:/DICOM/")

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM = load_dicom(filename = "C:/DICOM/", crop_in = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM = load_dicom(filename = "C:/DICOM/", replace_in = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM = load_dicom(filename = "C:/DICOM/", center_in = FALSE)

## -----------------------------------------------------------------------------
DICOM_codes

## ---- eval=FALSE--------------------------------------------------------------
#  add <- as.data.frame(array(c("Manufacturer", "0008", "0070"), dim = c(1,3)))
#  colnames(add) <- c("Name", "Group", "Element")
#  DICOM = load_dicom(filename = "C:/DICOM/", header_add = add)

## ---- eval=FALSE--------------------------------------------------------------
#  exclude <- as.data.frame(DICOM_codes[3:6,])
#  DICOM = load_dicom(filename = "C:/DICOM/", header_exclude = exclude)

## -----------------------------------------------------------------------------
DICOM <- NRS

## -----------------------------------------------------------------------------
DICOM = discretize(RIA_data_in = DICOM, bins_in = 2)

## -----------------------------------------------------------------------------
DICOM = discretize(RIA_data_in = DICOM, bins_in = 2, equal_prob = TRUE)

## -----------------------------------------------------------------------------
DICOM$log$events

## -----------------------------------------------------------------------------
DICOM$log$cuts_es_2
DICOM$log$cuts_ep_2

## -----------------------------------------------------------------------------
DICOM = discretize(RIA_data_in = DICOM, bins_in = c(4,8,16,32))

## -----------------------------------------------------------------------------
DICOM = discretize(RIA_data_in = DICOM, bins_in = 2^(2:5), equal_prob = TRUE)

## -----------------------------------------------------------------------------
DICOM$log$events
names(DICOM$discretized)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM)
RIA:::list_to_df(DICOM$stat_fo$orig)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE)
RIA:::list_to_df(DICOM$stat_fo$ep_32)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE, use_slot = "discretized$es_2")
RIA:::list_to_df(DICOM$stat_fo$es_2)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE, use_slot = "discretized$es_2", save_name = c("equaly_sized_2bins"))
RIA:::list_to_df(DICOM$stat_fo$equaly_sized_2bins)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "discretized")
names(DICOM$stat_fo)

## -----------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "discretized",
                    save_name = c("Name_1", "Name_2", "Name_3", "Name_4",
                                  "Name_5", "Name_6", "Name_7", "Name_8",
                                  "Name_9", "Name_10"))
names(DICOM$stat_fo)

## ---- error = TRUE------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "discretized",
                    save_name = c("Name_1", "Name_2", "Name_3", "Name_4"))

## -----------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, off_right = 1, off_down = 2, off_z = 2)
dim(DICOM$glcm$ep_32)

## -----------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_slot = "discretized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = FALSE, normalize = FALSE)
DICOM$glcm$ep_2
DICOM = glcm(RIA_data_in = DICOM, use_slot = "discretized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = FALSE)
DICOM$glcm$ep_2

## -----------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_slot = "discretized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = FALSE)
DICOM$glcm$ep_2
DICOM = glcm(RIA_data_in = DICOM, use_slot = "discretized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = TRUE)
DICOM$glcm$ep_2

## -----------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_type = "discretized",
             off_right = 1, off_down = 2, off_z = 2)
names(DICOM$glcm)

## -----------------------------------------------------------------------------
DICOM = glcm_all(DICOM, use_type = "single")
names(DICOM$glcm)

## -----------------------------------------------------------------------------
DICOM = glcm_all(DICOM, use_type = "discretized")
names(DICOM$glcm)

## -----------------------------------------------------------------------------
DICOM = glcm_all(DICOM, use_type = "discretized", distance = c(1:3))
names(DICOM$glcm)

## -----------------------------------------------------------------------------
DICOM = glcm_stat(DICOM, use_type = "single", use_slot = "glcm$ep_2_122")

## -----------------------------------------------------------------------------
DICOM = glcm_stat(DICOM, use_type = "glcm")
names(DICOM$stat_glcm)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM <- glcm_stat_all(DICOM)

## -----------------------------------------------------------------------------
DICOM = glrlm(RIA_data_in = DICOM, off_right = 1, off_down = 0, off_z = 1)
dim(DICOM$glrlm$ep_32)

## -----------------------------------------------------------------------------
DICOM = glrlm(RIA_data_in = DICOM, use_type = "discretized",
              off_right = 1, off_down = 0, off_z = 1)
names(DICOM$glrlm)

## -----------------------------------------------------------------------------
DICOM = glrlm_all(DICOM, use_type = "single")
names(DICOM$glrlm)

## -----------------------------------------------------------------------------
DICOM = glrlm_all(DICOM, use_type = "discretized")
names(DICOM$glrlm)

## -----------------------------------------------------------------------------
DICOM = glrlm_stat(DICOM, use_type = "single", use_slot = "glrlm$ep_2_101")

## -----------------------------------------------------------------------------
DICOM = glrlm_stat(DICOM, use_type = "glrlm")
names(DICOM$stat_glrlm)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM <- glrlm_stat_all(DICOM)

## -----------------------------------------------------------------------------
DICOM = geometry(RIA_data_in = DICOM, use_orig = TRUE, calc_sub = FALSE)
RIA:::list_to_df(DICOM$stat_geometry$orig)

## -----------------------------------------------------------------------------
DICOM = geometry(RIA_data_in = DICOM, use_slot = "discretized$es_2")
RIA:::list_to_df(DICOM$stat_geometry$es_2)

## -----------------------------------------------------------------------------
DICOM = geometry(DICOM, use_type = "discretized")
names(DICOM$stat_geometry)

## ---- eval=FALSE--------------------------------------------------------------
#  DICOM <- radiomics_all(DICOM, bins_in = c(8, 16, 32), equal_prob = "both")

## ---- eval=FALSE--------------------------------------------------------------
#  save_RIA(DICOM, save_to = "C:/Test/", save_name = "My_first_radiomics", group_name = "Case")

## ---- eval=FALSE--------------------------------------------------------------
#  folder <- "/Images/" #Location of folder containing individual folders per patient which contain nrrd files for the image and mask
#  out <- "/CSV/" #Location of folder where the results will be dumped
#  
#  patients <- list.dirs(folder, recursive = FALSE, full.names = FALSE) #Name of patient folders
#  patients_full <- list.dirs(folder, recursive = FALSE, full.names = TRUE) #Name of patient folders with full file path name
#  
#  library(foreach); library(doParallel) #Load required packages
#  doParallel::registerDoParallel(7) #Define how many threads to use, usually use the number of threads-1
#  
#  #Use parallelized for cycle to cycle through all the patients
#  data_out_paral <- foreach (i = 1:length(patients), .combine="rbind", .inorder=FALSE,
#                             .packages=c('RIA'), .errorhandling = c("pass"), .verbose=FALSE) %dopar% {
#  
#                               files <- list.files(patients_full[i]) #Names of the files in the current patient folder
#                               image <- grep("heart", files, ignore.case = T, value = T) #Full name of the image file
#                               masks <- grep("plaque", files, ignore.case = T, value = T) #Full name of the mask files
#  
#                               #RUN RIA
#                               IMAGE <- RIA::load_nrrd(filename = paste0(patients_full[i], "/", image),
#                                                       mask_filename = paste0(patients_full[i], "/", masks), switch_z = FALSE) #Load image and mask files
#                               IMAGE <- RIA::radiomics_all(IMAGE,  equal_prob = "both", bins_in= c(4, 8, 16)) #Calculate radiomic features
#                               RIA::save_RIA(IMAGE, save_to = out, save_name = patients[i], group_name = patients[i]) #Export results into csv
#                             }

