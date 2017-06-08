## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- results = "hide"---------------------------------------------------
library(RIA)

## ---- eval=FALSE---------------------------------------------------------
#  DICOM = load_dicom(filename = "C:\\DICOM\\")

## ---- eval=FALSE---------------------------------------------------------
#  DICOM = load_dicom(filename = "C:\\DICOM\\", crop_in = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  DICOM = load_dicom(filename = "C:\\DICOM\\", replace_in = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  DICOM = load_dicom(filename = "C:\\DICOM\\", center_in = FALSE)

## ------------------------------------------------------------------------
DICOM_codes

## ---- eval=FALSE---------------------------------------------------------
#  add <- as.data.frame(array(c("Manufacturer", "0008", "0070"), dim = c(1,3)))
#  colnames(add) <- c("Name", "Group", "Element")
#  DICOM = load_dicom(filename = "C:\\DICOM\\", header_add = add)

## ---- eval=FALSE---------------------------------------------------------
#  exclude <- as.data.frame(DICOM_codes[3:6,])
#  DICOM = load_dicom(filename = "C:\\DICOM\\", header_exclude = exclude)

## ------------------------------------------------------------------------
DICOM <- NRS

## ------------------------------------------------------------------------
DICOM = dichotomize(RIA_data_in = DICOM, bins_in = 2)

## ------------------------------------------------------------------------
DICOM = dichotomize(RIA_data_in = DICOM, bins_in = 2, equal_prob = TRUE)

## ------------------------------------------------------------------------
DICOM$log$events

## ------------------------------------------------------------------------
DICOM$log$cuts_es_2
DICOM$log$cuts_ep_2

## ------------------------------------------------------------------------
DICOM = dichotomize(RIA_data_in = DICOM, bins_in = c(4,8,16,32))

## ------------------------------------------------------------------------
DICOM = dichotomize(RIA_data_in = DICOM, bins_in = 2^(2:5), equal_prob = TRUE)

## ------------------------------------------------------------------------
DICOM$log$events
names(DICOM$dichotomized)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM)
RIA:::list_to_df(DICOM$stat_fo$orig)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE)
RIA:::list_to_df(DICOM$stat_fo$ep_32)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE, use_slot = "dichotomized$es_2")
RIA:::list_to_df(DICOM$stat_fo$es_2)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_orig = FALSE, use_slot = "dichotomized$es_2", save_name = c("equaly_sized_2bins"))
RIA:::list_to_df(DICOM$stat_fo$equaly_sized_2bins)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "dichotomized")
names(DICOM$stat_fo)

## ------------------------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "dichotomized",
                    save_name = c("Name_1", "Name_2", "Name_3", "Name_4",
                                  "Name_5", "Name_6", "Name_7", "Name_8",
                                  "Name_9", "Name_10"))
names(DICOM$stat_fo)

## ---- error = TRUE-------------------------------------------------------
DICOM = first_order(RIA_data_in = DICOM, use_type = "dichotomized",
                    save_name = c("Name_1", "Name_2", "Name_3", "Name_4"))

## ------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, off_right = 1, off_down = 2, off_z = 2)
dim(DICOM$glcm$ep_32)

## ------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_slot = "dichotomized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = FALSE, normalize = FALSE)
DICOM$glcm$ep_2
DICOM = glcm(RIA_data_in = DICOM, use_slot = "dichotomized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = FALSE)
DICOM$glcm$ep_2

## ------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_slot = "dichotomized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = FALSE)
DICOM$glcm$ep_2
DICOM = glcm(RIA_data_in = DICOM, use_slot = "dichotomized$ep_2", off_right = 1, off_down = 2, off_z = 2, symmetric = TRUE, normalize = TRUE)
DICOM$glcm$ep_2

## ------------------------------------------------------------------------
DICOM = glcm(RIA_data_in = DICOM, use_type = "dichotomized",
             off_right = 1, off_down = 2, off_z = 2)
names(DICOM$glcm)

## ------------------------------------------------------------------------
DICOM = glcm_stat(DICOM, use_type = "single", use_slot = "glcm$ep_2_122")

## ------------------------------------------------------------------------
DICOM = glcm_stat(DICOM, use_type = "glcm")
names(DICOM$stat_glcm)

## ------------------------------------------------------------------------
DICOM = glrlm(RIA_data_in = DICOM, right = 1, down = 0, forward = 1)
dim(DICOM$glrlm$ep_32)

## ------------------------------------------------------------------------
DICOM = glrlm(RIA_data_in = DICOM, use_type = "dichotomized",
             right = 1, down = 0, forward = 1)
names(DICOM$glrlm)

## ------------------------------------------------------------------------
DICOM = glrlm_stat(DICOM, use_type = "single", use_slot = "glrlm$ep_2_101")

## ------------------------------------------------------------------------
DICOM = glrlm_stat(DICOM, use_type = "glrlm")
names(DICOM$stat_glrlm)

## ------------------------------------------------------------------------
DICOM = geometry(RIA_data_in = DICOM, use_orig = TRUE, calc_sub = FALSE)
RIA:::list_to_df(DICOM$stat_geometry$orig)

## ------------------------------------------------------------------------
DICOM = geometry(RIA_data_in = DICOM, use_slot = "dichotomized$es_2")
RIA:::list_to_df(DICOM$stat_geometry$es_2)

## ------------------------------------------------------------------------
DICOM = geometry(DICOM, use_type = "dichotomized")
names(DICOM$stat_geometry)

