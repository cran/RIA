<!-- README.md is generated from README.Rmd. Please edit that file -->
Radiomics Image Analysis (RIA)
==============================

Radiomics Image Analysis (RIA) package was developed to facilitate radiomic analysis of grayscale images. RIA can calculate hundreds of different statistics on both 2D and 3D images. Almost all calculations are vectorized and therefore are super-efficient. The package is developed by Márton Kolossváry a medical doctor not an engineer, therefore all functionalities of the software package are developed in a way that can be learnt by non-professionals. RIA is constantly updated with new functionalities and wrap-around functions to make the calculation of radiomic metrics even simpler.

Its as easy as 1, 2 ...
=======================

RIA allows users to take control of each and every aspect of radiomic analysis using specific functions. However, for most users 2 lines of simple code: one loading the data and one calculating the statistics, is enough:

``` r
#Load the data by providing the location of the DICOM files
DICOM <- load_dicom(filename = "C:/DICOM/")
#Calculate first-order, GLCM, GLRLM and geometry based statistics
DICOM <- radiomics_all(DICOM, equal_prob = FALSE, bins_in = c(2,8,32,128), distance = c(1:3), fo_discretized = TRUE, geometry_discretized = TRUE)
```

These two simple lines of code result in thousands of radiomic parameters calculated for the given image! For a more detailed introduction to RIA please read the [vignette](https://CRAN.R-project.org/package=RIA/vignettes/RIA.html). If you wish to better understand Radiomics I would suggest reading ["Cardiac Computed Tomography Radiomics: A Comprehensive Review on Radiomic Techniques"](https://www.ncbi.nlm.nih.gov/pubmed/28346329) and ["Radiomic Features Are Superior to Conventional Quantitative Computed Tomographic Metrics to Identify Coronary Plaques With Napkin-Ring Sign"](http://circimaging.ahajournals.org/content/10/12/e006843) which describes the calculation and each statistic in detail in the supplementary files.
