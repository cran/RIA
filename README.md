<!-- README.md is generated from README.Rmd. Please edit that file -->

Radiomics Image Analysis (RIA)
==============================

Radiomics Image Analysis (RIA) package was developed to facilitate
radiomic analysis of medical images. `RIA` can calculate hundreds of
different statistics on both 2D and 3D images. `RIA` supports analysis
of `DICOM`, `NIfTI` and `nrrd` images. Almost all calculations are
vectorized and therefore are super-efficient. The package is developed
by Márton Kolossváry a medical doctor not an engineer, therefore all
functionalities of the software package are developed in a way that can
be learnt by non-professionals. `RIA` is constantly updated with new
functionalities and wrap-around functions to make the calculation of
radiomic metrics even simpler.

Its as easy as 1, 2, 3
======================

RIA allows users to take control of each and every aspect of radiomic
analysis using specific functions. However, for most users 3 lines of
simple code: one loading the data and one calculating the statistics,
and one exporting the results is enough:

``` r
#Load the data by providing the location of the DICOM, NIfTI or nrrd file(s)
DICOM <- load_dicom(filename = "C:/Image/")

#Calculate first-order, GLCM, GLRLM and geometry based statistics
DICOM <- radiomics_all(DICOM, equal_prob = "both")

#Save output to csv
save_RIA(DICOM, save_to = "C:/Test/", save_name = "My_first_radiomics", group_name = "Case")
```

These three simple lines of code result in thousands of radiomic
parameters calculated for the given image! For a more detailed
introduction to RIA please read the
[vignette](https://CRAN.R-project.org/package=RIA/vignettes/RIA.html).
If you wish to better understand Radiomics I would suggest reading
[“Cardiac Computed Tomography Radiomics: A Comprehensive Review on
Radiomic Techniques”](https://pubmed.ncbi.nlm.nih.gov/28346329/) and
[“Radiomic Features Are Superior to Conventional Quantitative Computed
Tomographic Metrics to Identify Coronary Plaques With Napkin-Ring
Sign”](https://pubmed.ncbi.nlm.nih.gov/29233836/) which describes the
calculation and each statistic in detail in the supplementary files.
