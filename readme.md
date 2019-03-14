Supplementary material of the paper "Blood is thicker than bloodshed"
================
Diego Alburez-Gutierrez

About this repository
---------------------

This repository provides the code needed to reproduce the results of the study:

Alburez-Gutierrez, Diego (2019) Blood is thicker than bloodshed: A genealogical approach to reconstruct populations after armed conflicts, *Demographic Research* 40(23): 627-656. doi: 10.4054/DemRes.2019.40.23.

The repository contains the source code of the `EGM` R package (EGM stands for Extended Genealogy Method). The package contains the R functions used to analyse the data collected by the author in the village of Rio Negro (2015-2016). It can be downloaded from RStudio using the `devtools` package in R:

``` r
  library(devtools)
  install_github("alburezg/EGM", dep = FALSE)
  library(EGM)    
```

The directory [1_reproducible_results](1_reproducible_results) of this repository includes the full code needed to reproduce the results, figures, and tables presented in the paper. 
The data used in this study is not included in the `EGM` package or in this repository in order to protect the privacy of the respondents.
However, this directory gives details about the structure of the data to the degree needed to understand how the analysis was conducted. 

The directory [2_appendix_a](2_appendix_a) gives more detail about how the participants for the study were selected and the genealogical data managed. This document can be built in RStudio and requires the `EGM` package described above.

The directory [3_questionnaires](3_questionnaires) contains the questionnaires used to collect the genealogical data in Rio Negro. The translated questionnaires have been edited for style and can be downloaded as pdf files.

Get involved!
-------------

Researchers interested in using the data for academic analysis or learning more about the Extended Genealogy Method (EGM) can contact Diego Alburez at: alburezgutierrez\[at\]demogr.mpg.de.