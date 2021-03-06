ProteoModlR 
=========

R Suite for Quantitative Proteomics Pathway Modeling.

Overview of pipeline
--------------------

ProteoModlR is an open-sourced R suite for quantitative mass spectrometry analysis of the relative concentration of proteins and the stoichiometry of post-translational chemical modifications. Due to its modular design and flexible analysis pipeline, ProteoModlR allows for seamless integration with existing proteomics software, such as MaxQuant and Skyline, as well as with statistical and pathway analysis tools. It facilitates analysis and visualization of quantitative proteomics data enabling researchers to assess differential activation of functional cellular processes. 

Manifest
--------

QC.R - Script for format verification and data classification

Normalize.R - Script for data normalization

Analyze.R - Script for calculating abundance and stoichiometry

Documentation
--------
Please refer to the PDF file in the vignettes folder

Installation
------------

Enter the following commands directly into R console (assumes devtools is installed):
    
    library(devtools)
    install_github('kentsisresearchgroup/ProteoModlR')

Dependencies
------------

'reshape2', 'plyr, and 'ggplot2' from CRAN

Authors
-------

* Mojdeh Shakiba | mojdeh.shakiba24@gmail.com
* Paolo Cifani | cifanip@mskcc.org
* Alex Kentsis | kentsisresearchgroup@gmail.com
