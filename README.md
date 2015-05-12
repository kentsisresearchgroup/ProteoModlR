ProteoModlR 
=========


R Suite for Quantitative Proteomics Pathway Modeling.

Online documentation
--------------------
Pending.

Authors
-------

* Mojdeh Shakiba | mojdeh.shakiba24@gmail.com
* Paolo Cifani | cifanip@mskcc.org
* Alex Kentsis | kentsisresearchgroup@gmail.com

Overview of pipeline
--------------------

ProteoModlR is a R suite for analysis of quantitative mass spectrometry proteomics data, specifically designed to investigate variation in expression and post-translational modifications in biological pathway-centric manner. The program requires a tabular format input file that includes unique protein identifiers, corresponding peptide sequences with annotated chemical modifications, and abundance values, such as extraction ion currents and chromatographic peak alignments. Biologic pathway annotation is specified by an accompanying file that includes either user-specified pathways of interest, or publically available annotations such as the NCI Pathway Interaction Database. After basic quality control and normalization (if applicable), the data is analyzed for differential abundance and site occupancy of the various modifications, providing functional insight into pathway activity and regulation. 

Manifest
--------

QC.R - Script for format verification and parameter initialization. 

Kinetics.R - main code; dependent on QC.R.

Stats.R - statistical analysis; dependent on MSstats.R

Installation
------------

From source:

    git clone https://github.com/kentsisresearchgroup/ProteoModlR.git

Dependencies
------------

* R/Bioconductor version 3.0.2 - source("http://bioconductor.org/biocLite.R")
* MSstats - http://msstats.org/


