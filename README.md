# Chromatin_Modification_scripts
This repository contains the code used to perform the analysis, data processing and models as descibed in the project 'Multi-omic characterization of novel histone post-translational modifications', conducted by Agustín Vera Enguídanos and supervised by Gonzalo Millán Zambrano.

The main data processing pipeline included in this project is the AVE_chipseq.sh shell script. This pipelines requires as an input a single .txt parameters  file where the input files and pipelines specifications must be indicated. Please see the example file. The AVE_chipseq.sh main script can call to other auxilliary scripts as it runs, AVE_read_mapping.sh and AVE_peak_calling.sh. AVE_ChromHMM.sh can be runned to generate a chromatin segmentation model.

Once the ChIP-seq reads have been mapped into BAM files and the coverage calculated into bigWig files, other shell scripts can be used for data representation, mainly AVE_bigwig_processing, AVE_metaplot or AVE_library_complexity. Selected commands can be copy/pasted for a personalizad analysis. 

This repository also includes several R scripts that have been used for further data analysis and representation of the counts data and the chromHMM model.

Please note that not all of the scripts and functions included in this repository have been used in the final manuscript of the project. Instead, this repository intends to be useful for a variety of applications for ChIP-seq data analysis, both within the GMZ Chromatin Modifications group and for external users.
